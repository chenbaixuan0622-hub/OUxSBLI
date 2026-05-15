# Plan: Fix LES Viscous Kernels in calc_visc2_curv.f90

## Context

The NS curvilinear viscous kernels (`calc_Ev2_curv`, `calc_Fv2_curv`, `calc_Gv2_curv`) have been verified correct and the solver runs with `id_visc = integer(4)`. However, deep audit of the three LES kernels (`calc_Ev_LES2_curv`, `calc_Fv_LES2_curv`, `calc_Gv_LES2_curv`) against the Cartesian reference (`calc_visc2.f90`) and the NS curvilinear kernels reveals systematic bugs: SGS gradients are not projected to physical (x,y) space via the curvilinear chain rule, and some stencil cross-references use the wrong direction.

Note: `Q(:,2,...) = u` (primitive velocity), not ρu — `calc_quantities_3D` converts before these kernels are called. The load_smem routines are correct.

---

## Bugs (all in LES kernels, id_visc=8 not yet wired)

### Bug 1 — `calc_Ev_LES2_curv` line 411: `mwxsgs` uses `mwzsgs` (z-stencil for w) for the η-contribution

The physical μ_sgs*∂w/∂x = μ_sgs*(∂w/∂ξ·ξ_x + ∂w/∂η·η_x). The ξ term is `mut_f*dw_dxi*xi_x_f`; the η term needs `mwysgs` (η-stencil for w with SGS μ), which is **never computed**. Currently `mwzsgs` (z-stencil for w) is used in its place.

### Bug 2 — `calc_Ev_LES2_curv` lines 425-427: SGS stress projection is incomplete

The molecular part correctly does `txx = (txx_p*nxx + txy_p*nxy)/S` (full 2-term projection). The SGS additions only apply one normal component per stress, and use raw η-stencil values (`muysgs`, `mvysgs`) as physical y-gradients (they are μ_sgs*∂/∂η, not μ_sgs*∂/∂y).

### Bug 3 — `calc_Fv_LES2_curv` lines 564-566: same class as Bug 2

SGS additions use ξ-stencil values (`muysgs` = μ_sgs*∂u/∂ξ, `mvysgs` = μ_sgs*∂v/∂ξ) directly as physical gradients. The projection is also incomplete.

### Bug 4 — `calc_Gv_LES2_curv` line 665-666: `muzsgs` is η-stencil for u, not μ_sgs*∂u/∂z

`muzsgs` is used in `tzx += (mwxsgs + muzsgs)` for τ_sgs_xz = μ_sgs*(∂w/∂x + **∂u/∂z**). At the z-face, ∂u/∂z comes from the direct z-stencil (shared memory), not the η cross-stencil.

### Bug 5 — `calc_Gv_LES2_curv` lines 671-672: `mvzsgs` is a copy of `mvysgs` (Q(:,3,...) with η-stencil)

`mvzsgs` should be μ_sgs*∂v/∂z (z-face direct), but is identical to `mvysgs`. Used in `tzy += (mwysgs + mvzsgs)` for τ_sgs_yz = μ_sgs*(∂w/∂y + **∂v/∂z**).

---

## Fix Plan

### Pattern (same for all three LES kernels)

The correct fix mirrors the molecular kernel exactly:
1. Compute all needed cross-stencils (ξ and η for Ev; η for Gv; ξ for Fv)
2. Build physical SGS gradients via chain rule (same structure as `mux = mu_f*du_dxi*xi_x_f + mu_eta_u*eta_x_f`)
3. Assemble stress tensor components `txx_sgs, txy_sgs, tyy_sgs, txz_sgs, tyz_sgs`
4. Project onto face normal: `txx += (txx_sgs*nxx + txy_sgs*nxy)/S`, etc.

### `calc_Ev_LES2_curv` — after line ~393, before "Physical gradients" block

Add η-stencil for w with SGS μ (currently missing):
```fortran
mwysgs_stencil = my1sgs*(-Q(i,4,j-1,k)-Q(i+1,4,j-1,k)) + (my1sgs-my2sgs)*(w(idx)+w(idx+1)) &
               + my2sgs*(Q(i,4,j+1,k)+Q(i+1,4,j+1,k))
```

Fix `mwxsgs` (line 411):
```fortran
! Before (wrong):  mwxsgs = ... + mwzsgs*eta_x_f
! After (correct): mwxsgs = ... + mwysgs_stencil*eta_x_f
```

Add physical SGS y-gradients (after the existing muxsgs/mvxsgs/mwxsgs block):
```fortran
real(8) mut_f_xi, muysgs_phys, mvysgs_phys, mwysgs_phys
mut_f_xi    = 0.5d0*(mut(i,j,k)+mut(i+1,j,k))
muysgs_phys = mut_f_xi*du_dxi*xi_y_f + muysgs*eta_y_f
mvysgs_phys = mut_f_xi*dv_dxi*xi_y_f + mvysgs*eta_y_f
mwysgs_phys = mut_f_xi*dw_dxi*xi_y_f + mwysgs_stencil*eta_y_f
```

Replace lines 425-427 with full projection:
```fortran
! SGS stress components
block
  real(8) txx_sgs, txy_sgs, tyy_sgs, txz_sgs, tyz_sgs
  txx_sgs = two_third*(2.d0*muxsgs   - mvysgs_phys - mwzsgs)
  txy_sgs = muysgs_phys + mvxsgs
  tyy_sgs = two_third*(2.d0*mvysgs_phys - muxsgs   - mwzsgs)
  txz_sgs = mwxsgs + muzsgs
  tyz_sgs = mwysgs_phys + mvzsgs
  txx = txx + (txx_sgs*nxx + txy_sgs*nxy) / S
  txy = txy + (txy_sgs*nxx + tyy_sgs*nxy) / S
  txz = txz + (txz_sgs*nxx + tyz_sgs*nxy) / S
end block
```

### `calc_Fv_LES2_curv` — replace lines 564-566

Add physical SGS gradients and fix projection:
```fortran
block
  real(8) mut_f_eta, muxsgs_phys, mvxsgs_phys, mwxsgs_phys
  real(8) muysgs_phys, mvysgs_phys, mwysgs_phys
  real(8) txx_sgs, txy_sgs, tyy_sgs, txz_sgs, tyz_sgs
  mut_f_eta   = 0.5d0*(mut(i,j,k)+mut(i,j+1,k))
  muxsgs_phys = muysgs*xi_x_f + mut_f_eta*du_deta*eta_x_f
  mvxsgs_phys = mvysgs*xi_x_f + mut_f_eta*dv_deta*eta_x_f
  mwxsgs_phys = mwysgs*xi_x_f + mut_f_eta*dw_deta*eta_x_f
  muysgs_phys = muysgs*xi_y_f + mut_f_eta*du_deta*eta_y_f
  mvysgs_phys = mvysgs*xi_y_f + mut_f_eta*dv_deta*eta_y_f
  mwysgs_phys = mwysgs*xi_y_f + mut_f_eta*dw_deta*eta_y_f
  txx_sgs = two_third*(2.d0*muxsgs_phys - mvysgs_phys - mwzsgs)
  txy_sgs = muysgs_phys + mvxsgs_phys
  tyy_sgs = two_third*(2.d0*mvysgs_phys - muxsgs_phys - mwzsgs)
  txz_sgs = mwxsgs_phys + muzsgs
  tyz_sgs = mwysgs_phys + mvzsgs
  tyx = tyx + (txx_sgs*nex + txy_sgs*ney) / S
  tyy = tyy + (txy_sgs*nex + tyy_sgs*ney) / S
  tyz = tyz + (txz_sgs*nex + tyz_sgs*ney) / S
end block
```

### `calc_Gv_LES2_curv` — replace lines 663-689

Add missing ξ-stencil for v, rename current `muzsgs`/`mvzsgs` to descriptive names, add correct z-face direct stencils, compute physical gradients, fix stress:
```fortran
! SGS gradients
muxsgs = (mx1sgs*(-Q(i-1,2,j,k)-Q(i-1,2,j,k+1)) + (mx1sgs-mx2sgs)*(u(idx)+u(idx+1)) &
        + mx2sgs*(Q(i+1,2,j,k)+Q(i+1,2,j,k+1)))
mu_eta_usgs = (my1sgs*(-Q(i,2,j-1,k)-Q(i,2,j-1,k+1)) + (my1sgs-my2sgs)*(u(idx)+u(idx+1)) &
             + my2sgs*(Q(i,2,j+1,k)+Q(i,2,j+1,k+1)))
mvysgs = (my1sgs*(-Q(i,3,j-1,k)-Q(i,3,j-1,k+1)) + (my1sgs-my2sgs)*(v(idx)+v(idx+1)) &
        + my2sgs*(Q(i,3,j+1,k)+Q(i,3,j+1,k+1)))
mvxsgs = (mx1sgs*(-Q(i-1,3,j,k)-Q(i-1,3,j,k+1)) + (mx1sgs-mx2sgs)*(v(idx)+v(idx+1)) &
        + mx2sgs*(Q(i+1,3,j,k)+Q(i+1,3,j,k+1)))
mwxsgs = (mx1sgs*(-Q(i-1,4,j,k)-Q(i-1,4,j,k+1)) + (mx1sgs-mx2sgs)*(w(idx)+w(idx+1)) &
        + mx2sgs*(Q(i+1,4,j,k)+Q(i+1,4,j,k+1)))
mwysgs = (my1sgs*(-Q(i,4,j-1,k)-Q(i,4,j-1,k+1)) + (my1sgs-my2sgs)*(w(idx)+w(idx+1)) &
        + my2sgs*(Q(i,4,j+1,k)+Q(i,4,j+1,k+1)))

block
  real(8) mut_f_z, muzsgs, mvzsgs, mwzsgs_direct
  real(8) muxsgs_phys, mvysgs_phys, mwxsgs_phys, mwysgs_phys
  real(8) txx_sgs, txy_sgs, tyy_sgs, txz_sgs, tyz_sgs
  mut_f_z      = 0.5d0*(mut(i,j,k)+mut(i,j,k+1))
  muzsgs       = mut_f_z*(u(idx+1)-u(idx))/dz         ! μ_sgs*∂u/∂z
  mvzsgs       = mut_f_z*(v(idx+1)-v(idx))/dz         ! μ_sgs*∂v/∂z
  mwzsgs_direct= mut_f_z*(w(idx+1)-w(idx))/dz         ! μ_sgs*∂w/∂z
  muxsgs_phys  = muxsgs*xi_x_f + mu_eta_usgs*eta_x_f  ! μ_sgs*∂u/∂x
  mvysgs_phys  = mvxsgs*xi_y_f + mvysgs*eta_y_f       ! μ_sgs*∂v/∂y
  mwxsgs_phys  = mwxsgs*xi_x_f + mwysgs*eta_x_f       ! μ_sgs*∂w/∂x
  mwysgs_phys  = mwxsgs*xi_y_f + mwysgs*eta_y_f       ! μ_sgs*∂w/∂y
  txx_sgs = two_third*(2.d0*muxsgs_phys - mvysgs_phys - mwzsgs_direct)
  txy_sgs = 0.d0  ! τ_xy not projected at z-face (normal = (0,0,1))
  txz_sgs = mwxsgs_phys + muzsgs
  tyz_sgs = mwysgs_phys + mvzsgs
  tzz_sgs = two_third*(2.d0*mwzsgs_direct - muxsgs_phys - mvysgs_phys)
  tzx = tzx + txz_sgs
  tzy = tzy + tyz_sgs
  tzz = tzz + tzz_sgs
end block
```

---

## Critical Files to Modify

| File | Lines |
|------|-------|
| `3D_solver_curv/src/calc_visc2_curv.f90` | Ev_LES2: ~line 411, 425-427; Fv_LES2: 564-566; Gv_LES2: 663-689 |

---

## Variable Declarations to Add

- `calc_Ev_LES2_curv`: add `mwysgs_stencil, mut_f_xi, muysgs_phys, mvysgs_phys, mwysgs_phys` to declarations
- `calc_Fv_LES2_curv`: the new variables are inside `block` so no declaration changes needed
- `calc_Gv_LES2_curv`: add `mu_eta_usgs, mvxsgs` to module-level declarations; new variables are inside `block`

---

## Verification

1. `cd 3D_solver_curv/CORN && make clean && make` — must compile clean (LES kernels compile even if not wired)
2. LES kernels correctness is verified when `id_visc = integer(8)` is later wired in `calc_flux_base_curv.f90`
3. The NS run (`id_visc = integer(4)`) should be unaffected — no changes to NS kernels
