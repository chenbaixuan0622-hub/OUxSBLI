# Plan: Curvilinear Viscous Terms (calc_visc2_curv + load_smem_visc2_curv)

## Context

The curvilinear solver works for Euler (CORN/NACA verified). To support NS (id_visc kind=4), the viscous flux kernels must be completed. Two stub files exist with partial signatures but no math:
- `3D_solver_curv/src/calc_visc2_curv.f90` — only `calc_Ev2_curv` partial, placeholder call
- `3D_solver_curv/src/load_smem_visc2_curv.f90` — `load_smem_visc2_x_curv` stub, empty body

`calc_flux_base_curv.f90` also needs a viscous dispatch block for kind=4.

---

## Key Conventions (confirmed by reading source)

**Flux array scaling** (from `calc_steps_curv.f90` + `calc_keep_kernel_curv.f90`):
- E, F are **area-scaled**: `E(:,i,...) = flux_per_unit_area × S_xi`; step uses scalar `dt*dz`
- G is **not area-scaled**: step uses `dt_Szeta(i,j) = dt × J_2D(i+1,j+1)` per cell
- Therefore: `E(2..5,...) -= t* × S` (stub is correct), `G(2..5,...) -= t*` (no S)

**Face normals** (from `set_coordinate.f90:383-413`):
- `n_xi_x(m,n) = +0.25*(y_{m,n+2}+y_{m+1,n+2}−y_{m,n}−y_{m+1,n})` ≈ ∂y/∂η (face area vector)
- `n_xi_y(m,n) = −0.25*(x...)` ≈ −∂x/∂η
- `S = sqrt(n_xi_x²+n_xi_y²)` = face area per unit ζ-length

**Cell-center metrics** (from `set_coordinate.f90:414-426`):
- `xi_x = yeta/J2`, `xi_y = −xeta/J2`, `eta_x = −yxi/J2`, `eta_y = xxi/J2`
- Filled interior only (m=2..nx-1, n=2..ny-1); ghost cells filled by set.f90

**Physical derivative chain rule** at ξ-face (i+1/2, j):
```
∂u/∂x = ξ_x_f × δu/δξ + η_x_f × δu/δη      ξ_x_f = 0.5*(xi_x(i,j)+xi_x(i+1,j))
∂u/∂y = ξ_y_f × δu/δξ + η_y_f × δu/δη
∂u/∂z = (1/dz) × (u_{k+1}−u_{k−1})/2        (uniform z, from global Q)
```

**Projected flux variables** at ξ-face (what `txx`, `txy`, `txz` must be in the stub):
```
txx = (τ_xx×nxx + τ_xy×nxy) / S    → E(2,...) −= txx × S
txy = (τ_xy×nxx + τ_yy×nxy) / S    → E(3,...) −= txy × S
txz = (τ_xz×nxx + τ_yz×nxy) / S    → E(4,...) −= txz × S
```
Note: τ_yz = μ(∂w/∂y + ∂v/∂z) contributes to txz **only in curvilinear** (not in Cartesian E-kernel) because the ξ-face has both x and y normal components.

---

## File 1: `load_smem_visc2_curv.f90`

Near-identical copy of `3D_solver/src/load_smem_visc2.f90`. Changes only:
- Module name: `load_smem_visc2_curv`
- Function names: `load_smem_visc2_curv_x`, `load_smem_visc2_curv_y`, `load_smem_visc2_curv_z`
- `public` list updated to the three curvilinear names

The pipeline loading logic (`pipelineMemcpyAsync`, `pipelineCommit`, `pipelineWaitPrior`, `syncthreads`) is identical.

---

## File 2: `calc_visc2_curv.f90`

Replace the stub entirely. Implement all 6 subroutines (3 NS + 3 LES).

### `calc_Ev2_curv` (ξ-face, NS)

**Fix stub issues**: replace `sx=*, sy=*, sz=*` with `threadsEv%x+1, threadsEv%y, threadsEv%z`. Replace `call load_smem_visc2_curv_x(*)` with proper call.

**Thread/face mapping** (same as Cartesian Ev2):
```fortran
i = (blockIdx%x-1)*blockDim%x + it          ! i = 1..nx-1 (ξ-face index)
j = (blockIdx%y-1)*blockDim%y + jt + 1      ! j = 2..ny-1
k = (blockIdx%z-1)*blockDim%z + kt + 1      ! k = 2..nz-1
idx = (it-1) + (jt-1)*sx + (kt-1)*sx*sy
```

After shared memory load and guard:
```fortran
nxx = n_xi_x(i,j-1);  nxy = n_xi_y(i,j-1);  S = sqrt(nxx²+nxy²)
```

**Block A — metrics at face**:
```fortran
xi_x_f  = 0.5*(xi_x(i,j)+xi_x(i+1,j));    xi_y_f  = 0.5*(xi_y(i,j)+xi_y(i+1,j))
eta_x_f = 0.5*(eta_x(i,j)+eta_x(i+1,j));  eta_y_f = 0.5*(eta_y(i,j)+eta_y(i+1,j))
mu_f    = 0.5*(mu(i,j,k)+mu(i+1,j,k))
du_dxi  = u(idx+1)-u(idx);  dv_dxi = v(idx+1)-v(idx);  dw_dxi = w(idx+1)-w(idx)
dT_dxi  = T(i+1,j,k)-T(i,j,k)
```

**Block B — η-tangential stencil (4-corner mu, unit Δη)**:
```fortran
my1 = 0.0625*(mu(i,j-1,k)+mu(i,j,k)+mu(i+1,j-1,k)+mu(i+1,j,k))
my2 = 0.0625*(mu(i,j,k)+mu(i,j+1,k)+mu(i+1,j,k)+mu(i+1,j+1,k))
! pattern: my1*(-Q_jm1_i−Q_jm1_ip1) + (my1−my2)*(u_face+u_face) + my2*(Q_jp1+Q_jp1)
mu_eta_u = my1*(-Q(i,2,j-1,k)-Q(i+1,2,j-1,k)) + (my1-my2)*(u(idx)+u(idx+1)) &
           + my2*(Q(i,2,j+1,k)+Q(i+1,2,j+1,k))
mu_eta_v = (same for Q(:,3,:,:))
mu_eta_w = (same for Q(:,4,:,:))
dT_deta  = 0.25d0*((T(i,j+1,k)+T(i+1,j+1,k))-(T(i,j-1,k)+T(i+1,j-1,k)))
```

**Block C — z-tangential stencil (4-corner mu × 1/dz)**:
```fortran
mz1 = 0.0625*(mu(i,j,k-1)+mu(i,j,k)+mu(i+1,j,k-1)+mu(i+1,j,k))
mz2 = 0.0625*(mu(i,j,k)+mu(i,j,k+1)+mu(i+1,j,k)+mu(i+1,j,k+1))
muz = (mz1*(-Q(i,2,j,k-1)-Q(i+1,2,j,k-1))+(mz1-mz2)*(u(idx)+u(idx+1)) &
       +mz2*(Q(i,2,j,k+1)+Q(i+1,2,j,k+1))) / dz
mvz = (same for Q(:,3,:,:)) / dz    ! NEW vs Cartesian — needed for τ_yz = μ(∂w/∂y+∂v/∂z)
mwz = (same for Q(:,4,:,:)) / dz
```

**Block D — physical gradients via chain rule**:
```fortran
mux = mu_f*du_dxi*xi_x_f + mu_eta_u*eta_x_f   ! μ∂u/∂x
muy = mu_f*du_dxi*xi_y_f + mu_eta_u*eta_y_f   ! μ∂u/∂y
mvx = mu_f*dv_dxi*xi_x_f + mu_eta_v*eta_x_f   ! μ∂v/∂x
mvy = mu_f*dv_dxi*xi_y_f + mu_eta_v*eta_y_f   ! μ∂v/∂y
mwx = mu_f*dw_dxi*xi_x_f + mu_eta_w*eta_x_f   ! μ∂w/∂x
mwy = mu_f*dw_dxi*xi_y_f + mu_eta_w*eta_y_f   ! μ∂w/∂y  (for τ_yz)
! muz, mvz, mwz from Block C
```

**Block E — stress tensor projection and output**:
```fortran
txx_p = two_third*(2*mux - mvy - mwz)   ! τ_xx
txy_p = muy + mvx                         ! τ_xy
tyy_p = two_third*(2*mvy - mux - mwz)   ! τ_yy
txz_p = mwx + muz                         ! τ_xz
tyz_p = mwy + mvz                         ! τ_yz
txx = (txx_p*nxx + txy_p*nxy) / S
txy = (txy_p*nxx + tyy_p*nxy) / S
txz = (txz_p*nxx + tyz_p*nxy) / S

dTdx = xi_x_f*dT_dxi + eta_x_f*dT_deta
dTdy = xi_y_f*dT_dxi + eta_y_f*dT_deta
viscous_work = Cp_over_Pr*mu_f*(dTdx*nxx+dTdy*nxy)/S &
             + 0.5d0*((u(idx)+u(idx+1))*txx+(v(idx)+v(idx+1))*txy+(w(idx)+w(idx+1))*txz)

E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx*S
E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy*S
E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz*S
E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - viscous_work*S
```

---

### `calc_Fv2_curv` (η-face, NS)

Same structure as Ev2_curv but with η as the normal direction:
- Shared: `sx=threadsFv%x, sy=threadsFv%y+1, sz=threadsFv%z`
- Thread: `i=...+it+1, j=...+jt, k=...+kt+1`; `idx=(jt-1)+(it-1)*sy+(kt-1)*sy*sx`
- Load: `call load_smem_visc2_curv_y(...)` (loads along j-direction)
- Face normal: `nex=n_eta_x(i-1,j), ney=n_eta_y(i-1,j), S=sqrt(nex²+ney²)`
- Metrics at η-face: `xi_x_f=0.5*(xi_x(i,j)+xi_x(i,j+1))` etc.
- "Normal" derivative from shared: `du_deta=u(idx+1)-u(idx)`, etc.
- **ξ-cross stencil** (mx1, mx2 at ξ-face corners, unit Δξ):
  ```fortran
  mx1 = 0.0625*(mu(i-1,j,k)+mu(i,j,k)+mu(i-1,j+1,k)+mu(i,j+1,k))
  mx2 = 0.0625*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j+1,k)+mu(i+1,j+1,k))
  mu_xi_u = mx1*(-Q(i-1,2,j,k)-Q(i-1,2,j+1,k))+(mx1-mx2)*(u(idx)+u(idx+1)) &
            +mx2*(Q(i+1,2,j,k)+Q(i+1,2,j+1,k))
  ! same for v, w
  ```
- Physical gradients: `mux = mu_xi_u*xi_x_f + mu_f*du_deta*eta_x_f`, etc.
- Projection uses (nex, ney):
  ```fortran
  tyx = (τ_xx*nex + τ_xy*ney) / S   ! x-momentum at η-face
  tyy = (τ_xy*nex + τ_yy*ney) / S   ! y-momentum
  tyz = (τ_xz*nex + τ_yz*ney) / S   ! z-momentum
  ```
- Output: `F(2..5, i-1, j, k-1) -= t* × S`

---

### `calc_Gv2_curv` (z-face, NS)

z is Cartesian; face normal = (0,0,1) so no area-scaling on output.
- Shared: `sx=threadsGv%x, sy=threadsGv%y, sz=threadsGv%z+1`
- Thread: `i=...+it+1, j=...+jt+1, k=...+kt`; `idx=(kt-1)+(jt-1)*sz+(it-1)*sz*sy`
- Load: `call load_smem_visc2_curv_z(...)`
- No face normal args needed; same signature shape as Cartesian `calc_Gv2` but with xi_x/y, eta_x/y added
- Metrics at z-face: cell-center values suffice since z is uniform: `xi_x_f=xi_x(i,j)` etc. (no z-averaging)
- **z-normal derivatives** from shared:
  ```fortran
  mu_f = 0.5*(mu(i,j,k)+mu(i,j,k+1))
  muz = mu_f*(u(idx+1)-u(idx))/dz;  mvz = mu_f*(v(idx+1)-v(idx))/dz;  mwz = mu_f*(w(idx+1)-w(idx))/dz
  ```
- **ξ-cross** (at z-face corners: mu at (i-1/2, j, k+1/2)):
  ```fortran
  mx1 = 0.0625*(mu(i-1,j,k)+mu(i,j,k)+mu(i-1,j,k+1)+mu(i,j,k+1))
  mx2 = 0.0625*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j,k+1)+mu(i+1,j,k+1))
  mu_xi_u = mx1*(-Q(i-1,2,j,k)-Q(i-1,2,j,k+1))+(mx1-mx2)*(u(idx)+u(idx+1))+mx2*(Q(i+1,2,j,k)+Q(i+1,2,j,k+1))
  mu_xi_w = (same for Q(:,4,:,:))
  ```
- **η-cross** (at z-face corners: mu at (i, j-1/2, k+1/2)):
  ```fortran
  my1 = 0.0625*(mu(i,j-1,k)+mu(i,j,k)+mu(i,j-1,k+1)+mu(i,j,k+1))
  my2 = 0.0625*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j,k+1)+mu(i,j+1,k+1))
  mu_eta_v = my1*(-Q(i,3,j-1,k)-Q(i,3,j-1,k+1))+(my1-my2)*(v(idx)+v(idx+1))+my2*(Q(i,3,j+1,k)+Q(i,3,j+1,k+1))
  mu_eta_w = (same for Q(:,4,:,:))
  ```
- Physical gradients:
  ```fortran
  mux = mu_xi_u*xi_x_f + mu_eta_...*eta_x_f   ! μ∂u/∂x  (for τ_zz divergence term)
  mvy = mu_xi_v*xi_y_f + mu_eta_v*eta_y_f     ! μ∂v/∂y
  mwx = mu_xi_w*xi_x_f + mu_eta_w*eta_x_f     ! μ∂w/∂x  (for τ_zx)
  mwy = mu_xi_w*xi_y_f + mu_eta_w*eta_y_f     ! μ∂w/∂y  (for τ_zy)
  ```
- Stress (z-face, no projection):
  ```fortran
  tzx = mwx + muz;  tzy = mwy + mvz;  tzz = two_third*(2*mwz - mux - mvy)
  viscous_work = Cp_over_Pr*mu_f*(T(i,j,k+1)-T(i,j,k))/dz &
               + 0.5d0*((u(idx)+u(idx+1))*tzx+(v(idx)+v(idx+1))*tzy+(w(idx)+w(idx+1))*tzz)
  G(2,i-1,j-1,k) -= tzx;  G(3,i-1,j-1,k) -= tzy
  G(4,i-1,j-1,k) -= tzz;  G(5,i-1,j-1,k) -= viscous_work
  ```

---

### LES versions

Follow Cartesian LES pattern exactly: add `mut(nx,ny,nz)` and `qc2(nx,ny,nz)` arguments. Compute SGS versions of each stencil (muxsgs etc.) using `mut` in place of `mu`. Add SGS contribution on top of molecular stress before projection.

---

## File 3: `calc_flux_base_curv.f90`

Add NS viscous dispatch:
1. Add `use calc_visc2_curv` and add `blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv` to `use mod_globals`
2. Add `mu` as `intent(inout), device, contiguous :: mu(nx,ny,nz)` to `calc_EFG_curv` signature
3. Add overloaded `calc_visc_curv` interface (kind=2 = no-op, kind=4 = NS):
   - NS overload: call `calc_physical_quantities_mu` (or equivalent) to fill `mu` from `T`, then call the 3 viscous kernels with `<<<blocksEv, threadsEv>>>` etc.
4. Call `calc_visc_curv(id_visc, ...)` after `calc_conv_curv` in `calc_EFG_curv`
5. Update `calc_time_dev_curv.f90` to pass `mu` in the `calc_EFG_curv` call

---

## Critical Files

| File | Role |
|------|------|
| `3D_solver_curv/src/calc_visc2_curv.f90` | Primary — implement 6 kernels |
| `3D_solver_curv/src/load_smem_visc2_curv.f90` | Primary — implement 3 loaders |
| `3D_solver_curv/src/calc_flux_base_curv.f90` | Add NS dispatch |
| `3D_solver_curv/src/calc_time_dev_curv.f90` | Pass `mu` argument |
| `3D_solver/src/calc_visc2.f90` | Reference for LES stencil pattern |
| `3D_solver/src/load_smem_visc2.f90` | Reference — near-identical copy |
| `src/set_coordinate.f90:383-427` | Metric/normal definitions |
| `3D_solver_curv/src/calc_steps_curv.f90:3` | Confirms E/F area-scaled, G not |
| `3D_solver_curv/src/calc_keep_kernel_curv.f90:47` | Confirms `E = KEEP2*S` convention |

---

## Verification

1. In CORN `mod_globals.f90`: change `id_visc` to `integer(4)` (NS)
2. `cd 3D_solver_curv/CORN && make clean && make` — should compile
3. `bash calc.sh` — run; check `nohup.out` for no NaN/segfault
4. Compare density/pressure profiles with Euler run: at high Re (μ→0), NS→Euler
5. At low Re: wall boundary layer should appear (velocity decreases toward wall)
