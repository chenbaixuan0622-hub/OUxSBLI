# Bug Survey Report

## Files Surveyed

| Area | Files |
|------|-------|
| `src/` | `mod_constant.f90.fypp`, `print.f90`, `print_curv.f90`, `set_coordinate.f90`, `calc_muscl.f90`, `calc_physical_quantities.f90`, `cpu_gpu_mpi.f90` |
| `3D_solver_curv/src/` | `calc_flux_base_curv.f90.fypp`, `calc_time_dev_curv.f90.fypp`, `preprocess_curv.f90.fypp`, `main_curv.f90.fypp`, `calc_visc2_curv.f90`, `load_smem_visc2_curv.f90`, `calc_para_curv.f90`, `calc_steps_curv.f90`, `calc_keep_kernel_curv.f90`, `calc_slau_kernel_curv.f90`, `calc_hybrid_kernel_curv.f90` |
| `3D_solver/src/` | `calc_flux_base.f90.fypp`, `calc_time_dev.f90.fypp`, `preprocess.f90.fypp`, `set_bc_tbl_sbli.f90`, `set_bc_common.f90`, `calc_steps.f90`, `calc_hybrid.f90` |
| `2D_solver/src/` | `load_smem_visc4.f90`, `calc_visc_me4_base.f90`, `set_bc_common.f90` |
| Case configs | `3D_solver_curv/CORN/`, `3D_solver_curv/NACA/`, `2D_solver/OS/`, `2D_solver/BL/`, partial 3D cases |

## Not Yet Fully Surveyed

- `3D_solver/src/calc_visc4.f90.fypp`, `calc_visc4_internal.f90.fypp`
- `3D_solver/src/calc_keep_kernel*.f90.fypp`, `calc_slau_kernel*.f90.fypp`, `calc_roe_kernel*.f90.fypp`, `calc_hybrid_kernel*.f90.fypp`
- `3D_solver_curv/src/calc_visc2_curv.f90` LES section (lines 300+)
- Full `set.f90` files for each 3D solver case (ETGV, DHIT, SBLI, TBL, etc.)
- `src/set_compressible_bl.f90`

---

## Confirmed Bugs

---

### BUG 1 — `src/set_coordinate.f90` `set_Jacobian_xy3`: redundant k-loop gives wrong Jacobian for non-uniform z

**Lines:** 104–109

```fortran
do k = 2, nz-1
  do j = 2, ny-1
    do i = 2, nx-1
      Jacobian(i,j) = 8.d0 / &
      & ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)) * (dz(k-1) + dz(k)))
enddo;enddo;enddo
```

`Jacobian` is declared `real(8) :: Jacobian(nx,ny)` — a 2D array. The k-loop overwrites every `(i,j)` cell `nz-2` times. The surviving value uses `dz(nz-2) + dz(nz-1)` (from the last iteration k = nz-1).

- **Uniform z** (all current production cases): all dz values are equal, so every iteration gives the same number — harmless.
- **Non-uniform z** (e.g. a future wall-normal stretched z-grid): the Jacobian ends up wrong everywhere, silently producing incorrect flux scaling.

Compare with `set_Jacobian_xy2` (lines 83–86) which has no k-loop and correctly gives `4 / (sum_dx * sum_dy)`.

**Fix:** Remove the k-loop; use a representative dz:
```fortran
do j = 2, ny-1
  do i = 2, nx-1
    Jacobian(i,j) = 8.d0 / ((dx(i-1)+dx(i)) * (dy(j-1)+dy(j)) * 2.d0*dz(1))
enddo;enddo
```

---

### BUG 2 — `2D_solver/src/load_smem_visc4.f90`: cross-derivative shared arrays `uy/vy` and `ux/vx` left uninitialized at boundary threads

**`load_smem_visc4_x` lines 37–44; `load_smem_visc4_y` lines 78–85**

```fortran
! load_smem_visc4_x: only initialises uy, vy when 3 <= j <= ny-2
if (1 <= i .and. i <= nx .and. 3 <= j .and. j <= ny-2) then
  uy(idx) = (two_third * (-Q(i,2,j-1) + Q(i,2,j+1)) - one_twelfth * (...)) * inv_dy(j)
  vy(idx) = ...
endif
! NO else branch — uy(idx), vy(idx) are garbage when j < 3 or j > ny-2
```

The first pipeline loop (lines 28–35) initialises `u/v` with a proper else branch. The second loop for gradient arrays has no else. Threads serving boundary rows/columns carry uninitialized shared-memory gradients into the viscous-flux kernel. Identical issue in `load_smem_visc4_y` for `ux/vx` when `i < 3` or `i > nx-2`.

**Fix:** Add a zero-initialisation else branch in both loops:
```fortran
else
  uy(idx) = 0.d0
  vy(idx) = 0.d0
endif
```
(and analogously for `ux/vx` in `load_smem_visc4_y`).

---

### BUG 3 — `src/print.f90` and `src/print_curv.f90`: three non-blocking sends/receives all share the same MPI tag

**`send_recv_for_print_even2` (print.f90 ~line 256), `send_recv_for_print_even3` (~line 279), `send_recv_for_print_even_curv` (print_curv.f90 ~line 122), and their odd counterparts:**

```fortran
! Even rank (sender):
call MPI_ISEND(rho1d, ..., myrank+1, myrank+1, ...)  ! tag = myrank+1
call MPI_ISEND(p1d,   ..., myrank+1, myrank+1, ...)  ! tag = myrank+1  ← same
call MPI_ISEND(v1d,   ..., myrank+1, myrank+1, ...)  ! tag = myrank+1  ← same

! Odd rank (receiver):
call MPI_IRECV(rho1d, ..., myrank-1, myrank, ...)    ! tag = myrank
call MPI_IRECV(p1d,   ..., myrank-1, myrank, ...)    ! tag = myrank    ← same
call MPI_IRECV(v1d,   ..., myrank-1, myrank, ...)    ! tag = myrank    ← same
```

Also `print0` / `pre_calc_curv` use blocking `MPI_SEND` with the same tag for both `ke0` and `entropy0`.

The MPI standard guarantees FIFO ordering for same (source, destination, tag) messages, so the three fields arrive in the correct order today. The fragility: silently wrong data if any send or receive is ever reordered, or if a wildcard receive (`MPI_ANY_TAG`) is introduced.

**Fix:** Use distinct tags per field:
```fortran
call MPI_ISEND(rho1d, ..., myrank+1, 3*(myrank+1)-2, ...)
call MPI_ISEND(p1d,   ..., myrank+1, 3*(myrank+1)-1, ...)
call MPI_ISEND(v1d,   ..., myrank+1, 3*(myrank+1),   ...)
! and matching MPI_IRECV with 3*myrank-2, 3*myrank-1, 3*myrank
```

---

### BUG 4 — `3D_solver_curv/src/calc_visc2_curv.f90` `calc_Gv2_curv`: z-derivative variables named as η-direction derivatives

**Lines 234–256 (declaration and assignment):**

```fortran
real(8) duz_dz, dmvz_deta, dmwz_deta   ! misleading names
...
duz_dz    = mu_f * (u(idx+1) - u(idx)) / dz   ! actually mu * ∂u/∂z
dmvz_deta = mu_f * (v(idx+1) - v(idx)) / dz   ! actually mu * ∂v/∂z
dmwz_deta = mu_f * (w(idx+1) - w(idx)) / dz   ! actually mu * ∂w/∂z
```

- `duz_dz` could be read as `d(u_z)/dz` — ambiguous.
- `dmvz_deta` and `dmwz_deta` both say "η-direction" (`deta`) when they are z-direction.

The stress tensor computation (lines 285–287) is numerically correct — `mux`/`mvy` and `duz_dz`/`dmvz_deta`/`dmwz_deta` all carry the `μ` factor consistently. But the names are a maintenance hazard: a future developer reading `tzx = mwx + duz_dz` and seeing the name `duz_dz` might wrongly conclude that `duz_dz` is a pure gradient and add a `mu_f` multiplication, double-counting it.

**Fix:** Rename to reflect mu-weighting:
```fortran
real(8) mu_u_z, mu_v_z, mu_w_z   ! mu * ∂/∂z applied to u, v, w
mu_u_z = mu_f * (u(idx+1) - u(idx)) / dz
mu_v_z = mu_f * (v(idx+1) - v(idx)) / dz
mu_w_z = mu_f * (w(idx+1) - w(idx)) / dz
```

---

## Needs Deeper Verification

### POSSIBLE BUG A — `2D_solver/src/calc_visc_me4_base.f90` `calc_tau_straight`: coefficient `2.25` may be 2× too large

**Lines 22–27:**
```fortran
tmp1 = two_third * mu(1) * ((2.25d0 * (-u(2) + u(3)) - (-u(1) + u(4)) * one_twelfth) * d - ...)
```

The standard 4th-order half-point derivative stencil is `(9/8*(u3-u2) - 1/24*(u4-u1)) / dx`. Expanded: coefficient `9/8 = 1.125`.

The code uses `2.25 = 9/4`. This is correct **only if** `d = 1/(2·dx)`. If `d = 1/dx`, the coefficient is off by ×2.
**Action needed:** Read the call site in `calc_visc4.f90.fypp` to confirm what value is passed for `d`.

**Author's Comment**
You should ignore this issue now.

### POSSIBLE BUG B — `3D_solver/src/preprocess.f90.fypp`: NaN check before file read

When `RESTART = True`, the code checks `if (Qm_cpu(j) /= Qm_cpu(j))` (NaN detection) before reading from the restart file. If `Qm_cpu` is allocated but not yet initialised, this reads garbage. Need to verify execution order in context.

---

## False Positives (Investigated and Ruled Out)

| Reported | Verdict |
|----------|---------|
| `load_smem_visc2_curv_z` missing `1 <= i/j` lower bound | Not a bug — `i` and `j` are fixed inputs from CUDA 1-based thread indexing, so they are always ≥ 1 |
| `set_bc_tbl_sbli.f90` p_wall missing Jacobian multiplication | Not a bug — QJ stores Q/J already; formula is dimensionally consistent without extra J factor |
| `print.f90` velocity index `m` calculation | Correct — `m = 1 + 3*(i-1) + 3*nx*(j-1) + 3*nx*ny*(k-1)` matches the increment pattern |
| `calc_Gv2_curv` stress tensor dimensional inconsistency | Not a bug — `mux`, `mvy` are computed from mu-weighted stencils (`mu_xi_u`, `mu_eta_v`), consistent with `duz_dz` = `mu_f * du/dz` |
| `set_Jacobian_xy3` out-of-bounds array access | Not a bug — `Jacobian(nx,ny)` is 2D; loop overwrites values (wrong answer for non-uniform z) but no OOB access |
