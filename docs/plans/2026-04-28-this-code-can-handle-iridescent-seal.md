# Plan: Curvilinear Grid Support for NACA 0012 Wing Simulation

## Status (as of 2026-04-28)

### Completed
| File | Notes |
|------|-------|
| `src/set_coordinate.f90` | Added `set_naca0012`, `set_grid_c_wing`, `set_metrics_curv` |
| `3D_solver/src/calc_steps_curv.f90` | 5 step kernels (calc_R_curv + 4 RK variants) |
| `3D_solver/src/calc_keep_kernel_curv.f90` | KEEP xi, eta, z kernels (2nd-order) |
| `3D_solver/src/calc_slau_kernel_curv.f90` | SLAU xi, eta, z kernels (2nd-order) |
| `3D_solver/src/calc_hybrid_curv.f90` | Ducros sensor (`calc_Ducros_curv`) + Albada/sigmoid/wiggle helpers |
| `3D_solver/src/calc_hybrid_kernel_curv.f90` | Hybrid xi, eta, z kernels |
| `3D_solver/src/preprocess.f90` | `allocate_device_mem_curv` + `pre_calc_curv` added (has bugs, see §4) |
| `3D_solver/src/calc_time_dev_curv.f90` | `RungeKutta_curv` TVD-RK3 (has bugs, see §2) |
| `src/print_curv.f90` | `print_vtk_curv` (VTK StructuredGrid) |
| `3D_solver/WING/mod_globals.f90` | nx=256, ny=96, nz=4, Ma=0.3, Euler, Hybrid |
| `3D_solver/WING/set.f90` | `set_grid`, `set_metrics`, `set_init`, `set_bc` (has bug, see §5) |
| `3D_solver/WING/main.f90` | Program main |
| `3D_solver/WING/Makefile` | Build config (has bugs, see §6) |
| `3D_solver/WING/calc.sh` | Run script |

---

## Remaining Work

### 1. Create `3D_solver/WING/calc_flux_base_curv.f90` **(CRITICAL — blocking)**

`calc_time_dev_curv.f90` calls `calc_EFG_curv` which does not exist anywhere.
The Makefile already lists `calc_flux_base_curv.o` but the source is missing.

```fortran
module calc_flux_base_curv
  use mod_globals, only : id_scheme, id_accuracy, id_slau, sp, blocks, threads, &
                          blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
  use calc_physical_quantities
  use calc_keep_kernel_curv
  use calc_slau_kernel_curv
  use calc_hybrid_kernel_curv
  use calc_hybrid_curv
  implicit none
  private
  public init_sensor_curv, calc_EFG_curv

  real(sp), allocatable, device, save :: sensor(:,:,:)
  real(8),  allocatable, device, save :: Q_prim(:,:,:,:)   ! primitive variables workspace
contains

  subroutine init_sensor_curv(nx, ny, nz)
    integer, intent(in) :: nx, ny, nz
    allocate(sensor(nx,ny,nz), Q_prim(nx,5,ny,nz))
  end subroutine init_sensor_curv

  !> Euler curvilinear flux (id_visc kind=2).
  !> QJ = conservative vars / Jac_code; Jacobian = 1/(J_2D*dz).
  subroutine calc_EFG_curv(id_visc, nx, ny, nz, dz, &
      n_xi_x, n_xi_y, n_eta_x, n_eta_y, &
      xi_x, xi_y, eta_x, eta_y, Jacobian, &
      QJ, T, E, F, G)
    integer(2), intent(in), value :: id_visc
    integer,    intent(in), value :: nx, ny, nz
    real(8),    intent(in), value :: dz
    real(8),    intent(in), device, contiguous :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8),    intent(in), device, contiguous :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8),    intent(in), device, contiguous :: xi_x(nx,ny), xi_y(nx,ny)
    real(8),    intent(in), device, contiguous :: eta_x(nx,ny), eta_y(nx,ny)
    real(8),    intent(in), device, contiguous :: Jacobian(nx,ny)
    real(8),    intent(in), device, contiguous :: QJ(nx,5,ny,nz)
    real(8),    intent(out),device, contiguous :: T(nx,ny,nz)
    real(8),    intent(out),device, contiguous :: E(5,nx-1,ny-2,nz-2)
    real(8),    intent(out),device, contiguous :: F(5,nx-2,ny-1,nz-2)
    real(8),    intent(out),device, contiguous :: G(5,nx-2,ny-2,nz-1)

    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q_prim, T)
    call calc_conv_curv(id_scheme, nx, ny, nz, dz, &
        n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, &
        Q_prim, T, E, F, G)
  end subroutine calc_EFG_curv

  ! Dispatch by kind of id_scheme using three overloaded procedures:
  !   integer(2) → calc_conv_curv_keep    (KEEP)
  !   real(2)    → calc_conv_curv_slau    (SLAU)
  !   real(8)    → calc_conv_curv_hybrid  (Hybrid)
  ! Pattern: same as calc_conv interface in calc_flux_base.f90 (lines 31-33)
  interface calc_conv_curv
    module procedure calc_conv_curv_keep, calc_conv_curv_slau, calc_conv_curv_hybrid
  end interface
```

For KEEP variant (no sensor needed):
```fortran
  call calc_keep_xi_curv<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, n_xi_x, n_xi_y, Q, T, E)
  call calc_keep_eta_curv<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, n_eta_x, n_eta_y, Q, T, F)
  call calc_keep_z_curv<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, T, G)
```

For SLAU variant:
```fortran
  call calc_Ducros_curv<<<blocks,threads>>>(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, Q, sensor)
  call calc_slau_xi_curv<<<blocksE,threadsE>>>(id_accuracy, nx, ny, nz, n_xi_x, n_xi_y, Q, sensor, E)
  call calc_slau_eta_curv<<<blocksF,threadsF>>>(id_accuracy, nx, ny, nz, n_eta_x, n_eta_y, Q, sensor, F)
  call calc_slau_z_curv<<<blocksG,threadsG>>>(id_accuracy, nx, ny, nz, Q, sensor, G)
```

For Hybrid variant: same as SLAU but use `calc_hybrid_xi/eta/z_curv` and pass T too.

---

### 2. Fix `3D_solver/src/calc_time_dev_curv.f90`

**Bug A** (line 9): `use calc_flux_base` → `use calc_flux_base_curv`  
Also add `use print_curv` to the use list.

**Bug B** (3 places, lines 92, 101, 111): Add `Jacobian` to the `calc_EFG_curv` call:
```fortran
! Change:
call calc_EFG_curv(id_visc, nx, ny, nz, dz_val, &
    n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, &
    QJ, T, E, F, G)
! To:
call calc_EFG_curv(id_visc, nx, ny, nz, dz_val, &
    n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, Jacobian, &
    QJ, T, E, F, G)
```

**Bug C** (lines 124, 126): `send_recv_for_print_even`/`send_recv_for_print_odd` do not exist.
Change to `send_recv_for_print_even_curv` / `send_recv_for_print_odd_curv` (added in §3).

---

### 3. Add print helpers to `src/print_curv.f90`

The existing `send_recv_for_print_even3`/`odd3` in [src/print.f90:288-332](src/print.f90#L288-L332) take 1D `x(nx), y(ny)` and call `print_vtk` (RectilinearGrid). The curvilinear case needs 2D coords.

Add to `src/print_curv.f90` (follow the `even3`/`odd3` pattern exactly):

```fortran
subroutine send_recv_for_print_even_curv(myrank, nranks, step, nx, ny, nz, &
    x_phys, y_phys, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
  ! Q = QJ (device→host); make_1d_for_print; MPI_ISEND ×3; MPI_WAITALL

subroutine send_recv_for_print_odd_curv(myrank, nranks, step, nx, ny, nz, &
    x_phys, y_phys, z, Jacobian_cpu, Q, ke0, entropy0)
  ! MPI_IRECV ×3; MPI_WAITALL; print_vtk_curv(step, nx, ny, nz, myrank, nranks, ...)
```

---

### 4. Fix `3D_solver/src/preprocess.f90` (inside `allocate_device_mem_curv` and `pre_calc_curv`)

**Bug A** ([preprocess.f90:213](3D_solver/src/preprocess.f90#L213)):
```fortran
! Change:
use calc_flux_base, only : init_sensor_curv
! To:
use calc_flux_base_curv, only : init_sensor_curv
```

**Bug B** ([preprocess.f90:340](3D_solver/src/preprocess.f90#L340)):
`call print_vtk(0, nx, ny, nz, myrank+1, nranks, x_phys, y_phys, z, ...)` passes 2D `x_phys`
to a function expecting 1D `x(nx)`. Fix:
```fortran
! Add local USE inside pre_calc_curv:
use print_curv, only : print_vtk_curv
! Change call to:
call print_vtk_curv(0, nx, ny, nz, myrank+1, nranks, x_phys, y_phys, z, rho1d, p1d, v1d, ke0, entropy0)
```

---

### 5. Fix `3D_solver/WING/set.f90` — Out-of-bounds in `set_bc_cyclic_x`

`Q` is `Q(nx,5,ny,nz)` but line 101 accesses `Q(nx+1,:,:,:)` (out of bounds).

Correct ghost-cell convention (interior cells are i=2..nx-1 per calc_step_curv):
```fortran
subroutine set_bc_cyclic_x(nx, ny, nz, Q)
  integer, intent(in), value :: nx, ny, nz
  real(8), intent(inout), device :: Q(nx,5,ny,nz)
  Q(1,:,:,:)  = Q(nx-1,:,:,:)   ! left ghost ← last interior
  Q(nx,:,:,:) = Q(2,:,:,:)      ! right ghost ← first interior
end subroutine set_bc_cyclic_x
```

Also remove line 129 in `set_bc`: `Q(:,:,:,nz+1) = Q(:,:,:,2)` (out of bounds for `Q(:,:,:,nz)`).

---

### 6. Fix `3D_solver/WING/Makefile`

`main_curv.o` is listed but source is `main.f90`. Add missing objects.

```makefile
OBJ = mod_globals.o \
      mod_constant.o \
      cpu_gpu_mpi.o \
      set_coordinate.o \
      calc_physical_quantities.o \
      calc_muscl.o \
      calc_hybrid.o \
      calc_hybrid_curv.o \
      calc_keep_kernel_curv.o \
      calc_slau_kernel_curv.o \
      calc_hybrid_kernel_curv.o \
      calc_steps_curv.o \
      set_bc_common.o \
      set.o \
      calc_flux_base_curv.o \
      print.o \
      print_curv.o \
      preprocess.o \
      calc_para.o \
      calc_time_dev_curv.o \
      main.o
```

Key dependency lines:
```makefile
calc_flux_base_curv.o: mod_globals.mod calc_physical_quantities.mod \
    calc_hybrid_curv.mod calc_keep_kernel_curv.mod \
    calc_slau_kernel_curv.mod calc_hybrid_kernel_curv.mod
preprocess.o: mod_globals.mod print.mod print_curv.mod calc_flux_base_curv.mod
calc_time_dev_curv.o: mod_globals.mod preprocess.mod set.mod \
    calc_flux_base_curv.mod calc_steps_curv.mod calc_para.mod print_curv.mod
main.o: mod_globals.mod set_coordinate.mod set.mod calc_time_dev_curv.mod
```

---

## Implementation Order

1. Create `3D_solver/WING/calc_flux_base_curv.f90`
2. Fix `calc_time_dev_curv.f90` (use/call changes)
3. Add `send_recv_for_print_even/odd_curv` to `src/print_curv.f90`
4. Fix `preprocess.f90` (2 targeted lines)
5. Fix `set.f90` out-of-bounds
6. Fix `Makefile`
7. `make clean && make` — resolve remaining errors

---

## Key Conventions

- **Area-scaled fluxes**: E = KEEP2(n̂) × |S_ξ|, F = KEEP2(n̂) × |S_η|; G is NOT area-scaled
- **dt factors**: `dt_xi = dt_eta = dt×dz` (scalars); `dt_Szeta(i,j) = dt / (Jac_code(i,j) × dz)`
- **Jacobian convention**: `Jac_code(i,j) = 1/(J_2D(i,j)×dz)` so `QJ = Q / Jac_code`
- **Face normals**: `n_xi_x = +∂y/∂η`, `n_xi_y = -∂x/∂η`; `n_eta_x = -∂y/∂ξ`, `n_eta_y = +∂x/∂ξ`
- **O-grid ghost cells**: i=1 ghost ← i=nx-1 (last interior); i=nx ghost ← i=2 (first interior)
- **Kind-based dispatch**: `integer(2)`→KEEP, `real(2)`→SLAU, `real(8)`→Hybrid (value ignored)

---

## Verification

1. `cd 3D_solver/WING && make clean && make` — should compile without errors
2. `bash calc.sh` — runs 2 MPI ranks; check `results/` for VTK output
3. Open `results/rank0/*.vts` in ParaView — should show structured O-grid around airfoil
4. Print `minval(Jac_cpu)` in `set_metrics` — must be positive everywhere
5. Euler inviscid AoA=0, Ma=0.3: pressure should be symmetric about chord line
