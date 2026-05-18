# Plan: Move curvilinear metric variable declarations to main_curv.f90

## Context

Currently, 9 allocatable arrays for 2D physical coordinates and grid metrics are declared as `save` variables inside the `set` module in `3D_solver_curv/CORN/set.f90`:

```fortran
real(8), allocatable, save :: x_phys_g(:,:), y_phys_g(:,:)
real(8), allocatable, save :: n_xi_x_cpu(:,:), n_xi_y_cpu(:,:)
real(8), allocatable, save :: n_eta_x_cpu(:,:), n_eta_y_cpu(:,:)
real(8), allocatable, save :: xi_x_cpu(:,:), xi_y_cpu(:,:)
real(8), allocatable, save :: eta_x_cpu(:,:), eta_y_cpu(:,:)
real(8), allocatable, save :: Jac_cpu(:,:)
```

These variables are case-independent: every curvilinear case (CORN, NACA, future cases) needs the same set. Moving them to `main_curv.f90` as local allocatables removes the per-case duplication and makes ownership explicit. The caller (`main_curv.f90`) already passes all of them by name to `RungeKutta_curv` — this plan simply makes the declaration match that pattern.

`set_coordinate.f90::set_metrics_curv` already accepts all metric arrays as dummy arguments (not via `use set`), so **no changes to `set_coordinate.f90`** are needed.

## Files to Modify

### 1. `3D_solver_curv/CORN/set.f90`

**a) Remove module-level declarations (lines 8–13):**
Delete the 9 `allocatable, save` declarations.

**b) Update `set_grid` signature:**
Add `x_phys_g` and `y_phys_g` as `intent(out), allocatable` dummy arguments so the subroutine allocates and fills them directly into the caller's variables.

Current:
```fortran
subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
  ...
  allocate(x_phys_g(nx,ny), y_phys_g(nx,ny))
  call set_grid_c_corner(..., x_phys_g, y_phys_g, zc, dz)
```

New:
```fortran
subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz, &
                    x_phys_g, y_phys_g)
  real(8), intent(out), allocatable :: x_phys_g(:,:), y_phys_g(:,:)
  ...
  allocate(x_phys_g(nx,ny), y_phys_g(nx,ny))
  call set_grid_c_corner(..., x_phys_g, y_phys_g, zc, dz)
```

**c) Update `set_metrics` signature:**
Add all 7 metric/Jacobian arrays as `intent(out), allocatable` dummy arguments.

Current:
```fortran
subroutine set_metrics(nx, ny)
  allocate(n_xi_x_cpu(nx-1,ny-2), ...)
  ...
  call set_metrics_curv(nx, ny, x_phys_g, y_phys_g, n_xi_x_cpu, ...)
```

New:
```fortran
subroutine set_metrics(nx, ny, x_phys_g, y_phys_g, &
                       n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
                       xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, Jac_cpu)
  real(8), intent(in)  :: x_phys_g(nx,ny), y_phys_g(nx,ny)
  real(8), intent(out), allocatable :: n_xi_x_cpu(:,:), n_xi_y_cpu(:,:)
  real(8), intent(out), allocatable :: n_eta_x_cpu(:,:), n_eta_y_cpu(:,:)
  real(8), intent(out), allocatable :: xi_x_cpu(:,:), xi_y_cpu(:,:)
  real(8), intent(out), allocatable :: eta_x_cpu(:,:), eta_y_cpu(:,:)
  real(8), intent(out), allocatable :: Jac_cpu(:,:)
  allocate(n_xi_x_cpu(nx-1,ny-2), ...)
  ...
  call set_metrics_curv(nx, ny, x_phys_g, y_phys_g, n_xi_x_cpu, ...)
```

### 2. `3D_solver_curv/src/main_curv.f90`

**a) Add local declarations** after the existing `real(8), allocatable :: x(:), ...` line:
```fortran
real(8), allocatable :: x_phys_g(:,:), y_phys_g(:,:)
real(8), allocatable :: n_xi_x_cpu(:,:), n_xi_y_cpu(:,:)
real(8), allocatable :: n_eta_x_cpu(:,:), n_eta_y_cpu(:,:)
real(8), allocatable :: xi_x_cpu(:,:), xi_y_cpu(:,:)
real(8), allocatable :: eta_x_cpu(:,:), eta_y_cpu(:,:)
real(8), allocatable :: Jac_cpu(:,:)
```

**b) Update `call set_grid(...)` to pass `x_phys_g, y_phys_g`:**
```fortran
call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz, &
              x_phys_g, y_phys_g)
```

**c) Update `call set_metrics(...)` to pass all metric arrays:**
```fortran
call set_metrics(nx, ny, x_phys_g, y_phys_g, &
                 n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
                 xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, Jac_cpu)
```

The `Jacobian` conversion loop (lines 38–42) and the `RungeKutta_curv` call already reference these variables by name — no changes needed there.

**d) Update the comment** on `use set` (line 7) to remove the stale note about which vars it provides.

## No Changes Needed

- `src/set_coordinate.f90` — `set_metrics_curv` already takes all metric arrays as dummy arguments.
- `3D_solver_curv/NACA/set.f90` — does not currently declare these variables (NACA's `set` module is incomplete); when it is built out, it will not include metric variable declarations.

## Verification

```bash
cd 3D_solver_curv/CORN
make clean && make
```

Confirm the build succeeds with no undefined-variable or type-mismatch errors.
