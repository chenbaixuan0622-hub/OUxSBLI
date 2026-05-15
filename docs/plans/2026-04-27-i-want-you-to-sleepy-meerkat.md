# Plan: Create `2D_solver` Directory

## Context

The project has a `3D_solver/` directory with GPU-accelerated CUDA Fortran solver kernels and per-case subdirectories. The user wants an analogous `2D_solver/` structure. The `src/` shared utilities already contain 2D-ready routines (`calc_quantities_2D`, `print_vtk_2D`, `set_block2`, `set_Jacobian_xy2`, `set_grid_cyclic*_2D`), so `src/` needs only one targeted change. The `3D_solver/src/` kernels are inherently 3D (G-flux, 5-component Q, MPI z-ghost cells) and must be adapted for 2D.

## Key Design Decisions

- **Q layout for 2D**: `Q(4,nx,ny)` — component index first, then x, then y. Matches existing `calc_quantities_2D` / `print_vtk_2D` convention.
- **Conservative variables**: [ρ, ρu, ρv, ρE] (4 components — no w).
- **Fluxes**: `E(4,nx-1,ny-2)`, `F(4,nx-2,ny-1)` only — no G flux.
- **dt metrics**: 1D arrays `dtdy(ny-2)` (for E-flux) and `dtdx(nx-2)` (for F-flux).
- **MPI buffer**: `overlap*(ny-2)*4` (no z-ghost cells).
- **Excluded initially**: Roe scheme (CLAUDE.md: not used), LES, visc4 internals, rescale.

## Directory Structure to Create

```
2D_solver/
├── src/          ← 2D-adapted kernel files
└── CASE/         ← empty skeleton test case
    ├── data/     ← output directory
    ├── recal/    ← restart file directory
    └── ...
```

## Files to Create in `2D_solver/src/`

### 1. `preprocess.f90`
Adapt from `3D_solver/src/preprocess.f90`:
- `allocate_device_mem`: allocate `E(4,nx-1,ny-2)`, `F(4,nx-2,ny-1)` (no G), `dtdy(ny-2)`, `dtdx(nx-2)` (no dtdydz/dtdzdx/zetaz), sensor(nx,ny) 2D.
- `pre_calc`: access `Q(m,i,j)` with m=1..4 (not `Q(i,l,j,k)` 3D layout), compute `dtdy_cpu(j) = dt*dy(j)` and `dtdx_cpu(i) = dt*dx(i)`, no zetaz/dtdydz/dtdzdx, call `make_1d_for_print2` / `print_vtk_2D`.

### 2. `calc_steps.f90`
Adapt from `3D_solver/src/calc_steps.f90`:
- `calc_R`: remove G-flux term; use `dtdy` (scalar) for E, `dtdx` (scalar) for F; reduce to 4 components.
- `calc_step1`, `calc_step2_3`: 2D CUDA kernels (2D block/thread indices, no k-index); Q(4,nx,ny) layout.
- `calc_step` (RK4 stages 1-3), `calc_step4`: same adaptations; `Rs(4,nx-2,ny-2)`.

### 3. `calc_para.f90`
Adapt from `3D_solver/src/calc_para.f90`:
- All flatten/reconstruct routines: buffer size `overlap*(ny-2)*4` (no nz-6 factor); Q(4,nx,ny) indexing `Q(l,overlap+i,j+1)` (j offset +1 for ghost).
- `exchange_cyclic`, `exchange_rescale`: same MPI pattern but with 2D buffer sizes.

### 4. `calc_hybrid.f90` (Ducros sensor 2D)
Adapt from `3D_solver/src/calc_hybrid.f90`:
- `calc_Ducros` kernel: 2D thread indices; 2D divergence `du/dx + dv/dy`; vorticity magnitude `(dv/dx - du/dy)²`; sensor(nx,ny) 2D array.
- Q primitives accessed as Q(2,i,j)=u, Q(3,i,j)=v (already primitive from `calc_quantities_2D`).

### 5. `calc_keep_kernel.f90` + `calc_keep_kernel_internal.f90` + `calc_keep_2d.f90`
Adapt from `3D_solver/src/` equivalents:
- Provide `calc_keep_x`, `calc_keep_x_in`, `calc_keep_y`, `calc_keep_y_in` kernels only (no z variants).
- Q(4,nx,ny), E(4,nx-1,ny-2), F(4,nx-2,ny-1).
- `calc_keep_2d.f90`: helper file with 2D KEEP flux formulas (analog of `calc_keep_3d.f90`).

### 6. `calc_slau_kernel.f90` + `calc_slau_kernel_internal.f90` + `calc_slau_2d.f90`
Adapt similarly: x/y kernels only, Q(4,nx,ny), no w-component in flux.

### 7. `calc_hybrid_kernel.f90` + `calc_hybrid_kernel_internal.f90`
Adapt: x/y kernels only, blending KEEP+SLAU via 2D sensor(nx,ny).

### 8. `calc_visc2.f90`
Adapt from `3D_solver/src/calc_visc2.f90`:
- Provide `calc_Ev2`, `calc_Fv2` only (no `calc_Gv2`).
- Viscous stress tensor in 2D: τ_xx, τ_yy, τ_xy (no z-shear); heat flux q_x, q_y.
- Q(4,nx,ny), mu(nx,ny).

### 9. `calc_visc4.f90`
Adapt from `3D_solver/src/calc_visc4.f90` similarly (x/y only). Include `calc_visc_me4_base.f90` and `load_smem_visc4.f90` equivalent for 2D if needed.

### 10. `calc_flux_base.f90`
Adapt from `3D_solver/src/calc_flux_base.f90`:
- Remove Roe and LES (not needed initially).
- `sensor`: `real(sp), allocatable, device :: sensor(:,:)` — 2D.
- `calc_conv_keep/slau/hybrid`: remove inv_dz, G parameters; call only x/y kernels.
- `calc_EFG_Euler/visc`: use `calc_quantities_2D` / `calc_quantities_T_2D`; no G; call only Ev2/Fv2.

### 11. `calc_time_dev.f90`
Adapt from `3D_solver/src/calc_time_dev.f90`:
- `RungeKutta_3rd(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, x, dx_cpu, y, dy_cpu, Jacobian_cpu, Q)` — remove nz/z/dz.
- `allocate_device_mem` call: no G, no zetaz, no dtdydz/dtdzdx.
- `calc_step1/step2_3` kernels: pass `dtdy`/`dtdx` instead of 2D metric arrays.
- `send_recv_for_print_even` calls `make_1d_for_print2` / `print_vtk_2D`.
- Include only 3rd-order TVD RK and 4th-order RK (no rescale variants initially).

### 12. `set_bc_common.f90`
Adapt from `3D_solver/src/set_bc_common.f90`:
- Keep 2D BC routines: `set_bc_cyclic` (for periodic x/y), `set_bc_cyclic_x`, `set_bc_cyclic_y`.
- Q(4,nx,ny); remove z-BC routines.

### 13. `set_init_common.f90`
Adapt from `3D_solver/src/set_init_common.f90`:
- Keep 2D-compatible common init utilities.

## File to Edit in `src/`

### `src/main.f90`
Add a branch so `dimension == 2` calls the 2D RungeKutta (which lacks nz/z/dz args):

```fortran
if (dimension == 3) then
  call RungeKutta(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, Q)
else
  call RungeKutta(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, x, dx, y, dy, Jacobian, Q)
endif
```

The `use calc_time_dev` is resolved at compile time via the Makefile's `vpath`, so `2D_solver/src/calc_time_dev.f90` provides the 2D `RungeKutta` when building from `2D_solver/CASE/`.

Also fix the data-saving loop for 2D (currently loops over nz; for 2D this should loop over 1..1 for the z dimension, which already works since nz is set by mod_globals).

## Files to Create in `2D_solver/CASE/`

### `mod_globals.f90` (skeleton)
```fortran
module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 2
  integer(4), parameter :: id_visc = 2      ! kind=4 → NS (or kind=2 → Euler)
  real(2), parameter :: id_scheme = 0       ! real(2) → KEEP
  integer(kind=8), parameter :: id_accuracy = 0  ! kind=2 → 2nd order
  integer(kind=8), parameter :: id_tvd = 0
  integer(kind=2), parameter :: id_rescale = 0
  integer(kind=2), parameter :: id_gpumpi = 0
  ! mesh — fill in actual values
  real(8), parameter :: Lx = 1.d0, Ly = 1.d0
  integer, parameter :: nx = 65, ny = 65     ! placeholder
  integer, parameter :: nz = 1               ! dummy (2D)
  logical, parameter :: id_bc_x = .false., id_bc_y = .false., id_bc_z = .false.
  integer, parameter :: nre1 = 1, nre2 = nx, rerank = 0
  ! GPU thread blocks (2D only — no G blocks)
  type(dim3), parameter :: threadsE = dim3(32,1,1)
  type(dim3), parameter :: threadsF = dim3(32,4,1)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threads = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksEv, blocksFv, blocks
  integer, parameter :: sp = kind(1.d0)
  ! time
  integer(kind=2), parameter :: id_recal = 0
  integer(kind=2), parameter :: id_RungeKutta = 0
  integer, parameter :: step_offset = 0, start_rescale = 0
  ! physical properties — fill in actual values
  real(8), parameter :: gamma = 1.4d0, Pr = 0.71d0, Prt = 0.9d0, R = 287.03d0
  real(8), parameter :: CFL = 0.5d0
  real(8), parameter :: dt = 0.d0    ! placeholder: set based on physics
  integer, parameter :: np = 10, nt = 100
end module mod_globals
```

### `set.f90` (empty skeleton)
```fortran
module set
  use mod_globals
  use set_coordinate
  use set_bc_common
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    ! TODO: fill in grid generation
  end subroutine set_grid

  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    real(8), intent(inout) :: Q(4,nx,ny)
    ! TODO: fill in initial conditions
    Q = 0.d0
  end subroutine set_init

  subroutine set_bc(myrank, nx, ny, Jacobian, Q)
    real(8), intent(in), device :: Jacobian(nx,ny)
    real(8), intent(inout), device :: Q(4,nx,ny)
    ! TODO: fill in boundary conditions
  end subroutine set_bc
end module set
```

### `Makefile`
Mirror the NSTGV Makefile but:
- Remove Roe, LES, visc4 entries
- Use `vpath %f90 ../src:../../src` (resolves `2D_solver/src/` then `src/`)
- Adjusted OBJ list (no `calc_roe_*.o`, no `calc_visc4*.o`, no `calc_les.o`, no `calc_rescale.o`)

### `calc.sh`
```bash
nohup mpirun -n 2 a.out &
cp ./set.f90 ./data
cp ./mod_globals.f90 ./data
```

### `data/` and `recal/` directories
Empty directories needed at runtime.

## Critical Files (for reference)

- Files to create: `2D_solver/src/` (13 files), `2D_solver/CASE/` (4 files + 2 dirs)
- Files to edit: `src/main.f90` (add dimension branch for RungeKutta call)
- Files reused unchanged from `src/`: `mod_constant.f90`, `cpu_gpu_mpi.f90`, `calc_muscl.f90`, `calc_physical_quantities.f90` (already has 2D), `print.f90` (already has 2D), `set_coordinate.f90` (already has 2D)

## Verification

1. `cd 2D_solver/CASE && make clean && make` — should compile without errors (once set.f90 is filled in minimally).
2. `bash calc.sh` — should run with mpirun (after setting a real IC in set.f90).
3. VTK output in `./data/` should be readable by ParaView.
4. For a trivial test: set Q = uniform state in set_init and verify solver runs for nt steps without crashing.
