# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

OUxSBLI is a GPU-accelerated CFD solver for compressible flows (Euler/Navier-Stokes), written in CUDA Fortran with MPI parallelization. It targets NVIDIA GPUs via the HPC SDK and solves test cases defined in `3D_solver/<CASE>/`.

## Build & Run

All work happens inside a specific test-case directory. There is no top-level build.

```bash
cd 3D_solver/NSTGV    # or ETGV, IVST, KHI, SBLI, TBL

make clean && make    # compile with mpif90 + CUDA flags
bash calc.sh          # run simulation (typically mpirun -n 2 a.out)
bash profile.sh       # (where available) nsys/ncu profiling
```

**Compiler requirement:** NVIDIA HPC SDK (`mpif90` with `-cuda -acc -fast -gpu=ptxinfo,rdc,lto`). Versions 24.* and 25.* are confirmed working.

**Output format:** XML VTK files, readable with ParaView.

## Architecture

Source is split across two directories; the Makefile's `vpath` merges them:

- `src/` — shared utilities: `main.f90`, `mod_constant.f90`, `cpu_gpu_mpi.f90`, `calc_muscl.f90`, `calc_physical_quantities.f90`, `print.f90`, `set_coordinate.f90`, `set_compressible_bl.f90`
- `3D_solver/src/` — solver kernels: all `calc_*.f90` scheme and time-integration modules
- `3D_solver/<CASE>/` — per-case configuration: `mod_globals.f90` (parameters), `set.f90` (grid/IC/BC), `Makefile`, `calc.sh`

### Data Flow

```
main.f90
  └─ set_coordinate()        # grid metrics & Jacobians
  └─ set() [set.f90]         # grid, IC, BC (case-specific)
  └─ calc_time_dev()         # time-loop entry point
        └─ calc_flux_base()  # dispatches to convective + viscous kernels
        └─ calc_steps()      # RK stage update (Q += dt * RHS)
        └─ calc_para()       # MPI ghost-cell exchange
        └─ print()           # VTK output
```

### Convective Schemes (dispatched from `calc_flux_base.f90`)

| Scheme | Files | Best for |
|--------|-------|----------|
| KEEP (Kinetic Energy & Entropy Preserving) | `calc_keep_kernel.f90`, `calc_keep_3d.f90` | Smooth vortical flows (TGV, KHI) |
| SLAU (Simple Low-dissipation AUSM) | `calc_slau_kernel.f90`, `calc_slau_3d.f90` | Compressible turbulence with shocks (SBLI, TBL) |
| Roe | `calc_roe_kernel.f90`, `calc_roe_3d.f90` | Supersonic/hypersonic discontinuities |
| Hybrid (KEEP↔SLAU via Ducros sensor) | `calc_hybrid_kernel.f90`, `calc_hybrid.f90` | Mixed smooth/shocked regions |

Viscous discretization: `calc_visc2.f90` (Gaitonde & Visbal 2nd-order) or `calc_visc4.f90` (4th-order compact).

## Configuration in `mod_globals.f90`

**Critical:** scheme/method selection uses the Fortran `kind` of constants as a compile-time flag—not the value:

```fortran
integer(4), parameter :: id_visc     = 2   ! kind=2→Euler, kind=4→NS, kind=8→LES
real(2),    parameter :: id_scheme   = 0   ! real(2)→KEEP, real(2) with threshold→SLAU, real(4)→Roe, real(8)→Hybrid
integer(8), parameter :: id_accuracy = 0   ! kind=2→2nd, kind=4→4th, kind=8→6th order
integer(8), parameter :: id_tvd      = 0   ! kind=2→no TVD, kind=4→minmod, kind=8→MUSCL 4th
integer(2), parameter :: id_rescale  = 0   ! kind=2→off, kind=4→on (SBLI reference state)
integer(2), parameter :: id_RungeKutta = 0 ! kind=2→TVD RK3, kind=4→RK4
integer(2), parameter :: id_recal    = 0   ! kind=2→initialize, kind=4→restart from file
```

The **value** of these parameters is ignored; only the **type kind** matters. For example, to switch from Euler to Navier-Stokes, change `integer(4)` to `integer(2)` for `id_visc`—not the value `2`.

**GPU thread blocks** are also set here as `dim3` constants (`threadsE`, `threadsF`, `threadsG`, etc.) and must be tuned to the GPU architecture and grid size.

## Conservative Variable Layout

`Q(5, nx, ny, nz)` = `[ρ, ρu, ρv, ρw, ρE]`. Flux arrays use trimmed index ranges: E-flux drops boundary points in y/z, F-flux in x/z, G-flux in x/y.

## MPI Decomposition

1D decomposition in the x-direction via `calc_para.f90`. Default is 2 MPI ranks (`mpirun -n 2 a.out`), with `mygpu = myrank / 2` (2 ranks per GPU). GPU-aware MPI is optional via `id_gpumpi`.

## Adding a New Test Case

Copy an existing case directory (e.g., `cp -r 3D_solver/NSTGV 3D_solver/MYCASE`), then edit:
1. `mod_globals.f90` — grid size, physical parameters, scheme flags
2. `set.f90` — grid generation, initial conditions, boundary condition calls
3. `calc.sh` — MPI rank count and any case-specific runtime args

## Notice
* From an occupancy perspective, the subroutines invoked within `calc_flux_base.f90` should not be executed on separate streams.
* `calc_flux_base.f90` and `calc_steps.f90` are the main bottleneck. You should optimize them.
* Roe scheme is not used. KEEP, SLAU, Hybrid schemes should be optimized.
