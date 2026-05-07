# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

OUxSBLI is a GPU-accelerated CFD solver for compressible flows (Euler/Navier-Stokes), written in CUDA Fortran with MPI parallelization. It targets NVIDIA GPUs via the HPC SDK and solves test cases defined in `3D_solver/<CASE>/`.

## Build & Run

All work happens inside a specific test-case directory. There is no top-level build.

```bash
module load nvhpc-openmpi3        # load NVIDIA HPC SDK + MPI (required before compiling)
cd 3D_solver/NSTGV                # or ETGV, IVST, KHI, SBLI, TBL

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

Viscous discretization: `calc_visc2.f90` (Gaitonde & Visbal 2nd-order) or `calc_visc4.f90` (4th-order compact), plus `load_smem_visc2.f90` / `load_smem_visc4.f90` for async shared-memory loaders.

### calc_visc2.f90 — 2nd-order viscous kernels

6 kernels: `calc_Ev2/Fv2/Gv2` (NS) and `calc_Ev_LES2/Fv_LES2/Gv_LES2` (LES).

**Shared memory:** `load_smem_visc2_{x,y,z}` fills `u/v/w` with `Q(:,2:4,:,:)` (momentum ρu, ρv, ρw) via `cp.async` pipeline. `syncthreads()` is called inside the loader. After the call, `u(idx)` and `u(idx+1)` are the two face-adjacent values in the primary direction.

**Index layout:**
- E-kernel: `idx = (it-1) + (jt-1)*sx + (kt-1)*sx*sy`, primary direction = x, `idx+1` = i+1
- F-kernel: `idx = (jt-1) + (it-1)*sy + (kt-1)*sy*sx`, primary direction = y, `idx+1` = j+1
- G-kernel: `idx = (kt-1) + (jt-1)*sz + (it-1)*sz*sy`, primary direction = z, `idx+1` = k+1

**Variable naming:**
- `mudx` = `μ * (1/Δx)` — pre-multiplied coefficient; reduces 4 separate `* dx` multiplications to 1
- `mux` = `mudx * Δ(ρu)` (not velocity gradient — works on momentum directly)
- `my1, my2` = face-averaged μ at j-½ and j+½ cross-direction stencil points (scalars, not arrays)
- `mxsgsdx` / `mysgsdy` / `mzsgsdz` = SGS viscosity pre-multiplied by grid spacing
- `mysgs1/mysgs2` = cross-direction SGS viscosity stencil scalars

**LES extras:** `mut` = turbulent (SGS) viscosity, `qc2` = turbulent kinetic energy. LES kernels compute both molecular (`my1/my2`) and SGS (`mysgs1/mysgs2`) cross-direction stresses in parallel blocks.

### calc_visc4.f90 — 4th-order viscous kernels

6 kernels with the same naming pattern. Each kernel has two branches:
- **Interior** (`3 ≤ i ≤ nx-3`, etc.): uses 4th-order stencil via `calc_tau_straight` / `calc_tau_cross` / `flux4` (defined in `calc_visc_me4_base.f90`, included at compile time). Uses `mu3(3)` array for 3-point face viscosity.
- **Boundary fallback**: uses 2nd-order stencil identical in structure to `calc_visc2`.

`load_smem_visc4_{x,y,z}` loads `u/v/w` plus cross-direction derivative arrays `uy/vy/uz/wz` (pre-computed y- and z-derivatives of momentum) into shared memory.

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

`Q(nx, 5, ny, nz)` = `[ρ, ρu, ρv, ρw, ρE]`. Flux arrays use trimmed index ranges: E-flux drops boundary points in y/z, F-flux in x/z, G-flux in x/y. AoS and SoA hybrid memory layout is used.

## MPI Decomposition

1D decomposition in the x-direction via `calc_para.f90`. Default is 2 MPI ranks (`mpirun -n 2 a.out`), with `mygpu = myrank / 2` (2 ranks per GPU). GPU-aware MPI is optional via `id_gpumpi`.

## Adding a New Test Case

Copy an existing case directory (e.g., `cp -r 3D_solver/NSTGV 3D_solver/MYCASE`), then edit:
1. `mod_globals.f90` — grid size, physical parameters, scheme flags
2. `set.f90` — grid generation, initial conditions, boundary condition calls
3. `calc.sh` — MPI rank count and any case-specific runtime args

## GPU Kernel Optimization Rules (NVHPC / CUDA Fortran)

### Stack frame (local memory) — the key bottleneck

NVHPC places variables in **local memory** (per-thread stack, slow) instead of registers when:
- The variable is a **Fortran array** of any size (`dimension(n)`) — CUDA registers are not addressable, so arrays always spill to local memory.
- The variable is declared inside a **`block` construct** — NVHPC typically assigns a stack slot to block-scoped variables even if they are scalar, because it treats each `block` scope separately during register allocation.

**Consequences for this codebase:**
- All `dimension(2)` cross-direction arrays (`mysgs(2)`, `my(2)`, `H(2)`, etc.) in LES kernels force 16-byte stack slots. Replace with scalar pairs `mysgs1, mysgs2`.
- Stress tensor variables (`txx/txy/txz` etc.) declared in a final `block` construct add 24 bytes of stack. Declare them at subroutine scope instead and remove the `block`/`end block` wrapper.
- Purely temporary scalars that are only needed inside one block (e.g., `mudx`, `my1/my2`, `mx1/mx2`) can remain block-scoped — the compiler reuses the stack slot across sequential non-overlapping blocks, keeping peak frame small.

**Diagnosis:** Compile with `-Minfo=accel -gpu=ptxinfo` (already in the Makefiles). Look for `Stack frame=N bytes` per kernel in the compiler output. Target: 0 bytes for non-LES kernels, minimal for LES.

### Intent(out) scalars passed to device subroutines

Fortran passes `intent(out)` scalars by reference. In CUDA Fortran, this forces the **caller** to allocate addressable local memory for each such scalar (8 bytes each). If a device subroutine has 4 `intent(out)` scalars, the caller gets a 32-byte stack frame regardless of what the subroutine does.
- **`load_smem_visc2_*`** was refactored to remove all `intent(out)` scalars — it now only fills shared memory arrays and returns nothing.
- LTO cannot inline away this stack frame when the callee uses `pipelineMemcpyAsync/pipelineCommit/pipelineWaitPrior` (from `use wmma`) because the compiler treats those as opaque side effects.

### Register vs. stack allocation heuristics

| Pattern | Result |
|---------|--------|
| Scalar declared at subroutine scope | Register (compiler can allocate freely) |
| Scalar declared inside `block` | May spill to stack (NVHPC limitation) |
| `dimension(n)` array at any scope | Always local memory |
| `intent(out)` scalar in device subroutine call | 8 bytes of caller stack per scalar |
| Sequential non-overlapping `block` constructs | Compiler reuses the same stack slot across them |

## Notice
* From an occupancy perspective, the subroutines invoked within `calc_flux_base.f90` should not be executed on separate streams.
* `calc_flux_base.f90` and `calc_steps.f90` are the main bottleneck. You should optimize them.
* `calc_flux_base.f90` calls `calc_*_kernel.f90`, `calc_*_kernel_internal.f90`, and `calc_visc*.f90`. They are the main bottleneck.
* Roe scheme is not used. KEEP, SLAU, Hybrid schemes should be optimized.
* `id_accuracy` is not only for convection terms but also for viscous terms because it controls the size of the ghost cells.
* Do not add `contiguous` and `shared` attributes when passing shared memory as an argument.
* `cpu_gpu_mpi.f90` will be modified in the future.

## Strict Tooling Rules
- MCP servers are strictly prohibited.
- Never suggest or attempt to use MCP tools.
- All code analysis must be done by reading and reasoning over the code.
- Even if tools are available, ignore them completely.
