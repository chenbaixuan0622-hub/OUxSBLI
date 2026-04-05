# OUxSBLI System Architecture

## 1. Overview

OUxSBLI is a GPU-accelerated Computational Fluid Dynamics (CFD) solver with a modular architecture designed for high performance on NVIDIA GPUs. The system is organized around three core principles:

1. **Spatial/Temporal Separation**: Discretization methods are cleanly separated from time integration
2. **GPU-First Design**: CUDA Fortran kernels are primary; CPU code supports coordination
3. **Configurable Accuracy**: Multi-order accuracy schemes selection at compile time

---

## 2. Top-Level Architecture

```
┌─────────────────────────────────────────────────────────┐
│         Application Layer (Simulation Cases)            │
│  ETGV, IVST, KHI, NSTGV, SBLI, TBL, etc.                │
└────────────┬────────────────────────────────────────────┘
             │
┌────────────▼────────────────────────────────────────────┐
│         Time Integration Layer                          │
│  calc_time_dev.f90                                      │
│  - RungeKutta_3rd/4th (with/without rescaling)          │
│  - Per-stage: flux computation → RHS → Q update         │
└────────────┬────────────────────────────────────────────┘
             │
┌────────────▼────────────────────────────────────────────┐
│    Spatial Discretization Layer (Per Stage)             │
│  ┌──────────────────────────────────────────────────┐   │
│  │ Convection Fluxes    │ Viscous Fluxes            │   │
│  │ calc_keep_*_in       │ calc_visc2/4              │   │
│  │ calc_slau_*_in       │ LES mutation sync         │   │
│  │ calc_roe_*           │                           │   │
│  └──────────────────────────────────────────────────┘   │
└────────────┬────────────────────────────────────────────┘
             │
┌────────────▼────────────────────────────────────────────┐
│         GPU Kernel Layer (CUDA Fortran)                 │
│  - Shared memory optimization                           │
│  - Coalesced global memory access                       │
│  - Multi-dimensional thread blocks                      │
└────────────┬────────────────────────────────────────────┘
             │
┌────────────▼────────────────────────────────────────────┐
│         Device Memory & Communication                   │
│  - Q, E, F, G arrays (device)                           │
│  - MPI ghost region exchange (GPU-aware)                │
│  - Boundary condition application                       │
└─────────────────────────────────────────────────────────┘
```

---

## 3. Module Organization

### 3.1 Tier 1: Problem Setup (Case-Specific)

**`mod_globals.f90`** (per case)
- Parameter declarations: grid size, scheme type, accuracy order
- GPU thread/block configuration
- Physical constants (depending on case)
- Role: Compile-time configuration

**`set.f90`** (per case)
- Grid generation and coordinate system setup
- Initial condition functions
- Boundary condition implementation
- Role: Problem specification

### 3.2 Tier 2: Core Physics Modules (Shared)

**`mod_constant.f90`**
- Physical constants (γ, R, Pr, etc.)
- Conversion factors and mathematical constants
- Role: Global physics parameters

**`calc_physical_quantities.f90`**
- Temperature from internal energy
- Dynamic viscosity (Sutherland's law)
- Speed of sound
- Role: Primitive variable computations

**`set_coordinate.f90`**
- Jacobian determinant for curvi linear grids
- Coordinate metrics (∂ξ/∂x, ∂η/∂y, etc.)
- Grid spacing (dx, dy, dz) handling
- Role: Grid metric computations

### 3.3 Tier 3: Convection Flux Layer

**`calc_keep_kernel_internal.f90`** (Optimized KEEP)
- Unified 2nd/4th/6th-order KEEP kernels
- `calc_keep_x_in`, `calc_keep_y_in`, `calc_keep_z_in`
- Calls include 'calc_keep_3d.f90' for flux functions
- Role: Kinetic energy preserving high-accuracy schemes

**`calc_slau_kernel_internal.f90`** (Optimized SLAU)
- Multi-accuracy interpolation (interp2/4/6)
- Shock-adaptive flux with wiggle detector
- `calc_slau_x_in`, `calc_slau_y_in`, `calc_slau_z_in`
- Calls include 'calc_slau_3d.f90' for flux functions
- Role: Low-dissipation shock-adaptive schemes

**`calc_roe_kernel.f90`**
- Classical Roe scheme
- Role: Robust flux splitting for shocks

**`calc_hybrid_kernel.f90`**
- Hybrid KEEP/SLAU/Roe switching
- Role: Adaptively combine schemes

### 3.4 Tier 4: Viscous & LES Layer

**`calc_visc2.f90`**
- 2nd-order viscous fluxes
- `calc_Ev2()`, `calc_Fv2()`, `calc_Gv2()`
- LES variants: `calc_Ev_LES2()`, etc.

**`calc_visc4.f90`**
- 4th-order (ME4-Base) viscous fluxes
- `calc_Ev4()`, `calc_Fv4()`, `calc_Gv4()`
- LES variants available

**`calc_les.f90`**
- Large Eddy Simulation turbulence model
- Selective mixed scale approach
- Role: SGS (subgrid-scale) dissipation

### 3.5 Tier 5: Boundary Conditions

**`set_bc_common.f90`**
- Periodic (cyclic) BC kernels for all accuracy orders
- `set_bc_cyclic2()`, `set_bc_cyclic4()`, `set_bc_cyclic6()`
- LES mutation synchronization: `set_bc_mut_common()`
- Role: BC application after RK stages

**`set_bc_tbl_sbli.f90`** (case-specific)
- Wall BC for turbulent boundary layer
- Oblique shock BC for SBLI
- Role: Non-trivial BC implementations

### 3.6 Tier 6: Time Integration

**`calc_time_dev.f90`**
- RK3-TVD: `RungeKutta_3rd`, `RungeKutta_3rd_rescale`
- RK4: `RungeKutta_4th`, `RungeKutta_4th_rescale`
- Role: Time stepping orchestration

**`calc_rescale.f90`**
- Mean flow rescaling for SBLI
- `calc_mean()`: Accumulate statistics
- `step_rescale()`: Apply transformation
- Role: Shock-turbulence equilibrium state tracking

### 3.7 Tier 7: Communication & Utilities

**`cpu_gpu_mpi.f90`**
- MPI send/recv with GPU memory
- CPU route (pinned memory): `CPU_MPI_SEND`, `CPU_MPI_ISEND`
- GPU route (CUDA-aware): `GPU_MPI_SEND`, `GPU_MPI_ISEND`
- Role: GPU-MPI abstraction layer

**`calc_muscl.f90`**
- MUSCL reconstruction for interpolation
- Role: High-order primitive reconstruction

**`print.f90`**
- Output formatting and I/O
- Role: Data serialization

---

## 4. Data Flow Architecture

### 4.1 State Representation

```fortran
! Primary state vector
Q(1, i, j, k) = ρ         (density)
Q(2, i, j, k) = ρu        (x-momentum)
Q(3, i, j, k) = ρv        (y-momentum)
Q(4, i, j, k) = ρw        (z-momentum)
Q(5, i, j, k) = E or p    (energy or pressure, config-dependent)
```

### 4.2 Flux Computation Pipeline

```
Input: Q (conservative), T (temperature)
  ↓
[Convection Fluxes]
  - KEEP/SLAU/Roe kernel
  - Output: E(5,nx-1,ny-2,nz-2), F(), G()
  ↓
[Apply Boundary Conditions]
  - set_bc_cyclic() kernel
  ↓
[Viscous Fluxes]
  - calc_visc2/4() kernels
  - Output: E_v(), F_v(), G_v()
  - LES: Add mut, qc2 contributions
  ↓
[Compute RHS]
  - ∂Q/∂t = -(∂E/∂x + ∂F/∂y + ∂G/∂z)
  ↓
[Update Q]
  - Runge-Kutta stage update
```

### 4.3 GPU Memory Hierarchy

```
┌─────────────────────────────────────┐
│  Host (CPU) Memory                  │
│  ├─ Q_cpu (pinned for MPI)          │
│  ├─ Jacobian_cpu                    │
│  ├─ Coordinate arrays               │
│  └─ Checkpoint I/O buffers          │
└────────────┬────────────────────────┘
             │ cudaMemcpy
┌────────────▼────────────────────────┐
│  Device (GPU) Memory                │
│  ├─ Q (state vector)                │
│  ├─ E, F, G (fluxes)                │
│  ├─ T, mu (temperature, viscosity)  │
│  ├─ mut, qc2 (LES)                  │
│  ├─ dx, dy, dz (device copies)      │
│  ├─ Jacobian (2D metric)            │
│  └─ Temp arrays (stage storage)     │
│                                     │
│  ┌──────────────────────────────┐   │
│  │  Shared Memory (per block)   │   │
│  │  ├─ Stencil data             │   │
│  │  ├─ Flux components          │   │
│  │  └─ Reduction results        │   │
│  └──────────────────────────────┘   │
│                                     │
│  ┌──────────────────────────────┐   │
│  │  Registers (per thread)      │   │
│  │  ├─ Loop indices             │   │
│  │  ├─ Local state              │   │
│  │  └─ Intermediate values      │   │
│  └──────────────────────────────┘   │
└─────────────────────────────────────┘
```

---

## 5. Design Patterns

### 5.1 Parameter Selection Pattern

**Problem**: How to support multiple accuracy orders without code duplication?

**Solution**: Use Fortran `kind()` intrinsic to encode accuracy

```fortran
integer(kind=2), parameter :: id_accuracy = 0  ! 2nd-order
integer(kind=4), parameter :: id_accuracy = 0  ! 4th-order
integer(kind=8), parameter :: id_accuracy = 0  ! 6th-order

! Then: io = kind(id_accuracy) / 3
! Results: io = 0 (2nd), io = 1 (4th), io = 2 (6th)
```

**Benefits**:
- Single source code, three accuracy options
- Compile-time optimization (constants eliminated)
- No runtime branching overhead

### 5.2 Shared Memory Tiling Pattern

**Problem**: GPU kernels need to load neighboring values for stencils

**Solution**: Divide thread blocks with overlap for stencil access

```fortran
! KEEP kernel with 6th-order stencil
integer, parameter :: sx = threadsE%x + 2*io + 1

! Shared memory indexed: -(io-1):sx*sy*sz-io
real(8), shared :: rho(-(io-1):sx*sy*sz-io)

! Load: each thread loads from its position ± stencil radius
do ii = it-io, threadsE%x+io+1, blockDim%x
  rho(ii + offset_yz) = Q(1, i_base+ii, j, k)
enddo
call syncthreads()

! Compute: only safe threads use full stencil
if (io+1 <= i .and. i <= nx-(io+1)) then
  flux = KEEP(rho(idx-io:idx+io+1), ...)
endif
```

**Benefits**:
- Efficient cache utilization
- Coalesced memory access
- Reduces global memory bandwidth by 10-100x

### 5.3 Multi-Phase Computation Pattern (SLAU)

**Problem**: SLAU interpolation and flux need coordinated memory use

**Solution**: Split into two synchronized phases

```fortran
Phase 1: Interpolation
  - Load stencil data: larger shared memory buffer
  - Compute left/right states
  - Store results: smaller dedicated shared memory
  - sync_threads()

Phase 2: Flux Computation
  - Reuse shared memory from Phase 1 (freed up)
  - Load shock sensor
  - Compute SLAU flux
  - Write output
```

**Benefits**:
- Minimizes shared memory footprint
- Avoids stack frame overhead
- Enables higher occupancy

### 5.4 Abstraction for GPU-Aware MPI

**Problem**: Not all systems have CUDA-aware MPI

**Solution**: Interface with runtime selection

```fortran
interface CPUGPU_MPI_SEND
  module procedure CPU_MPI_SEND, GPU_MPI_SEND
end interface

! Selected by id_gpu_mpi parameter
if (id_gpu_mpi == 2) then
  call CPU_MPI_SEND(...)   ! Pin memory, MPI from host
else
  call GPU_MPI_SEND(...)   ! Direct GPU MPI
endif
```

**Benefits**:
- Portability across HPC centers
- Optimization opportunity for newer systems

---

## 6. Execution Model

### 6.1 Per-Timestep Execution

```
DO step = 1, nt
  │
  ├─ [Stage 1] RK coefficient α₁
  │   ├─ calc_physical_quantities(Q → T, μ)
  │   ├─ calc_EFG_Euler(Q, T → E, F, G)
  │   ├─ calc_EFG_Visc(Q, T, μ → ΔE, ΔF, ΔG)
  │   ├─ Compute RHS: dQ = -∇·Flux
  │   ├─ Update: Q₁ = Qⁿ + α₁·dt·dQ
  │   └─ set_bc(), MPI_exchange() async
  │
  ├─ [Stage 2] RK coefficient α₂
  │   └─ (repeat for Q₂)
  │
  ├─ [Optional] SBLI Rescaling
  │   ├─ calc_mean(Q_rescaled)
  │   ├─ MPI_sync_mean()
  │   └─ step_rescale(Q_rescaled → Q)
  │
  └─ [Output] Every np steps
      └─ I/O to recal/ directory
```

### 6.2 MPI Rank Assignment

```
Rank 0  ──GPU 0──┐
Rank 1  ──CPU──  └─ Group 0
        ──GPU 0──┐
Rank 2  ──GPU──  │
Rank 3  ──CPU──  └─ Group 1
        ──GPU 1──┐
Rank 4  ──GPU──  │
Rank 5  ──CPU──  └─ Group 2

Pattern: mygpu = myrank / 2
  Even ranks compute
  Odd ranks gather for I/O
```

---

## 7. Performance Considerations

### 7.1 GPU Occupancy

**Target**: > 50% occupancy for mem-bound kernels

**Tuning**:
- Thread block dimensions (via `threadsE`, `threadsF`, `threadsG`)
- Shared memory usage (<96 KB per block typical)
- Register pressure (monitor with `pgprof`)

### 7.2 Memory Bandwidth

**Bottleneck**: E, F, G flux arrays each 5×(nx,ny,nz) reads/writes

**Optimization**:
- Shared memory reduces by 10-100x
- Coalesced access pattern (CUDA C++ best practices)
- Use lower precision where acceptable

### 7.3 Load Balancing

**Challenge**: Different schemes have different costs

**Balance Method**:
- Equal grid size per rank (z-direction decomposition)
- Assume KEEP ≈ SLAU ≈ Roe (similar flops)
- LES adds ~20% overhead

---

## 8. Extension Points

### 8.1 Adding New Spatial Schemes

1. Create `calc_mynewscheme_kernel.f90`
2. Implement 2D/3D GPU kernels (x, y, z directions)
3. Register in `calc_flux_base.f90` dispatcher
4. Add mod_globals parameter selector

### 8.2 Adding New Boundary Conditions

1. Implement in `set_bc_mynewbc.f90`
2. Called after each RK stage in `calc_time_dev.f90`
3. Ensure 2nd/4th/6th-order variants if accuracy-dependent

### 8.3 Adding New Post-Processing

1. Modify `print.f90` output formatting
2. Call from checkpoint I/O loop in `main.f90`
3. Consider VTK XML format for ParaView compatibility

---

## 9. Debugging & Profiling

### 9.1 Memory Access Patterns

```bash
# Check for misaligned access
nsys profile --trace cuda,mpi ./a.out

# Kernel-level metrics
ncu --set full ./a.out
```

### 9.2 MPI Communication

```bash
# Timeline of send/recv
nsys profile --trace mpi ./a.out
```

### 9.3 Algorithm Correctness

```bash
# Check intermediate values
pgdbg ./a.out  # NVIDIA debugger
```

---

**Document Version**: 1.0 (April 5, 2026)  
**Last Updated**: System architecture comprehensive review
