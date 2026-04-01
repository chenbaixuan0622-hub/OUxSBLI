# OUxSBLI Technical Documentation

## 1. Project Overview

**OUxSBLI** is a GPU-accelerated Computational Fluid Dynamics (CFD) solver written in **CUDA Fortran**. It solves the compressible Navier-Stokes and Euler equations on structured rectilinear grids using explicit high-order finite-difference schemes.

### Key Features
- **GPU Acceleration**: CUDA Fortran kernels run on NVIDIA GPUs
- **High-Order Schemes**: 2nd, 4th, 6th, 8th-order spatial discretization options
- **Multiple Schemes**: KEEP (Kinetic Energy Preserving), SLAU (Low Dissipation AUSM), Roe
- **Parallel Scalability**: MPI domain decomposition with GPU-aware communication
- **LES Support**: Large Eddy Simulation with selective mixed scale models
- **Reference State Control**: Special rescaling module for shock-turbulence interaction (SBLI) simulations

### Dependencies
- HPC SDK (version 24.* or 25.*)
- CUDA toolkit (compatible with HPC SDK)
- MPI library
- ParaView (VTK output visualization)
- Python interface (development): pybind11

---

## 2. Code Architecture

### Directory Structure

```
OUxSBLI/
├── 1D_solver/          # 1D solvers (under development)
│   ├── NS/             # Navier-Stokes 1D
│   ├── ST/             # Shock tube test case
│   └── src/            # Common 1D routines
├── 2D_solver/          # 2D solvers
│   ├── 2D/, DSL/, KHI/, NS/, SBLI/, etc.  # Various test cases
│   └── src/            # Common 2D source routines
├── 3D_solver/          # Main 3D solver (most complete)
│   ├── ETGV/           # Euler Taylor-Green Vortex
│   ├── IVST/           # Isentropic Vortex
│   ├── KHI/            # Kelvin-Helmholtz Instability
│   ├── NSTGV/          # Navier-Stokes Taylor-Green Vortex
│   ├── SBLI/           # Shock-Boundary Layer Interaction
│   ├── TBL/            # Turbulent Boundary Layer
│   └── src/            # Core 3D source routines
├── src/                # Top-level utilities
│   ├── main.f90        # Main entry point
│   ├── sbli.f90        # Alternative SBLI-specific main
│   ├── mod_globals_c.f90      # C bindings for Python interface
│   ├── mod_constant.f90       # Physical constants
│   ├── cpu_gpu_mpi.f90        # CPU-GPU MPI abstractions
│   ├── print.f90              # Output/printing utilities
│   ├── calc_muscl.f90         # MUSCL reconstruction
│   ├── calc_physical_quantities.f90  # Primitive-conservative conversion
│   ├── set_compressible_bl.f90      # Boundary layer setup
│   └── set_coordinate.f90           # Coordinate transformations & Jacobian
└── docs/
    ├── api.md          # Auto-generated API reference
    └── technical_doc.md # This file
```

### Configuration Per Case

Each case folder (3D_solver/ETGV, SBLI, etc.) contains:

```
case_name/
├── mod_globals.f90     # Global parameters (equation type, scheme, grid size, block dimensions)
├── set.f90             # Grid, initial conditions, boundary conditions
├── Makefile            # Build configuration
├── calc.sh             # Execution script
├── profile.sh          # Profiling script (when available)
└── data/               # Output directory
```

---

## 3. Spatial Discretization

### 3.1 Convection Schemes

#### KEEP (Kinetic Energy and Entropy Preserving)
- **File**: `calc_keep_kernel.f90` in src/
- **Order**: Up to 6th-order
- **Features**: 
  - Splits convective forms to preserve kinetic energy and entropy
  - High accuracy on smooth flows
  - Less dissipative than traditional schemes
- **Kernels**:
  - `calc_keep_x[2,4,6]()`: X-direction flux
  - `calc_keep_y[2,4,6]()`: Y-direction flux  
  - `calc_keep_z[2,4,6]()`: Z-direction flux

#### SLAU (Simple Low Dissipation AUSM)
- **File**: `calc_slau_kernel.f90` in src/
- **Order**: Up to 6th-order with shock sensor
- **Features**:
  - Low dissipation away from shocks
  - Shock-adaptive via sensor function
  - Smoother than Roe for compressible turbulence
- **Kernels**:
  - `calc_slau_x[2,4,6]()`: X-direction with sensor
  - `calc_slau_y[2,4,6]()`: Y-direction with sensor
  - `calc_slau_z[2,4,6]()`: Z-direction with sensor

#### Roe Scheme
- **File**: `calc_roe_kernel.f90` in src/
- **Features**:
  - Classical Roe flux splitting
  - Shock-adaptive with sensor
  - Good for supersonic/hypersonic flows

#### Hybrid Schemes
- **File**: `calc_hybrid_kernel.f90`
- **Features**: Blends schemes based on flow conditions

### 3.2 Viscous Terms

#### 2nd-Order Scheme
- **File**: `calc_visc2.f90`
- **Kernels**:
  - `calc_Ev2()`: X-direction viscous flux
  - `calc_Fv2()`: Y-direction viscous flux
  - `calc_Gv2()`: Z-direction viscous flux
  - `calc_Ev_LES2()`, `calc_Fv_LES2()`, `calc_Gv_LES2()`: With LES turbulence model

#### 4th-Order Scheme (Sandham's Laplacian Form)
- **File**: `calc_visc4.f90`
- **Features**: Recently added, validation in progress
- **Kernels**:
  - `calc_Ev4()`, `calc_Fv4()`, `calc_Gv4()`: Without LES
  - `calc_Ev_LES4()`, `calc_Fv_LES4()`, `calc_Gv_LES4()`: With LES

#### Viscous Boundary Conditions
- **File**: `set_bc_common.f90`
- Cyclic boundary conditions with multiple accuracy orders (2nd, 4th, 6th)
- LES mutation synchronization: `set_bc_mut_common()`

### 3.3 Turbulence Models

#### LES (Large Eddy Simulation)
- **File**: `calc_les.f90`
- **SGS Model**: Selective mixed scale model (under development)
- Added through `mut` (turbulent viscosity) and `qc2` (filtered kinetic energy) fields

---

## 4. Temporal Discretization

### Time Integration Methods

Implemented in `calc_time_dev.f90`:

#### Explicit Runge-Kutta 3rd-Order (3-3 TVD)
```fortran
subroutine RungeKutta_3rd(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)
```
- Three-stage TVD Runge-Kutta
- Used for stability with high-order spatial schemes

#### Explicit Runge-Kutta 4th-Order (4-4)
```fortran
subroutine RungeKutta_4th(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)
```
- Four-stage classical RK4
- Higher accuracy temporal integration

#### Rescaled Variants
- `RungeKutta_3rd_rescale()`: 3rd-order with reference state rescaling
- `RungeKutta_4th_rescale()`: 4th-order with reference state rescaling
- Used for shock-boundary layer interaction (SBLI) cases

#### Gauss-Legendre Runge-Kutta
- Under development

---

## 5. CUDA Fortran GPU Implementation

### 5.1 GPU Kernel Patterns

#### Thread Organization
- **3D Thread Blocks**: Each kernel uses 3D thread blocks for multi-dimensional stencils
  - E-flux kernels: `threadsE = (x, y, z)` 
  - F-flux kernels: `threadsF = (x, y, z)`
  - G-flux kernels: `threadsG = (x, y, z)`
  - Viscous kernels: Similar structure with different tile naming

- **Shared Memory**: 
  - Loaded per thread block to avoid global memory reads in inner loops
  - Padded with extra elements based on stencil width
  - Example for 6th-order: `sx = threadsE%x + 5` (3 points on each side + center)

#### Kernel Structure Example (calc_keep_x6)

```fortran
attributes(global) subroutine calc_keep_x6(id_accuracy, nx, ny, nz, Q, T, E)
  ! Declarations: thread indices (it, jt, kt), global coords (i, j, k)
  integer, parameter :: sx = threadsE%x + 5  ! Tile size with padding
  real(8), dimension(-1:sx*sy*sz-2), shared :: rho, u, v, w, p, tmp
  
  ! Cooperative loading into shared memory
  do ii = it-2, threadsE%x+3, blockDim%x
    i = i_base + ii
    if (within_bounds) then
      idx = ii + offset_yz
      rho(idx) = Q(1,i,j,k)     ! Load with padding
      ...
    endif
  enddo
  call syncthreads()  ! Synchronize all threads
  
  ! Local flux calculation using 6-point stencil
  if (3 <= i .and. i <= nx-3) then
    E(:,i,j-1,k-1) = KEEP6(rho(idx-2:idx+3), u(...), ...)
  elseif (2 <= i .and. i <= nx-2) then
    E(:,i,j-1,k-1) = KEEP4(rho(idx-1:idx+2), ...)  ! Reduce to 4-pt at bounds
  else
    E(:,i,j-1,k-1) = KEEP2(rho(idx:idx+1), ...)    ! 2-pt at edges
  endif
end subroutine
```

#### Key Optimizations
- **Shared Memory**: Reduces global memory bandwidth by 10-100x
- **Coalesced Access**: Ensures memory transactions are efficient
- **Boundary Fallback**: Seamless degradation to lower-order schemes
- **Register Usage**: Carefully managed to maximize occupancy

### 5.2 Device vs. Host Memory

#### State Vector Q
- **Host (CPU)**: `real(8), allocatable :: Q(5,nx,ny,nz)`
- **Device (GPU)**: Allocated implicitly in kernels via `intent(in/out), device`
- **Transfer**: Between MPI ranks (see cpu_gpu_mpi module)

#### Intermediate Fields
- **E, F, G Fluxes**: Device arrays, computed in kernels
- **T (Temperature)**: Device array, computed from Q
- **mu (Dynamic Viscosity)**: Device array
- **mut, qc2**: Device arrays for LES turbulence

#### Persistent Device Data
- Coordinate derivatives: `dx, dy, dz` on device
- Metrics: `xix, etay, zetaz` (coordinate transformation)
- Jacobian: `Jacobian(nx,ny)` on device

---

## 6. MPI Parallelization

### 6.1 Domain Decomposition

- **Decomposition Strategy**: One-dimensional domain split along z-direction in 3D
- **Rank Assignment**: 
  - 2 MPI ranks per GPU: `mygpu = myrank / 2`
  - Even ranks (0, 2, 4, ...): Compute on GPU
  - Odd ranks (1, 3, 5, ...): Data gathering for output

### 6.2 CPU-GPU MPI Communication (cpu_gpu_mpi module)

Provides abstraction for MPI with GPU memory:

```fortran
interface CPUGPU_MPI_SEND
  module procedure CPU_MPI_SEND, GPU_MPI_SEND
end interface

interface CPUGPU_MPI_ISEND
  module procedure CPU_MPI_ISEND, GPU_MPI_ISEND
end interface
```

#### CPU_MPI_* (Pinned Host Memory)
- Copies device → host memory
- Performs MPI operation
- Requires extra memory copies but works on all systems

#### GPU_MPI_* (CUDA-Aware MPI)
- Direct MPI from GPU memory
- Requires CUDA-aware MPI library (NVIDIA OpenMPI)
- Lower latency and higher bandwidth

#### Non-Blocking Operations
- `CPU_MPI_ISEND/IRECV`: Asynchronous operations
- `GPU_MPI_ISEND/IRECV`: GPU-direct async
- Allow overlap of communication and computation

### 6.3 Synchronization

Boundary data exchange happens at each time step:
1. Compute interior fluxes
2. Async send/recv ghost regions
3. Wait for communication to complete
4. Update boundary states

---

## 7. State Representation

### Conservative Variable Vector Q

```
Q(1, i, j, k) = ρ      (density)
Q(2, i, j, k) = ρu     (x-momentum)
Q(3, i, j, k) = ρv     (y-momentum)
Q(4, i, j, k) = ρw     (z-momentum)
Q(5, i, j, k) = E      (total energy) or p (pressure, configuration-dependent)
```

### Flux Outputs

**E-flux** (x-direction):
- Shape: `(5, nx-1, ny-2, nz-2)` after stencil computation
- Index bounds: x in [2, nx-2], y in [2, ny-1], z in [2, nz-1]

**F-flux** (y-direction):
- Shape: `(5, nx-2, ny-1, nz-2)` after stencil computation
- Index bounds: x in [2, nx-1], y in [2, ny-2], z in [2, nz-1]

**G-flux** (z-direction):
- Shape: `(5, nx-2, ny-2, nz-1)` after stencil computation
- Index bounds: x in [2, nx-1], y in [2, ny-1], z in [2, nz-2]

### Derived Fields

**T (Temperature)**: `real(8), device :: T(nx,ny,nz)`
- Computed from Q using equation of state

**mu (Dynamic Viscosity)**: `real(8), device :: mu(nx,ny,nz)`
- Computed from T using Sutherland's law or similar

**mut, qc2 (LES)**: `real(8), device :: mut(nx,ny,nz), qc2(nx,ny,nz)`
- Turbulent viscosity and filtered kinetic energy
- Used in LES viscous term calculations

---

## 8. Core Modules Reference

### 8.1 mod_globals (Parameter Module)

**Typical Contents** (in each case's mod_globals.f90):

```fortran
implicit none
! Grid parameters
integer, parameter :: nx = 256, ny = 128, nz = 256
real(8), parameter :: Lx = 2.d0*pi, Ly = 1.d0*pi, Lz = 2.d0*pi

! Equation type: 0=Euler, 1=NS
integer(kind=2), parameter :: id_visc = 1

! Scheme selection
integer(kind=2), parameter :: id_scheme   = 1      ! 1=KEEP, 2=SLAU, 3=Roe
integer(kind=2), parameter :: id_accuracy = 8      ! 2,4,6,8 (points)
integer(kind=2), parameter :: id_tvd      = 2      ! Time stepping type
integer(kind=2), parameter :: id_slau     = ?      ! SLAU variant

! GPU thread/block configuration
type(dim3) :: threads, threadsE, threadsF, threadsG
type(dim3) :: blocks, blocksE, blocksF, blocksG
! Similar for viscous: threadsEv, threadsFv, threadsGv, etc.

! Rescaling (SBLI cases)
integer(kind=4), parameter :: id_rescale = 10      ! Rescaling interval or 0

! Restart flag
integer(kind=2), parameter :: id_recal = 2         ! 2=fresh, 4=restart
```

### 8.2 set.f90 (Problem Setup)

**Typical Implementation**:

```fortran
module set
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    ! Create coordinate arrays and grid spacings
  end subroutine

  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    ! Set initial conditions
    ! Q(:,:,:,:) = IC at each grid point
  end subroutine

  subroutine set_Jacobian_xy(nx, ny, nz, dx, dy, dz, Jacobian)
    ! Compute Jacobian determinant for coordinate transformation
  end subroutine
end module set
```

### 8.3 calc_flux_base.f90 (Flux Wrapper)

Orchestrates flux computation:

```fortran
subroutine calc_EFG_Euler(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, ...)
  ! Dispatch to appropriate convection scheme (KEEP/SLAU/Roe)
  ! Convert Euler fluxes to physical fluxes
end subroutine

subroutine calc_EFG_visc(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, ...)
  ! Compute viscous fluxes and add to total
end subroutine

subroutine calc_EFG_LES(id_visc, nx, ny, nz, dx, dy, dz, ...)
  ! LES variant with turbulent dissipation
end subroutine
```

### 8.4 set_bc_common.f90 (Boundary Conditions)

**Features**:
- **Cyclic BC** (periodic): `set_bc_cyclic2/4/6()`
- **Accuracy Variants**: 2-point, 4-point, 6-point stencils for interior copies
- **Kernel Variants**: CPU and GPU versions
- **LES Mutation**: Synchronize `mut` and `qc2` at boundaries

**Flow**:
```
set_bc_cyclic2_init(Q) -> CPU version for initial condition
set_bc_cyclic2(Q)      -> GPU kernel during time stepping
```

### 8.5 calc_time_dev.f90 (Time Stepping)

**Main Interface**:

```fortran
subroutine RungeKutta(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)
  ! Dispatch to specific RK scheme
  if (id_RungeKutta == 2) then
    call RungeKutta_3rd(...)
  else if (id_RungeKutta == 4) then
    call RungeKutta_4th(...)
  endif
end subroutine

subroutine RungeKutta_3rd(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, ...)
  ! Three-stage TVD RK3:
  ! Q1 = Q^n + Δt * RHS(Q^n)
  ! Q2 = (3*Q^n + Q1 + Δt*RHS(Q1)) / 4
  ! Q^{n+1} = (Q^n + 2*Q2 + 2*Δt*RHS(Q2)) / 3
  
  ! Within each stage:
  do step = 1, nt
    call calc_physical_quantities(Q, T, mu)           ! Temp, viscosity
    call calc_EFG_Euler(..., E, F, G)                ! Convection fluxes
    call calc_EFG_visc(..., E, F, G)                 ! Add viscous
    call compute_RHS(...)                            ! Finite-diff div(flux)
    call UpdateQ(...)                                ! Update Q
    call set_bc_cyclic(Q)                            ! Boundary conditions
    call MPI_exchange(Q)                             ! Ghost region sync
    if (mod(step, id_rescale) == 0) then
      call step_rescale(...)                         ! Rescale if needed
    endif
  enddo
end subroutine
```

### 8.6 set_coordinate.f90 (Metrics)

**Functions**:
```fortran
subroutine set_Jacobian_xy(nx, ny, nz, dx, dy, dz, Jacobian)
  ! For 2D/3D curvilinear grids: J = ∂(ξ,η)/∂(x,y)
end subroutine

subroutine set_coordinate_derivative(dx, dy, dz, xix, etay, zetaz)
  ! ∂ξ/∂x = 1/dx, ∂η/∂y = 1/dy, etc.
end subroutine
```

### 8.7 calc_physical_quantities.f90 (Properties)

**Functions**:
- Convert between conservative and primitive variables
- Compute temperature from internal energy
- Compute dynamic viscosity from temperature (Sutherland's law)
- Compute speed of sound

### 8.8 cpu_gpu_mpi.f90 (Communication)

**Interfaces**:
- `CPUGPU_MPI_SEND(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)`
- `CPUGPU_MPI_RECV(id_gpu_mpi, buf, count, src, tag, comm, ireq, ierr)`
- `CPUGPU_MPI_ISEND(...)`
- `CPUGPU_MPI_IRECV(...)`

**Dispatch Logic**:
- `id_gpu_mpi == 2`: Pinned memory route (CPU_MPI_*)
- `id_gpu_mpi == 4`: CUDA-aware route (GPU_MPI_*)

### 8.9 calc_rescale.f90 (Reference State)

Used for SBLI simulations to preserve background flow while evolving turbulence:

**Routines**:
- `calc_mean(step, ireq, flag_re, nx, ny, nz, Jacobian, QJ, Qm)`: Accumulate mean
- `step_rescale(...)`: Apply rescaling transformation
- `rescale_recv_send(...)`: MPI communication for mean data
- `set_rescale(flag_re, step, nx, ny, nz, y, Jacobian, Qm, bltre, Qre)`: Transform Q back

---

## 9. Main Execution Flow

### 9.1 Program Entry Point (src/main.f90)

```fortran
program main
  ! 1. MPI initialization
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  mygpu = myrank / 2
  
  ! 2. Setup GPU blocks/threads (even ranks only)
  if (mod(myrank,2) == 0) then
    call set_block3(nx, ny, nz, threads, threadsE, ..., blocks, blocksE, ...)
  endif
  
  ! 3. Allocate memory
  allocate(Q(5,nx,ny,nz), x(nx), dx(nx-1), y(ny), dy(ny-1), z(nz), dz(nz-1), Jacobian(nx,ny))
  
  ! 4. Setup grid
  call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
  call set_Jacobian_xy(nx, ny, nz, dx, dy, dz, Jacobian)
  
  ! 5. Initial condition (or restart)
  if (id_recal == 4) then
    ! Read from checkpoint file
    open(10, file="recal/Q" // formatted_rank // ".dat", ...)
    read(10) header, Q
  elseif (id_recal == 2) then
    call set_init(myrank, nx, ny, nz, x, y, z, Q)
  endif
  
  ! 6. Time stepping
  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, Q)
  call cpu_time(t_end)
  
  ! 7. Save output (even ranks)
  if (mod(myrank,2) == 0) then
    do l = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q(:,i,j,l) = Jacobian(i,j) * Q(:,i,j,l)    ! Apply Jacobian
        enddo
      enddo
    enddo
    write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
    open(10, file=filename, ...)
    write(10) 'SEQFMT01'
    write(10) Q
    close(10)
  endif
  
  call MPI_FINALIZE(ierr)
end program main
```

### 9.2 SBLI Variant (src/sbli.f90)

- Supports dual regions (boundary layer + shock interaction)
- Different grid sizes for each region: `nx1/ny1/nz1` vs `nx2/ny2/nz2`
- Ranks 0-1 handle region 1, ranks 2+ handle region 2
- Uses rescaling module for flux matching

---

## 10. Build & Execution Workflow

### 10.1 Compilation

```bash
cd 3D_solver/ETGV    # Navigate to case directory
vi mod_globals.f90   # Edit parameters (grid, scheme, accuracy, blocks/threads)
vi set.f90           # Edit grid, IC, BC
make                 # Compile with HPC SDK compiler
```

**Makefile typically uses**:
- `pgfortran` or `nvfortran` (HPC SDK compilers)
- CUDA Fortran compilation flags
- Link to MPI library

### 10.2 Execution

```bash
bash calc.sh                 # Run simulation
# Output appears in data/ and recal/ directories

# Optional profiling
bash profile.sh             # Generate nsys/ncu profiles
```

### 10.3 Checkpoint Files

- **Input**: `recal/Q00001.dat`, `recal/Q00002.dat`, ... (one per rank pair)
- **Output**: Same structure; overwritten each output interval
- **Format**: Sequential unformatted Fortran or stream binary
- **Data**: Q array with Jacobian applied

### 10.4 Performance Monitoring

- **nsys** (system profiler): GPU utilization, memory bandwidth
- **ncu** (kernel profiler): Per-kernel metrics (occupancy, bandwidth, etc.)
- Running `profile.sh` generates `.nsys-rep` and `.ncu-rep` files

---

## 11. Configuration Guide

### Parameter Selection

#### Equation Type (id_visc)
```
0 = Euler (inviscid)
1 = Navier-Stokes (viscous)
```

#### Convection Scheme (id_scheme)
```
1 = KEEP        (energy-preserving, recommended for smooth flows)
2 = SLAU        (low-dissipation with shock sensor)
3 = Roe         (classical, good for shocks)
```

#### Spatial Accuracy (id_accuracy)
```
2 = 2nd-order  (2-point stencil)
4 = 4th-order  (4-point stencil)
6 = 6th-order  (6-point stencil)
8 = 8th-order  (8-point stencil, not all schemes)
```

#### Time Integration (id_RungeKutta)
```
2 = 3rd-order TVD Runge-Kutta (recommended with high-order spatial)
4 = 4th-order Runge-Kutta
```

#### Rescaling (id_rescale)
```
0 = No rescaling
N > 0 = Rescale every N time steps (SBLI cases)
```

#### Restart (id_recal)
```
2 = Fresh start with set_init()
4 = Restart from checkpoint file
```

#### Block/Thread Tuning

**Example for NVIDIA A100**:
```fortran
type(dim3) :: threads = dim3(8, 16, 4)        ! 512 threads/block
type(dim3) :: threadsE = dim3(32, 8, 4)       ! E-flux layout
type(dim3) :: threadsF = dim3(8, 32, 4)       ! F-flux layout
type(dim3) :: threadsG = dim3(8, 8, 8)        ! G-flux layout
```

General heuristic:
- Total threads/block: 128-512 (avoid exceeding SM capacity)
- Balance x,y,z for stencil computation patterns
- More threads in direction with more memory reuse

---

## 12. Validation Cases

### 12.1 Euler Vortex Convection (1D_solver/NS, 3D_solver/ETGV)

- **Purpose**: Grid convergence study
- **Setup**: Isentropic vortex advecting in uniform flow
- **Metrics**: L² error vs grid resolution
- **Expected**: 2nd-order KEEP ≈ O(h²), 6th-order KEEP ≈ O(h⁶)

### 12.2 Sod Shock Tube (1D_solver/ST)

- **Purpose**: Shock capturing validation
- **Setup**: Riemann problem with discontinuity
- **Metrics**: Visual comparison with reference solution
- **Schemes**: SLAU, KEEP all capture shock with appropriate dissipation

### 12.3 3D Taylor-Green Vortex (3D_solver/ETGV, NSTGV)

- **Purpose**: High-order accuracy on smooth periodic flow
- **Setup**: Viscous decay on [0, 2π]³
- **Metrics**: Kinetic energy, enstrophy evolution
- **Expected**: Excellent agreement with analytical decay rates

### 12.4 M=1.9 Supersonic Boundary Layer (3D_solver/TBL)

- **Purpose**: Wall-bounded turbulent flow
- **Setup**: Turbulent boundary layer at Mach 1.9
- **Metrics**: Van Driest transformed mean velocity, Reynolds stress profiles
- **Validation**: Against DNS/experiments

### 12.5 Shock-Boundary Layer Interaction (3D_solver/SBLI)

- **Purpose**: Complex compressible turbulence
- **Setup**: Incident shock on supersonic boundary layer
- **Features**: Uses rescaling module to maintain reference state
- **Metrics**: Unsteady shock motion, separation bubble dynamics

---

## 13. Known Issues & Limitations

### Under Development
- 1D and 2D solvers: Feature-complete but less tested than 3D
- LES module (`calc_les.f90`): SGS model validation ongoing
- Python interface (`pyETGV`, `pyNSTGV`): Not yet functional
- Gauss-Legendre RK time stepping: Experimental

### Current Constraints
- **Grid**: Only rectilinear (no body-fitted curvilinear)
- **BC**: Primarily cyclic periodic; other BCs case-specific
- **Restart**: ASCII format only (slow I/O at high resolution)
- **Visualization**: Manual VTK output generation (not real-time)

### Performance Considerations
- **Memory**: Entire Q, E, F, G arrays in device memory; watch for OOM on small GPUs
- **Bandwidth**: Communication overhead increases with rank count; test scaling
- **Block/Thread Tuning**: May require optimization per GPU model and case

---

## 14. References

1. Kuya, Y., Totani, K., & Kawai, S. (2018). *Kinetic energy and entropy preserving schemes for compressible flows by split convective forms*. Journal of Computational Physics, 372, 359-386.

2. Kuya, Y., & Kawai, S. (2021). *High-order accurate kinetic-energy and entropy preserving (KEEP) schemes on curvilinear grids*. Journal of Computational Physics, 432, 110134.

3. Jacobs, C. T., Jammy, S. P., & Sandham, N. D. (2017). *OpenSBLI: A framework for the automated derivation and parallel execution of finite difference solvers on a range of computer architectures*. Computer Physics Communications, 220, 1-12.

4. Allaneau, Y., & Jameson, A. (2012). *Direct numerical simulations of a two-dimensional viscous flow in a shocktube using kinetic energy preserving scheme*. AIAA Paper 2009-3797.

5. Tamaki, Y., & Kawai, S. (2023). *Wall-modeled LES of transonic buffet over NASA-CRM using Cartesian-grid-based flow solver FFVHC-ACE*. AIAA Paper 2023-0429.

6. NVIDIA. (2017). *CUDA Fortran Programming Guide and Reference*. https://docs.nvidia.com/hpc-sdk/pgi-compilers/2017/pgi17cudaforug.pdf

---

## Appendix A: Quick Reference

### Environment Setup
```bash
module load hpc-sdk/24.x     # or appropriate HPC SDK module
source activate hpc_env      # if using conda
export CUDA_VISIBLE_DEVICES=0  # Select GPU
```

### Essential Commands
```bash
make clean && make           # Full rebuild
make -j4                     # Parallel compilation
bash calc.sh                 # Run with script settings
mpirun -np 4 ./a.out         # Direct MPI launch (4 ranks, 2 GPUs)
```

### Debugging
```bash
pgdbg ./a.out                # NVIDIA debugger
cuda-memcheck ./a.out        # Memory error detection (older CUDA)
compute-sanitizer ./a.out    # Memory error detection (newer CUDA)
```

### Profiling
```bash
nsys profile ./a.out         # System profiler
ncu ./a.out                  # Kernel profiler
ncu --config full ./a.out    # Full metrics
```

---

## Appendix B: File Type Conventions

| Extension | Purpose | Example |
|-----------|---------|---------|
| `.f90` | Fortran 90 source | `calc_keep_kernel.f90` |
| `.mod` | Compiled module (generated) | `calc_keep_kernel.mod` |
| `.dat` | Binary checkpoint data | `recal/Q00001.dat` |
| `.sh` | Bash script | `calc.sh`, `profile.sh` |
| `.md` | Markdown documentation | `api.md`, `technical_doc.md` |
| `Makefile` | Build automation | Configure compilation |

---

## Appendix C: Glossary

**Convection**: Hyperbolic part of Euler/NS equations (flux divergence)
**KEEP**: Kinetic Energy and Entropy Preserving scheme
**SLAU**: Simple Low Dissipation AUSM (approximate Riemann solver)
**TVD**: Total Variation Diminishing (stability property)
**LES**: Large Eddy Simulation (turbulence model)
**SGS**: Subgrid-Scale (modeled turbulence in LES)
**Jacobian**: Determinant of coordinate transformation; geometric factor
**Stencil**: Set of neighboring grid points used in finite-difference
**Shared Memory**: Fast on-chip CUDA memory, shared among thread block
**SBLI**: Shock-Boundary Layer Interaction
**ETGV**: Euler Taylor-Green Vortex
**TBL**: Turbulent Boundary Layer
**Rescaling**: Flux rescaling to maintain reference mean state

---

**Document Version**: 1.0 (March 31, 2026)  
**Generated From**: OUxSBLI source code analysis  
**Relevant Source Docs**: `api.md`, `README.md`
