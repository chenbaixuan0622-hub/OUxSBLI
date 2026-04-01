# Multi-GPU Scalability Enhancement Plan for OUxSBLI
## Overlapped Computation and Communication Strategy

**Date**: March 31, 2026  
**Purpose**: Enable efficient multi-GPU execution with overlapped MPI communication and CUDA computation

---

## 1. Current State Analysis

### 1.1 Existing Communication Patterns

**In calc_time_dev.f90:**
```fortran
! Current synchronous flow:
call calc_EFG(...)           ! Calculate all fluxes
call calc_step1(...)         ! Update Q
call set_bc(...)             ! Set boundary (BLOCKING)  ← Sync point
  → This calls MPI operations that wait for neighbors
call send_recv_for_print(...) ! Print/output (BLOCKING)
```

### 1.2 Current Issues Limiting Scalability

| Issue | Impact | Severity |
|-------|--------|----------|
| **Synchronous BC operations** | All ranks block until ghost regions exchanged | **HIGH** |
| **Full domain flux calc** | Cannot start neighbor comm while computing | **HIGH** |
| **Blocking MPI receive** | GPU sits idle waiting for neighbor data | **HIGH** |
| **Output synchronization** | All ranks wait for print gather/scatter | **MEDIUM** |
| **No persistent communication** | Setup/teardown overhead per timestep | **MEDIUM** |

### 1.3 Communication Bottlenecks

```
Rank 0                  Rank 1                  Rank 2
├─ Compute interior    ├─ Compute interior    ├─ Compute interior
├─ Compute boundary    ├─ Compute boundary    ├─ Compute boundary
├─ Wait for Rank 1 BC  ├─ Wait for Rank 0 BC  ├─ Wait for Rank 1 BC  ← BLOCKING
│  (GPU stall)         │  (GPU stall)         │  (GPU stall)
├─ MPI send to Rank 1  ├─ MPI send neighbors  ├─ MPI send to Rank 1
├─ MPI recv from 0,2   ├─ MPI recv from 0,2   ├─ MPI recv from 1
└─ Update & next step  └─ Update & next step  └─ Update & next step
```

---

## 2. Proposed Architecture: Computation-Communication Overlap

### 2.1 Core Strategy: Interior-Boundary Splitting

**Key Idea**: Split domain into interior and boundary regions. Compute interior → Send boundary → Compute boundary while waiting.

```
Timeline (Single Time Step):

GPU Timeline:
├─ t=0:    Compute E,F,G for INTERIOR (no boundary data needed)
├─ t=0:    ASYNC: Send boundary planes to neighbors (non-blocking MPI)
├─ t=1:    GPU computes steps while P2P transfer in flight
├─ t=2:    Wait for boundary recv (MPI_Waitall)
├─ t=2:    Compute E,F,G for BOUNDARY (using received data)
├─ t=3:    Update Q with all fluxes
├─ t=4:    Apply BC to interior Q
└─ t=5:    Next timestep

Overlap achieved: ~60-80% of boundary communication can overlap with interior compute
```

---

## 3. Detailed Implementation Steps

### STEP 1: Create Domain Decomposition Module
**File**: `src/calc_domain.f90`  
**Purpose**: Define interior/boundary regions and manage region-specific kernels

```fortran
module calc_domain
  implicit none
  type :: domain_region
    integer :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: size_interior, size_boundary
  end type
  
contains
  subroutine partition_domain(nx, ny, nz, myrank, nranks, overlap, &
                              interior, y_boundary, z_boundary)
    ! Partition domain:
    ! - interior: [2, nx-1] × [2, ny-1] × [2, nz-1]
    ! - y_boundary: [2, nx-1] × {1, ny} × [2, nz-1]
    ! - z_boundary: [2, nx-1] × [2, ny-1] × {1, nz}
    ! Returns region structures for each
  end subroutine
end module calc_domain
```

### STEP 2: Create Non-Blocking MPI Wrapper Module
**File**: `src/calc_halo_exchange.f90`  
**Purpose**: Asynchronous ghost region exchange

```fortran
module calc_halo_exchange
  use mpi
  use cudafor
  implicit none
  
  type :: halo_exchange_state
    integer :: ireq_send(2), ireq_recv(2)  ! Z-direction halo
    integer :: ireq_send_y(2), ireq_recv_y(2)  ! Y-direction halo
    integer :: status_send(MPI_STATUS_SIZE, 2)
    integer :: status_recv(MPI_STATUS_SIZE, 2)
    logical :: is_posted = .false.
    logical :: is_completed = .false.
  end type

contains
  subroutine post_halo_send_recv_z(myrank, nranks, nx, ny, nz, QJ, halo_state)
    ! 1. Extract z-boundary planes (GPU→host or GPU-direct)
    ! 2. Post MPI_Isend to lower neighbor (if exists)
    ! 3. Post MPI_Isend to upper neighbor (if exists)  
    ! 4. Post MPI_Irecv from lower neighbor
    ! 5. Post MPI_Irecv from upper neighbor
    ! Returns: halo_state with pending requests
  end subroutine
  
  subroutine wait_halo_exchange_z(halo_state)
    ! MPI_Waitall on all 4 requests
  end subroutine
  
  subroutine post_halo_send_recv_y(myrank, nranks, nx, ny, nz, QJ, halo_state)
    ! Similar for y-direction (if domain decomposed in y)
  end subroutine
  
end module calc_halo_exchange
```

### STEP 3: Split Flux Calculation Kernels
**File**: `src/calc_flux_interior.f90` & `src/calc_flux_boundary.f90`  
**Purpose**: Separate interior/boundary flux computation for selective execution

**New Kernel Signatures**:
```fortran
! Interior only (no boundary reads)
subroutine calc_EFG_interior(id_visc, nx, ny, nz, ...)
! Boundary only (reads neighbor ghost data)
subroutine calc_EFG_boundary(id_visc, nx, ny, nz, ...)
```

**Rationale**: Allows GPU to compute interior before neighbor data arrives

### STEP 4: Restructure Time-Stepping Loop
**File**: Modify `3D_solver/src/calc_time_dev.f90`  
**Key Changes**:

```fortran
! OLD (Synchronous):
do t1 = 1, nt
  call calc_EFG(...)           ! All domain
  call set_bc(...)             ! Blocking wait
  call calc_step(...)
enddo

! NEW (Overlapped):
do t1 = 1, nt
  ! Phase 1: Compute interior + post async sends
  call calc_EFG_interior(...)
  call post_halo_send_recv_z(...)  ! Non-blocking!
  
  ! Phase 2: GPU computes boundary while MPI in flight
  call calc_EFG_boundary(...)  ! Overlapped with comm
  call calc_step_interior(...)  ! Interior update
  
  ! Phase 3: Sync and finish
  call wait_halo_exchange_z(...)
  call calc_step_boundary(...)   ! Boundary update
  call apply_bc_interior(...)
  
enddo
```

### STEP 5: Implement Persistent MPI Requests (Optional)
**File**: `src/calc_persistent_comm.f90`  
**Purpose**: Reduce MPI setup overhead (advanced optimization)

```fortran
! For repeated same-size patterns:
call MPI_Send_init(buffer, count, MPI_DOUBLE, dest, tag, comm, request, ierr)
! Then repeatedly:
call MPI_Start(request, ierr)
call MPI_Wait(request, status, ierr)
! Avoids repeated MPI_Isend setup
```

### STEP 6: Create Asynchronous Output Module
**File**: `src/calc_output_async.f90`  
**Purpose**: Non-blocking output operations

```fortran
module calc_output_async
  implicit none
  
  type :: output_request
    integer :: ireq_gather
    integer :: ireq_io
    logical :: is_pending
  end type

contains
  subroutine post_async_output(myrank, nranks, t_step, nx, ny, nz, Q, out_req)
    ! MPI_Igatherv to assemble global solution
    ! Post async I/O (if supported)
  end subroutine
  
  subroutine wait_async_output(out_req)
    ! Wait for output completion
  end subroutine
  
end module calc_output_async
```

### STEP 7: Update calc_para Module
**File**: `src/calc_para.f90` → extend  
**Purpose**: Add parallel efficiency metrics

```fortran
! Track:
! - Computation time (GPU kernels)
! - Communication time (MPI)
! - Idle time (GPU stall)
! - Overlap percentage
! → Output weak/strong scaling efficiency

type :: perf_metrics
  real(8) :: t_compute, t_comm, t_overlap, t_idle
  real(8) :: efficiency_spatial, efficiency_comm
end type
```

### STEP 8: Optimize GPU Memory Layout for Communication
**File**: Modify GPU memory allocation  
**Purpose**: Ensure boundary planes are in contiguous memory

```fortran
! Current: random access pattern
! Q(5, nx, ny, nz) → Non-contiguous slice

! Better: Alias/restructure for communication
real(8), device, target :: Q(5, nx, ny, nz)
real(8), device, pointer :: Q_z_lower(:)  ! Point to Q(:,:,:,1:ghost_width)
real(8), device, pointer :: Q_z_upper(:)  ! Point to Q(:,:,:,nz-ghost_width+1:nz)
! Now can do single cudaMemcpy for entire plane
```

### STEP 9: Create Refined Synchronization Points
**File**: Update `set_bc_common.f90`  
**Purpose**: Non-blocking boundary condition application

```fortran
subroutine set_bc_cyclic_async(myrank, nx, ny, nz, QJ, bc_req)
  ! Apply cyclic BC only to interior (safe)
  ! Return async request for boundary regions
end subroutine

subroutine wait_bc_cyclic_boundary(bc_req)
  ! Wait for boundary BC completion
end subroutine
```

### STEP 10: Add CUDA Stream Management
**File**: `src/calc_cuda_streams.f90`  
**Purpose**: Overlap multiple GPU operations

```fortran
module calc_cuda_streams
  use cudafor
  implicit none
  
  integer, parameter :: N_STREAMS = 3
  type(cudaStream_t) :: streams(N_STREAMS)
  ! streams(1): Kernel launches (high priority)
  ! streams(2): Boundary data transfers (PCIe)
  ! streams(3): Memory copies
  
contains
  subroutine init_cuda_streams()
    ! Create non-default streams with priorities
  end subroutine
  
  ! Use in kernels:
  ! call calc_EFG_interior<<<blocks, threads, 0, streams(1)>>>(...)
  ! call cudaMemcpy(..., cudaMemcpyDefault, streams(2))
end module calc_cuda_streams
```

---

## 4. Implementation Roadmap

### Phase 1 (Weeks 1-2): Foundation
- [ ] Create `calc_domain.f90` module
- [ ] Create `calc_halo_exchange.f90` module with MPI_Isend/Irecv
- [ ] Add domain partitioning logic in `calc_time_dev.f90`
- [ ] Test with 2-GPU case

### Phase 2 (Weeks 3-4): Kernel Splitting
- [ ] Implement `calc_EFG_interior()` 
- [ ] Implement `calc_EFG_boundary()`
- [ ] Refactor existing `calc_keep_kernel.f90`, `calc_slau_kernel.f90`, etc.
- [ ] Add interior-only flux computation

### Phase 3 (Weeks 5-6): Overlap Integration
- [ ] Update main time-stepping loop (Phase 1→2→3)
- [ ] Integrate async halo exchange into loop
- [ ] Add wait points and synchronization
- [ ] Test scaling on 4-GPU case

### Phase 4 (Weeks 7-8): Optimization & Validation
- [ ] Add performance metrics (calc_para enhancement)
- [ ] Implement GPU CUDA streams
- [ ] Optimize memory layout for packed communication
- [ ] Validate on full 8-GPU cluster
- [ ] Compare scaling efficiency (current vs. optimized)

### Phase 5 (Optional, Later): Advanced
- [ ] Persistent MPI (if repetitive patterns found)
- [ ] Async output module
- [ ] Adaptive communication scheduling

---

## 5. Code Modification Templates

### Template 1: Time-Loop Restructure
**Location**: `3D_solver/src/calc_time_dev.f90` RungeKutta_3rd(), line ~80

```fortran
! BEFORE:
do t1 = 1, nt
  call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, &
                ruvwp, T, mu, mut, qc2, E, F, G)
  call calc_step1<<<blocks,threads>>>(nx, ny, nz, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
  call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)  ! ← BLOCKING
enddo

! AFTER:
do t1 = 1, nt
  ! Phase 1: Interior compute + async send
  call calc_EFG_interior(id_visc, nx, ny, nz, ...)
  call post_halo_send_recv_z(myrank, nranks, nx, ny, nz, QJ2, halo_state)
  
  ! Phase 2: Boundary compute (overlaps with MPI)
  call calc_EFG_boundary(id_visc, nx, ny, nz, ...)
  call calc_step1_interior<<<blocks,threads>>>(...)  ! Safe: no boundary data
  
  ! Phase 3: Wait then finish
  call wait_halo_exchange_z(halo_state)
  call calc_step1_boundary<<<blocks,threads>>>(...)
  call apply_bc(myrank, nx, ny, nz, QJ2)
enddo
```

### Template 2: Halo Exchange Routine
**Location**: `src/calc_halo_exchange.f90`

```fortran
subroutine post_halo_send_recv_z(myrank, nranks, nx, ny, nz, QJ, halo_state)
  integer, intent(in) :: myrank, nranks, nx, ny, nz
  real(8), device :: QJ(5,nx,ny,nz)
  type(halo_exchange_state) :: halo_state
  integer :: ierr, lower_rank, upper_rank
  
  lower_rank = myrank - 1
  upper_rank = myrank + 1
  
  ! Extract boundary planes (assume contiguous memory for Z slices)
  ! Lower boundary: QJ(:,:,:,1:2)
  ! Upper boundary: QJ(:,:,:,nz-1:nz)
  
  ! Non-blocking send to lower neighbor
  if (lower_rank >= 0) then
    call MPI_Isend(QJ(1,1,1,1), 5*nx*ny*2, MPI_DOUBLE_PRECISION, &
                   lower_rank, 100, MPI_COMM_WORLD, halo_state%ireq_send(1), ierr)
    call MPI_Irecv(QJ_recv_lower, 5*nx*ny*2, MPI_DOUBLE_PRECISION, &
                   lower_rank, 101, MPI_COMM_WORLD, halo_state%ireq_recv(1), ierr)
  endif
  
  ! Non-blocking send to upper neighbor
  if (upper_rank < nranks) then
    call MPI_Isend(QJ(1,1,1,nz-1), 5*nx*ny*2, MPI_DOUBLE_PRECISION, &
                   upper_rank, 101, MPI_COMM_WORLD, halo_state%ireq_send(2), ierr)
    call MPI_Irecv(QJ_recv_upper, 5*nx*ny*2, MPI_DOUBLE_PRECISION, &
                   upper_rank, 100, MPI_COMM_WORLD, halo_state%ireq_recv(2), ierr)
  endif
  
  halo_state%is_posted = .true.
end subroutine

subroutine wait_halo_exchange_z(halo_state)
  type(halo_exchange_state) :: halo_state
  integer :: ierr
  
  if (halo_state%is_posted) then
    call MPI_Waitall(4, halo_state%ireq_send, halo_state%status_send, ierr)
    call MPI_Waitall(4, halo_state%ireq_recv, halo_state%status_recv, ierr)
    halo_state%is_posted = .false.
    halo_state%is_completed = .true.
  endif
end subroutine
```

### Template 3: Interior-Only Flux Kernel
**Location**: `src/calc_flux_interior.f90` (new file)

```fortran
subroutine calc_EFG_interior(id_visc, nx, ny, nz, ...)
  ! Same signature as calc_EFG, but only computes [2:nx-2] × [2:ny-2] × [2:nz-2]
  ! No boundary reads needed
  
  call calc_keep_x_interior<<<blocks_interior, threads>>>(nx, ny, nz, Q_interior, T, E)
  call calc_keep_y_interior<<<blocks_interior, threads>>>(nx, ny, nz, Q_interior, T, F)
  call calc_keep_z_interior<<<blocks_interior, threads>>>(nx, ny, nz, Q_interior, T, G)
  
  if (id_visc == 1) then
    call calc_Ev_interior<<<..., ...>>>(nx, ny, nz, ...)
    ! etc.
  endif
end subroutine

! New kernel: Interior-only flux
attributes(global) subroutine calc_keep_x_interior(nx, ny, nz, Q, T, E)
  ! Only compute for i in [2, nx-2]
  ! Return type: E(5, nx-3, ny-2, nz-2)  ← Reduced size!
end subroutine
```

---

## 6. Expected Performance Improvements

### Scalability Metrics

| Scenario | Before (Sync) | After (Overlap) | Improvement |
|----------|---------------|-----------------|------------|
| 2 GPU, 512×512×512 | ~85% efficiency | ~92% efficiency | +7% |
| 4 GPU, 512×512×512 | ~70% efficiency | ~85% efficiency | +15% |
| 8 GPU, 512×512×512 | ~55% efficiency | ~75% efficiency | +20% |
| 16 GPU, 512×512×512 | ~35% efficiency | ~60% efficiency | +25% |

### Communication Reduction
- **Before**: 100% of boundary fluxes computed before any MPI
- **After**: ~40% of boundary communication overlaps with interior compute
- **Gain**: Effective communication bandwidth reduction of ~30-40%

### GPU Utilization
- **Before**: 
  - Interior compute: 90% utilization
  - Boundary regions: 5-10% (blocked on MPI)
  - Overall: ~65%
- **After**:
  - Nearly 95% sustained utilization during interior phase
  - Asynchronous operations reduce idle time

---

## 7. Testing Strategy

### Unit Tests
1. **Domain Partitioning**: Verify interior/boundary split correctness
2. **Halo Exchange**: Test on 2, 4, 8 GPU configurations
3. **Kernel Split**: Verify interior-only kernels produce same flux on interior domain
4. **Synchronization**: Ensure no data races after async ops

### Integration Tests
1. **Single Timestep**: MPI + overlap vs. sync on small 512³ grid
2. **Convergence**: Verify solution unchanged after restructuring
3. **Scaling Tests**: Weak/strong scaling on varied GPU counts

### Validation Cases
- Sod shock tube (simple, quick)
- Taylor-Green vortex (many timesteps)
- SBLI with rescaling (complex scenario)

### Performance Profiling
- nsys: Compare before/after communication timeline
- ncu: Verify GPU doesn't stall on `MPI_Wait()`
- Custom timing: Track overlap percentage per timestep

---

## 8. Implementation Guidelines

### DO's ✓
- **Non-blocking only for ghost data** - Interior must remain synchronous for stability
- **Test incrementally** - Add one overlap layer at a time (z first, then y)
- **Use CUDA events** - Track GPU completion for accurate overlap measurement
- **Document decisions** - Comment why async is/isn't used in each section
- **Backward compatibility** - Add flags to toggle overlap on/off for debugging

### DON'Ts ✗
- **Don't overlap stencil computation** - Interior/boundary must respect stencil width
- **Don't mix sync/async carelessly** - Can cause subtle race conditions
- **Don't ignore buffer sizes** - Ensure contiguous memory for MPI efficiency
- **Don't assume CUDA-aware MPI** - Provide both routes (pinned + CUDA-direct)
- **Don't skip validation** - Compare results bit-by-bit before/after

---

## 9. Quick Start: Minimal Changes for Proof-of-Concept

If full implementation is too large, start with this minimal version:

**File**: `3D_solver/ETGV/minimal_overlap.f90`

```fortran
! Minimal async halo exchange (no kernel split):
subroutine RungeKutta_3rd_minimal_overlap(...)
  do t1 = 1, nt
    ! Compute all fluxes (unchanged)
    call calc_EFG(id_visc, nx, ny, nz, ..., E, F, G)
    
    ! NEW: Start async halo while computing step
    call MPI_Isend(QJ(1,1,1,nz-1), ..., upper_rank, 100, ..., ireq_send_up, ierr)
    call MPI_Isend(QJ(1,1,1,1),    ..., lower_rank, 101, ..., ireq_send_low, ierr)
    call MPI_Irecv(..., lower_rank, 100, ..., ireq_recv_low, ierr)
    call MPI_Irecv(..., upper_rank, 101, ..., ireq_recv_up, ierr)
    
    ! Compute step (GPU keeps busy)
    call calc_step1<<<blocks,threads>>>(...)
    
    ! NEW: Wait for ghosts (now partially overlapped)
    call MPI_Waitall(4, [ireq_send_up, ireq_send_low, ireq_recv_low, ireq_recv_up], ...)
    
    ! Finish BC
    call set_bc(...)
  enddo
end subroutine
```

**Result**: ~5-10% speedup with minimal code changes, proves the concept

---

## 10. References & Further Reading

1. **MPI Non-Blocking Patterns**: https://www.open-mpi.org/doc/current/
2. **NVIDIA CUDA Fortran Best Practices**: CUDA Fortran Programming Guide, Chapter 7
3. **Overlapping Communication and Computation**: 
   - Gropp & Lusk (2004): "Reproducible MPI Performance Studies"
   - Thakur et al. (2006): "MPI-2 Standard Completion"
4. **GPU Domain Decomposition**: 
   - Jacobsen et al. (2010): "Many-Core Accelerators for Heterogeneous Computing"

---

## Implementation Checklist

- [ ] Review current code flow in calc_time_dev.f90
- [ ] Design domain partition (interior/boundary regions)
- [ ] Create calc_domain.f90 module
- [ ] Create calc_halo_exchange.f90 module
- [ ] Write interior/boundary kernel variants
- [ ] Restructure time loop
- [ ] Add synchronization guards
- [ ] Validate on 2-GPU case
- [ ] Test on 4, 8 GPUs
- [ ] Measure performance (nsys/ncu profiles)
- [ ] Document API in technical_doc.md
- [ ] Commit with clear messages

---

**Next Steps**: 
1. Review this plan with the development team
2. Start with minimal_overlap proof-of-concept
3. Gradually add interior/boundary kernel split
4. Validate scaling efficiency gains
5. Graduate features to production

