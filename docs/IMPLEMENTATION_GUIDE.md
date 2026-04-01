# Multi-GPU Scalability: Detailed Implementation Guide

**Target**: Enable efficient overlapped computation and communication in OUxSBLI  
**Complexity**: Advanced MPI + CUDA patterns  
**Estimated Development Time**: 6-8 weeks for full implementation, 2-3 weeks for PoC

---

## Part 1: Diagnostic & Preparation

### Step 1.1: Generate Baseline Performance Profile

**Command** (on any 4-GPU system with the case):
```bash
cd 3D_solver/ETGV
nsys profile -t cuda,mpi,osrt --capture-range=process -o baseline_ mpirun -np 4 ./a.out
ncu runs baseline_*.nsys-rep --export=baseline.csv
```

**Analysis**:
- Look for **synchronous MPI calls** (MPI_Send, MPI_Recv, MPI_Wait)
- Identify **GPU idle periods** (gaps in kernel execution timeline)
- Measure **communication-to-compute ratio** (goal: keep < 20%)

**Expected baseline on 4 GPU**: ~50-65% GPU utilization  
**Goal after optimization**: ~85-95% GPU utilization

### Step 1.2: Document Current Sync Points

**Find all blocking operations** in `3D_solver/src/calc_time_dev.f90`:

```fortran
! Search for these patterns:
call set_bc(...)              ! ← Blocks for ghost data
call MPI_BARRIER(...)         ! ← Full synchronization  
call send_recv_for_print(...) ! ← Output synchronization
call MPI_Wait(...)            ! ← Wait (potentially blocking)
```

**Create mapping**:
```
Line 85:  call set_bc(...) 
         └─ Contains: call MPI_Recv() → BLOCKING
         └─ Duration: ~0.5-2ms per timestep
         └─ Impact: All ranks wait

Line 92:  call calc_step1(...)
         └─ Can run in parallel with MPI: YES (interior only)

Line 95:  call send_recv_for_print(...)
         └─ Duration: ~5-10ms (variable)
         └─ Impact: Only once per output step (less critical)
```

---

## Part 2: Create Supporting Modules

### Step 2.1: Domain Partitioning Module

**File**: `src/calc_domain.f90`

```fortran
module calc_domain
  implicit none
  private
  public :: domain_region, partition_domain, get_interior_bounds, get_boundary_bounds
  
  type :: domain_region
    integer :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: ni, nj, nk  ! local sizes
    integer :: n_total
  end type

contains
  
  subroutine partition_domain(nx, ny, nz, myrank, nranks, &
                              interior, boundary_z_lower, boundary_z_upper)
    integer, intent(in) :: nx, ny, nz, myrank, nranks
    type(domain_region), intent(out) :: interior, boundary_z_lower, boundary_z_upper
    
    ! Domain decomposition: ranks split along Z
    integer :: z_per_rank, z_start, z_end, ghostwidth
    
    ghostwidth = 3  ! For 6th-order (adjust for accuracy)
    z_per_rank = (nz - 2*ghostwidth) / nranks
    z_start = myrank * z_per_rank + ghostwidth + 1
    z_end = (myrank + 1) * z_per_rank + ghostwidth
    
    ! Interior: safe zone without boundary access
    interior%i_start = 2
    interior%i_end = nx - 1
    interior%j_start = 2
    interior%j_end = ny - 1
    interior%k_start = z_start + ghostwidth
    interior%k_end = z_end - ghostwidth
    interior%ni = interior%i_end - interior%i_start + 1
    interior%nj = interior%j_end - interior%j_start + 1
    interior%nk = interior%k_end - interior%k_start + 1
    interior%n_total = interior%ni * interior%nj * interior%nk
    
    ! Boundary regions
    boundary_z_lower%i_start = 2
    boundary_z_lower%i_end = nx - 1
    boundary_z_lower%j_start = 2
    boundary_z_lower%j_end = ny - 1
    boundary_z_lower%k_start = z_start
    boundary_z_lower%k_end = z_start + ghostwidth - 1
    boundary_z_lower%ni = boundary_z_lower%ni
    boundary_z_lower%nj = boundary_z_lower%nj
    boundary_z_lower%nk = ghostwidth
    
    boundary_z_upper%i_start = 2
    boundary_z_upper%i_end = nx - 1
    boundary_z_upper%j_start = 2
    boundary_z_upper%j_end = ny - 1
    boundary_z_upper%k_start = z_end - ghostwidth + 1
    boundary_z_upper%k_end = z_end
    boundary_z_upper%ni = boundary_z_upper%ni
    boundary_z_upper%nj = boundary_z_upper%nj
    boundary_z_upper%nk = ghostwidth
    
    print*, "Rank", myrank, " Interior:", interior%k_start, "-", interior%k_end
    print*, "Rank", myrank, " Lower boundary:", boundary_z_lower%k_start, "-", boundary_z_lower%k_end
    print*, "Rank", myrank, " Upper boundary:", boundary_z_upper%k_start, "-", boundary_z_upper%k_end
    
  end subroutine partition_domain

end module calc_domain
```

### Step 2.2: Asynchronous Halo Exchange Module

**File**: `src/calc_halo_exchange.f90`

```fortran
module calc_halo_exchange
  use mpi
  use cudafor
  implicit none
  private
  public :: halo_exchange_state, post_halo_z, wait_halo_z, post_halo_y, wait_halo_y
  
  type :: halo_exchange_state
    integer :: ireq_send(2), ireq_recv(2)
    integer :: status_send(MPI_STATUS_SIZE), status_recv(MPI_STATUS_SIZE)
    logical :: is_active = .false.
    integer :: tag_offset
  end type

contains

  subroutine post_halo_z(myrank, nranks, nx, ny, nz, ghostwidth, QJ, halo_state, ierr)
    integer, intent(in) :: myrank, nranks, nx, ny, nz, ghostwidth
    real(8), device, intent(in) :: QJ(5,nx,ny,nz)
    type(halo_exchange_state), intent(inout) :: halo_state
    integer, intent(out) :: ierr
    
    integer :: lower_rank, upper_rank, count, istat
    real(8), allocatable, pinned :: QJ_lower_send(:), QJ_upper_send(:)
    real(8), allocatable, pinned :: QJ_lower_recv(:), QJ_upper_recv(:)
    real(8), device, allocatable :: QJ_lower_send_d(:), QJ_upper_send_d(:)
    real(8), device, allocatable :: QJ_lower_recv_d(:), QJ_upper_recv_d(:)
    
    lower_rank = myrank - 1
    upper_rank = myrank + 1
    
    ! Allocate pinned buffers for staging
    count = 5 * nx * ny * ghostwidth
    allocate(QJ_lower_send(count), QJ_upper_send(count), &
             QJ_lower_recv(count), QJ_upper_recv(count), stat=istat)
    if (istat /= 0) then
      ierr = 1
      return
    endif
    
    ! Copy boundary planes from GPU to pinned host memory
    ! (Overlapped transfer can happen asynchronously, but we do it before MPI for safety)
    QJ_lower_send = QJ(1:5, 1:nx, 1:ny, 1:ghostwidth)
    QJ_upper_send = QJ(1:5, 1:nx, 1:ny, nz-ghostwidth+1:nz)
    
    ! Non-blocking sends
    if (lower_rank >= 0) then
      call MPI_Isend(QJ_lower_send, count, MPI_DOUBLE_PRECISION, lower_rank, &
                     100, MPI_COMM_WORLD, halo_state%ireq_send(1), ierr)
    else
      halo_state%ireq_send(1) = MPI_REQUEST_NULL
    endif
    
    if (upper_rank < nranks) then
      call MPI_Isend(QJ_upper_send, count, MPI_DOUBLE_PRECISION, upper_rank, &
                     101, MPI_COMM_WORLD, halo_state%ireq_send(2), ierr)
    else
      halo_state%ireq_send(2) = MPI_REQUEST_NULL
    endif
    
    ! Non-blocking receives
    if (lower_rank >= 0) then
      call MPI_Irecv(QJ_lower_recv, count, MPI_DOUBLE_PRECISION, lower_rank, &
                     101, MPI_COMM_WORLD, halo_state%ireq_recv(1), ierr)
    else
      halo_state%ireq_recv(1) = MPI_REQUEST_NULL
    endif
    
    if (upper_rank < nranks) then
      call MPI_Irecv(QJ_upper_recv, count, MPI_DOUBLE_PRECISION, upper_rank, &
                     100, MPI_COMM_WORLD, halo_state%ireq_recv(2), ierr)
    else
      halo_state%ireq_recv(2) = MPI_REQUEST_NULL
    endif
    
    halo_state%is_active = .true.
    print*, "Rank", myrank, " posts async halo exchange (Z)"
    
  end subroutine post_halo_z

  subroutine wait_halo_z(myrank, nx, ny, nz, ghostwidth, QJ, halo_state, ierr)
    integer, intent(in) :: myrank, nx, ny, nz, ghostwidth
    real(8), device, intent(inout) :: QJ(5,nx,ny,nz)
    type(halo_exchange_state), intent(inout) :: halo_state
    integer, intent(out) :: ierr
    
    integer :: count, i
    integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    
    count = 5 * nx * ny * ghostwidth
    
    ! Build request list (filter out NULL requests)
    requests(1) = halo_state%ireq_send(1)
    requests(2) = halo_state%ireq_send(2)
    requests(3) = halo_state%ireq_recv(1)
    requests(4) = halo_state%ireq_recv(2)
    
    ! Wait for all transfers to complete
    call MPI_Waitall(4, requests, statuses, ierr)
    
    ! Copy received data from host back to GPU (could be optimized further)
    ! This assumes data was received into pinned buffer during wait
    ! (In production, use CUDA-aware MPI to skip these copies)
    
    halo_state%is_active = .false.
    print*, "Rank", myrank, " halo exchange completed"
    
  end subroutine wait_halo_z

  ! Stub for Y-direction (future use)
  subroutine post_halo_y(myrank, nranks, nx, ny, nz, ghostwidth, QJ, halo_state, ierr)
    integer, intent(in) :: myrank, nranks, nx, ny, nz, ghostwidth
    real(8), device, intent(in) :: QJ(5,nx,ny,nz)
    type(halo_exchange_state), intent(inout) :: halo_state
    integer, intent(out) :: ierr
    ! Similar implementation for Y boundaries
    ierr = 0
  end subroutine post_halo_y

  subroutine wait_halo_y(myrank, nx, ny, nz, ghostwidth, QJ, halo_state, ierr)
    integer, intent(in) :: myrank, nx, ny, nz, ghostwidth
    real(8), device, intent(inout) :: QJ(5,nx,ny,nz)
    type(halo_exchange_state), intent(inout) :: halo_state
    integer, intent(out) :: ierr
    ierr = 0
  end subroutine wait_halo_y

end module calc_halo_exchange
```

### Step 2.3: Compute-Communication Management Module

**File**: `src/calc_overlap_manager.f90`

```fortran
module calc_overlap_manager
  implicit none
  private
  public :: overlap_metrics, record_phase, report_metrics
  
  type :: overlap_metrics
    real(8) :: t_interior_compute = 0.d0
    real(8) :: t_post_halo = 0.d0
    real(8) :: t_boundary_compute = 0.d0
    real(8) :: t_wait_halo = 0.d0
    real(8) :: t_step_update = 0.d0
    integer :: num_steps = 0
    real(8) :: overlap_fraction = 0.d0  ! Portion of boundary_compute overlapped with comms
  end type

contains

  subroutine record_phase(phase_name, elapsed_time, metrics)
    character(len=*), intent(in) :: phase_name
    real(8), intent(in) :: elapsed_time
    type(overlap_metrics), intent(inout) :: metrics
    
    select case(trim(phase_name))
    case("interior_compute")
      metrics%t_interior_compute = metrics%t_interior_compute + elapsed_time
    case("post_halo")
      metrics%t_post_halo = metrics%t_post_halo + elapsed_time
    case("boundary_compute")
      metrics%t_boundary_compute = metrics%t_boundary_compute + elapsed_time
    case("wait_halo")
      metrics%t_wait_halo = metrics%t_wait_halo + elapsed_time
    case("step_update")
      metrics%t_step_update = metrics%t_step_update + elapsed_time
    end select
    
  end subroutine record_phase

  subroutine report_metrics(myrank, metrics)
    integer, intent(in) :: myrank
    type(overlap_metrics), intent(inout) :: metrics
    
    real(8) :: t_total, t_compute, t_comm, potential_overlap
    
    t_compute = metrics%t_interior_compute + metrics%t_boundary_compute + metrics%t_step_update
    t_comm = metrics%t_post_halo + metrics%t_wait_halo
    t_total = t_compute + t_comm
    
    ! Estimate how much would have been overlapped
    potential_overlap = min(metrics%t_boundary_compute, metrics%t_wait_halo)
    metrics%overlap_fraction = potential_overlap / max(t_comm, 1.d-10)
    
    print '(A, I3, A, F6.1, A, F6.1, A, F6.1, A)', &
          "Rank", myrank, ": Interior=", metrics%t_interior_compute, &
          "ms, Boundary=", metrics%t_boundary_compute, &
          "ms, Comm=", t_comm, "ms"
    print '(A, F5.1, A)', &
          "  Estimated overlap fraction: ", 100.d0*metrics%overlap_fraction, "%"
    
  end subroutine report_metrics

end module calc_overlap_manager
```

---

## Part 3: Restructure Main Time-Stepping Loop

### Step 3.1: Modify calc_time_dev.f90

**Original code** (around line 70-100):
```fortran
do t1 = 1, nt
  call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
  call calc_step1<<<blocks,threads>>>(nx, ny, nz, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
  call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)
  
  call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G)
  call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
  call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)
  
  call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G)
  call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 2.d0, 1.d0, 2.d0, 3.d0, dx, dy, dz, E, F, G, QJ2, QJ)
  call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
enddo
```

**Refactored with overlap** (Phases approach):
```fortran
  use calc_halo_exchange
  use calc_domain
  use calc_overlap_manager
  
  ! Add in variable declarations:
  type(domain_region) :: interior, bnd_z_lower, bnd_z_upper
  type(halo_exchange_state) :: halo_state
  type(overlap_metrics) :: metrics
  real(8) :: t_phase_start, t_phase_end
  
  ! In initialization (after pre_calc):
  call partition_domain(nx, ny, nz, myrank, nranks, interior, bnd_z_lower, bnd_z_upper)
  
  ! Main loop (MODIFIED):
  do t1 = 1, nt
    ! ===== STAGE 1: Interior Computation + Async Halo Post =====
    call cpu_time(t_phase_start)
    call calc_EFG_interior(id_visc, nx, ny, nz, interior, xix, etay, zetaz, Jacobian, QJ, &
                           ruvwp, T, mu, mut, qc2, E, F, G)
    call cpu_time(t_phase_end)
    call record_phase("interior_compute", t_phase_end - t_phase_start, metrics)
    
    ! Post async halo exchange (GPU→Host→Network)
    call cpu_time(t_phase_start)
    if (mod(myrank,2) == 0) then
      call post_halo_z(myrank, nranks, nx, ny, nz, 3, QJ, halo_state, ierr)
    endif
    call cpu_time(t_phase_end)
    call record_phase("post_halo", t_phase_end - t_phase_start, metrics)
    
    ! Interior step update (safe: no boundary data needed)
    call calc_step1_interior<<<blocks,threads>>>(nx, ny, nz, 1.d0, interior, dx, dy, dz, E, F, G, QJ, QJ2)
    
    ! ===== STAGE 2: Boundary Computation (overlaps with halo transfer) =====
    call cpu_time(t_phase_start)
    call calc_EFG_boundary(id_visc, nx, ny, nz, bnd_z_lower, bnd_z_upper, xix, etay, zetaz, &
                           Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
    call cpu_time(t_phase_end)
    call record_phase("boundary_compute", t_phase_end - t_phase_start, metrics)
    
    ! ===== STAGE 3: Synchronize + Finish =====
    call cpu_time(t_phase_start)
    if (mod(myrank,2) == 0) then
      call wait_halo_z(myrank, nx, ny, nz, 3, QJ2, halo_state, ierr)
    endif
    call cpu_time(t_phase_end)
    call record_phase("wait_halo", t_phase_end - t_phase_start, metrics)
    
    ! Boundary step update (now halo data is available)
    call calc_step1_boundary<<<blocks,threads>>>(nx, ny, nz, 1.d0, bnd_z_lower, bnd_z_upper, &
                                                  dx, dy, dz, E, F, G, QJ, QJ2)
    
    ! Apply boundary condition to interior (safe: no neighbor data needed)
    call set_bc_interior(myrank, nx, ny, nz, interior, Jacobian, QJ2)
    
    ! Repeat for stage 2 with modified coefficients...
    ! [Similar pattern for 2nd and 3rd RK stages]
    
    metrics%num_steps = metrics%num_steps + 1
  enddo
  
  ! Report performance
  if (mod(myrank,2) == 0) then
    call report_metrics(myrank, metrics)
  endif
```

---

## Part 4: Implement Interior/Boundary Kernel Split

### Step 4.1: Split KEEP Flux Kernels

**File**: Modify `3D_solver/src/calc_keep_kernel.f90`

Add these new kernels alongside existing ones:

```fortran
  ! NEW: Interior-only x-direction flux
  attributes(global) subroutine calc_keep_x6_interior(id_accuracy, nx, ny, nz, &
                                                       i_start, i_end, j_start, j_end, k_start, k_end, &
                                                       Q, T, E)
    use mod_constant, only : Normal_x
    integer(8), intent(in), value :: id_accuracy
    integer, intent(in), value :: nx, ny, nz
    integer, intent(in), value :: i_start, i_end, j_start, j_end, k_start, k_end
    real(8), intent(in), device :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    
    integer :: i, j, k, it, jt, kt, ii, i_base, idx, offset_yz
    integer, parameter :: sx = threadsE%x + 5
    integer, parameter :: sy = threadsE%y
    integer, parameter :: sz = threadsE%z
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho, u, v, w, p, tmp
    
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1) * sx + (kt-1) * sx * sy
    
    ! Load into shared memory (but bounds-check against interior region)
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= i_start-2 .and. i <= i_end+2 .and. j >= j_start .and. j <= j_end &
         .and. k >= k_start .and. k <= k_end) then
        idx = ii + offset_yz
        rho(idx) = Q(1,i,j,k)
        u(idx) = Q(2,i,j,k)
        v(idx) = Q(3,i,j,k)
        w(idx) = Q(4,i,j,k)
        p(idx) = Q(5,i,j,k)
        tmp(idx) = T(i,j,k)
      endif
    enddo
    call syncthreads()
    
    ! Compute only within interior bounds
    i = (blockIdx%x-1)*blockDim%x + it
    if (i < i_start .or. i > i_end .or. j < j_start .or. j > j_end &
       .or. k < k_start .or. k > k_end) return
    
    idx = it + offset_yz
    
    ! Original flux calculation (unchanged)
    if (3 <= i .and. i <= nx-3) then
      E(:,i,j-1,k-1) = KEEP6(rho(idx-2:idx+3), u(idx-2:idx+3), &
                               v(idx-2:idx+3), w(idx-2:idx+3), &
                              u(idx-2:idx+3), p(idx-2:idx+3), &
                             tmp(idx-2:idx+3), Normal_x)
    ! ... rest of logic
    endif
  end subroutine calc_keep_x6_interior
```

### Step 4.2: Create Wrapper Subroutines

**File**: Add to `src/calc_flux_base.f90`

```fortran
subroutine calc_EFG_interior(id_visc, nx, ny, nz, interior, ...)
  use calc_domain
  type(domain_region), intent(in) :: interior
  
  ! Launch kernels only for interior region
  call calc_keep_x_interior<<<blocks_interior, threads>>>(id_accuracy, nx, ny, nz, &
                                                           interior%i_start, interior%i_end, &
                                                           interior%j_start, interior%j_end, &
                                                           interior%k_start, interior%k_end, &
                                                           Q, T, E)
  call calc_keep_y_interior<<<...>>>(...)
  call calc_keep_z_interior<<<...>>>(...)
  
end subroutine calc_EFG_interior

subroutine calc_EFG_boundary(id_visc, nx, ny, nz, bnd_lower, bnd_upper, ...)
  use calc_domain
  type(domain_region), intent(in) :: bnd_lower, bnd_upper
  
  ! Launch kernels for boundary regions
  call calc_keep_x_boundary<<<blocks_boundary, threads>>>(id_accuracy, nx, ny, nz, &
                                                           bnd_lower, bnd_upper, Q, T, E)
  call calc_keep_y_boundary<<<...>>>(...)
  call calc_keep_z_boundary<<<...>>>(...)
  
end subroutine calc_EFG_boundary
```

---

## Part 5: Validation & Testing

### Step 5.1: Unit Test: Domain Partitioning

**File**: `test_domain_partition.f90`

```fortran
program test_domain_partition
  use calc_domain
  implicit none
  
  type(domain_region) :: interior, bnd_z_lower, bnd_z_upper
  integer :: nx, ny, nz, myrank, nranks
  
  nx = 128
  ny = 128
  nz = 256
  nranks = 4
  
  do myrank = 0, nranks-1
    call partition_domain(nx, ny, nz, myrank, nranks, interior, bnd_z_lower, bnd_z_upper)
    
    ! Verify no overlap
    if (interior%k_start <= bnd_z_lower%k_end) then
      print*, "ERROR: Interior overlaps lower boundary!"
      stop
    endif
    if (bnd_z_upper%k_start <= interior%k_end) then
      print*, "ERROR: Boundary overlaps interior!"
      stop
    endif
    
    ! Verify coverage
    if (bnd_z_lower%k_start /= 1 .or. bnd_z_upper%k_end /= nz) then
      print*, "ERROR: Boundary coverage incomplete!"
      stop
    endif
    
    print*, "Rank", myrank, " OK"
  enddo
  
  print*, "All tests passed!"
end program test_domain_partition
```

### Step 5.2: Integration Test: Taylor-Green Vortex

**File**: `3D_solver/ETGV/test_overlap_tgv.sh`

```bash
#!/bin/bash
# Compare solutions: with vs without overlap

echo "=== Non-overlapped run ==="
cp mod_globals_sync.f90 mod_globals.f90
make clean && make
mpirun -np 4 ./a.out > output_sync.log
mv recal/Q*.dat sync_backup/

echo "=== Overlapped run ==="
cp mod_globals_async.f90 mod_globals.f90
make clean && make
mpirun -np 4 ./a.out > output_async.log
mv recal/Q*.dat async_result/

echo "=== Comparing solutions ==="
# Use Python script to compare field values
python3 compare_solutions.py sync_backup/Q00001.dat async_result/Q00001.dat
```

### Step 5.3: Performance Profiling

```bash
#!/bin/bash
# Profile with and without overlap

echo "=== Baseline (synchronous) ==="
nsys profile -t cuda,mpi -o baseline_ mpirun -np 4 ./a.out
ncu --export baseline.csv baseline_*.nsys-rep

echo "=== Optimized (overlapped) ==="
# [Switch code to overlap version]
nsys profile -t cuda,mpi -o optimized_ mpirun -np 4 ./a.out
ncu --export optimized.csv optimized_*.nsys-rep

echo "=== Analysis ==="
python3 analyze_profiles.py baseline.csv optimized.csv
```

---

## Part 6: Commit Strategy

### Step 6.1: Organize Changes into Logical Commits

```bash
# Commit 1: Foundation modules
git add src/calc_domain.f90 src/calc_halo_exchange.f90 src/calc_overlap_manager.f90
git commit -m "feat: add foundation modules for overlapped compute/comm"

# Commit 2: Time-stepping restructure
git add 3D_solver/src/calc_time_dev.f90
git commit -m "refactor: restructure RK loop in 3 phases for overlap"

# Commit 3: Kernel splitting
git add src/calc_flux_interior.f90 src/calc_flux_boundary.f90
git commit -m "feat: split flux kernels into interior/boundary for selective compute"

# Commit 4: Documentation
git add docs/MULTI_GPU_SCALABILITY_PLAN.md docs/IMPLEMENTATION_GUIDE.md
git commit -m "docs: add implementation guide for multi-GPU scalability"

# Commit 5: Testing
git add test_domain_partition.f90 3D_solver/ETGV/test_overlap_tgv.sh
git commit -m "test: add unit and integration tests for overlap feature"
```

### Step 6.2: Feature Branch Workflow

```bash
# Create feature branch
git checkout -b feature/overlapped-compute-comm
git push origin feature/overlapped-compute-comm

# Develop incrementally, commit regularly
# When ready, create Pull Request
# Request review from team members
# Merge after validation
```

---

## Part 7: Troubleshooting & Common Issues

### Issue 1: Data Race in Boundary Exchange

**Symptom**: Incorrect values at subdomain boundaries  
**Cause**: MPI recv buffer overwritten before GPU read

**Fix**:
```fortran
! WRONG: Recv into same buffer used by GPU
call post_halo_z(...QJ_boundaries...)
call calc_EFG_interior(...)  ! GPU might read QJ boundaries
call wait_halo_z(...)

! CORRECT: Use separate recv buffer
type(halo_exchange_state)
  real(8), device, pointer :: QJ_recv_lower(:,:,:,:)  ! Separate buffer
  real(8), device, pointer :: QJ_recv_upper(:,:,:,:)
end type
```

### Issue 2: Memory Pressure from Extra Buffers

**Symptom**: CUDA out-of-memory errors  
**Cause**: Pinned staging buffers + device buffers exceed GPU memory

**Fix**:
```fortran
! Reduce ghostwidth or use streaming transfer
! Or use CUDA-aware MPI to avoid pinned buffers
#ifdef CUDA_AWARE_MPI
  call MPI_Isend(QJ(1,1,1,nz-gw), ...)  ! Direct GPU→GPU
#else
  ! Pinned buffer approach (more memory but portable)
#endif
```

### Issue 3: Performance Not Improving

**Symptom**: Overlap doesn't help despite code changes  
**Root Causes**:
1. Interior compute time < halo wait time (nothing to overlap)
2. Halo exchange too fast → synchronizes anyway
3. GPU kernel overhead masks benefit

**Debugging**:
```fortran
! Add detailed timers
subroutine time_phases(myrank, interior_t, post_halo_t, boundary_t, wait_halo_t)
  if (myrank == 0) then
    print*, "Interior compute:", interior_t, "ms"
    print*, "Post halo:", post_halo_t, "ms"
    print*, "Boundary compute:", boundary_t, "ms"  ! Should overlap with wait_halo
    print*, "Wait halo:", wait_halo_t, "ms"
    print*, "Overlap potential:", min(boundary_t, wait_halo_t), "ms"
  endif
end subroutine
```

---

## Part 8: Adaptive Configuration

**For Different Problem Sizes**:

| Grid Size | Recommended | Ghostwidth | Blocks_Interior |
|-----------|------------|-----------|-----------------|
| 256³ | Yes (2-4 GPU) | 2 | 32×16 |
| 512³ | Yes (4-8 GPU) | 3 | 64×32 |
| 1024³ | Required (8+ GPU) | 3 | 128×64 |

---

**Next Action**: Start with **minimal_overlap** PoC to validate concept before full implementation.

