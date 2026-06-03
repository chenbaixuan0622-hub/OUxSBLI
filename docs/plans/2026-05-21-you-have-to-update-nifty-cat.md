# Plan: 3D_solver/STZ — Sod Shock Tube in Z (z-decomposition validation)

## Context

Validate the newly implemented `RungeKutta_4th_zdec` and the 4th-order viscous koff dispatch
added in the prior session. The test uses a 1D Sod shock tube propagating in the z-direction
with z-direction MPI domain decomposition on a single GPU (two MPI ranks, both on GPU 0).

---

## Key Findings from Exploration

### Grid layout (cell-centred, matching `set_grid_cyclic4_3D` pattern)
- For 2nd-order accuracy (`id_accuracy` kind=2): overlap=1, `nz_int = nz - 2` interior cells/rank.
- Global spacing: `dz1 = Lz / (nranks * nz_int)`.
- Rank `myrank` owns cells at global indices `[myrank*nz_int+1 .. (myrank+1)*nz_int]`.
- Ghost lo (k=1): centre = `(myrank*nz_int - 0.5)*dz1`; received from rank_lo via MPI.
- Ghost hi (k=nz): centre = `((myrank+1)*nz_int + 0.5)*dz1`; received from rank_hi via MPI.

### `set_bc_cyclic` MUST NOT be called in z-decomp
`set_bc_cyclic2/4/6` overwrites z ghost cells with the local periodic wrap-around, destroying
the MPI-exchanged values. A custom `set_bc` is required that:
- Applies periodic BC only in x and y (for all z including ghost layers).
- Applies zero-gradient (extrapolation) BC in z **only for the outermost ranks**
  (`myrank==0` for the lo face; `myrank==nranks-1` for the hi face).

### `main.f90` must be patched for z-decomp
Currently `set_block3` and `set_init` are guarded by `if (mod(myrank,2) == 0)`.
In z-decomposition ALL ranks are compute ranks, so both calls must happen for every rank.
Fix: change both guards to `if (mod(myrank,2) == 0 .or. kind(id_RungeKutta) == 8)`.
Same change for the final Q-save block. This is backwards-compatible (kind=8 is only true
for the zdec routines).

### `mygpu` assignment
`mygpu = myrank / 2` (already in main.f90). With 2 ranks on 1 GPU: rank 0→GPU 0,
rank 1→GPU 0. No change needed.

### VTK output directories
`print_vtk` for z-decomp writes to `data/<myrank>/`. Create `data/0/` and `data/1/` in
`calc.sh` with `mkdir -p`.

### Generic interface dispatch for `RungeKutta_4th_zdec`
Requires `integer(8)` for `id_RungeKutta` and `integer(4)` for `id_rescale`. Both are set
in `mod_globals.f90`.

### `set_Jacobian_xy3`
For uniform dx, dy, dz: `Jacobian(i,j) = 1/(dx*dy*dz_global)`.
`pre_calc_zdec` divides Q by Jacobian → `QJ = Q_physical * cell_volume`.

---

## Files to Modify

### 1. `3D_solver/src/main.f90`

Three `if (mod(myrank,2) == 0)` guards need to expand to also cover zdec ranks:

```fortran
! BEFORE:
if (mod(myrank,2) == 0) then
  call set_block3(...)
endif

! AFTER:
if (mod(myrank,2) == 0 .or. kind(id_RungeKutta) == 8) then
  call set_block3(...)
endif
```

Apply the same change to:
- The `set_block3` call (lines ~28–31)
- The `set_init` / `id_recal` block (lines ~36–66)
- The final Q-save block (lines ~72–98)

---

## Files to Create: `3D_solver/STZ/`

### 2. `3D_solver/STZ/mod_globals.f90`

```fortran
module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer(2), parameter      :: id_visc     = 2   ! kind=2 → Euler
  real(2),    parameter      :: id_scheme   = 0   ! real(2) → SLAU
  integer, parameter         :: sp          = kind(1.d0)
  real(sp),   parameter      :: threshold   = 0.4_sp
  integer(kind=2), parameter :: id_accuracy = 0   ! kind=2 → 2nd order (overlap=1)
  integer(kind=4), parameter :: id_tvd      = 0   ! kind=4 → minmod limiter
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=4), parameter :: id_rescale  = 0   ! kind=4 → selects RungeKutta_4th_zdec
  integer(kind=2), parameter :: id_gpumpi   = 0
  real(8),    parameter      :: blt         = 0.d0

  ! mesh — nz is LOCAL per rank; global nz = nranks*(nz-2)
  real(8), parameter :: Lx = 0.1d0
  real(8), parameter :: Ly = 0.1d0
  real(8), parameter :: Lz = 1.0d0
  integer, parameter :: nx = 8    ! small: 6 interior cells
  integer, parameter :: ny = 8
  integer, parameter :: nz = 130  ! 128 interior + 1 ghost each end

  logical, parameter :: id_bc_x = .false.  ! periodic x (uniform in x)
  logical, parameter :: id_bc_y = .false.  ! periodic y (uniform in y)
  logical, parameter :: id_bc_z = .true.   ! non-periodic z (extrapolation)

  integer, parameter :: nre1   = 1
  integer, parameter :: nre2   = nx
  integer, parameter :: rerank = 0

  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,4,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,4)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,4)
  type(dim3), parameter :: threads   = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time — RK4 z-decomposition
  ! id_RungeKutta kind=8 → zdec; id_rescale kind=4 → selects RungeKutta_4th_zdec
  integer(kind=2), parameter :: id_recal      = 0
  integer(kind=8), parameter :: id_RungeKutta = 0
  integer, parameter         :: step_offset   = 0
  integer, parameter         :: start_rescale = 0

  ! Physical (dimensionless Sod shock tube)
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0  ! unused for Euler, kept for NS switch

  real(8), parameter :: rho_L = 1.0d0          ! left state density
  real(8), parameter :: p_L   = 1.0d0          ! left state pressure
  real(8), parameter :: rho_R = 0.125d0         ! right state density
  real(8), parameter :: p_R   = 0.1d0           ! right state pressure

  ! CFL-based dt: c_max ~ sqrt(gamma*p_L/rho_L) + speed ~ 1.5
  ! dz_global = Lz/(nranks*(nz-2)) = 1/(2*128) ≈ 3.9e-3
  ! dt = CFL * dz / c_max ≈ 0.5 * 3.9e-3 / 1.5 ≈ 1.3e-3
  real(8), parameter :: dt = 1.0d-3
  integer, parameter :: np = 20   ! print 20 snapshots
  integer, parameter :: nt = 10   ! 10 steps per snapshot → T_end = 0.2
end module mod_globals
```

### 3. `3D_solver/STZ/set.f90`

```fortran
module set
  use mod_globals, only : nx, ny, nz, gamma, rho_L, p_L, rho_R, p_R, Lx, Ly, Lz
  use set_bc_common
  use set_coordinate
  implicit none
contains

  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    use mpi
    use mod_globals, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    real(8) dx1, dy1, dz1, xf(nx+1), yf(ny+1), zf(nz+1)
    integer nranks, ierr, nz_int, iz_offset, i, j, k
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! 2nd-order: overlap=1, interior per rank = nz-2
    nz_int    = nz - 2
    iz_offset = myrank * nz_int
    dx1 = Lx / dble(nx - 2)
    dy1 = Ly / dble(ny - 2)
    dz1 = Lz / dble(nranks * nz_int)
    dx(:) = dx1;  dy(:) = dy1;  dz(:) = dz1
    ! x: cell-face positions then cell centres
    do i = 2, nx;  xf(i) = dx1 * dble(i-2);  enddo
    xf(1) = xf(2) - dx1;  xf(nx+1) = xf(nx) + dx1
    do i = 1, nx;  x(i) = 0.5d0*(xf(i)+xf(i+1));  enddo
    ! y
    do j = 2, ny;  yf(j) = dy1 * dble(j-2);  enddo
    yf(1) = yf(2) - dy1;  yf(ny+1) = yf(ny) + dy1
    do j = 1, ny;  y(j) = 0.5d0*(yf(j)+yf(j+1));  enddo
    ! z: rank-local slab, with ghost face positions
    do k = 2, nz;  zf(k) = dz1 * dble(iz_offset + k - 2);  enddo
    zf(1) = zf(2) - dz1;  zf(nz+1) = zf(nz) + dz1
    do k = 1, nz;  z(k) = 0.5d0*(zf(k)+zf(k+1));  enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,5,ny,nz)
    integer i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (z(k) < Lz * 0.25d0) then
            Q(i,1,j,k) = rho_L
            Q(i,5,j,k) = p_L / (gamma - 1.d0)
          else
            Q(i,1,j,k) = rho_R
            Q(i,5,j,k) = p_R / (gamma - 1.d0)
          endif
          Q(i,2,j,k) = 0.d0   ! rho*u
          Q(i,3,j,k) = 0.d0   ! rho*v
          Q(i,4,j,k) = 0.d0   ! rho*w
        enddo
      enddo
    enddo
  end subroutine set_init


  ! Custom BC for z-decomposition:
  !   - periodic in x and y (for all z including ghost layers)
  !   - zero-gradient in z only for outermost ranks
  !   - does NOT call set_bc_cyclic (which would corrupt z ghost cells)
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ)
    use mpi
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,5,ny,nz)
    integer nranks, ierr, i, j, k, l
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! x periodic (2nd-order, 1 ghost each side, all z)
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do j = 2, ny-1
        do l = 1, 5
          QJ(1,l,j,k)  = QJ(nx-1,l,j,k)
          QJ(nx,l,j,k) = QJ(2,l,j,k)
        enddo
      enddo
    enddo
    ! y periodic
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 2, nx-1
        do l = 1, 5
          QJ(i,l,1,k)  = QJ(i,l,ny-1,k)
          QJ(i,l,ny,k) = QJ(i,l,2,k)
        enddo
      enddo
    enddo
    ! x/y corner ghost cells
    !$cuf kernel do(1)<<<*,*>>>
    do k = 1, nz
      do l = 1, 5
        QJ(1,l,1,k)   = QJ(nx-1,l,ny-1,k)
        QJ(nx,l,1,k)  = QJ(2,l,ny-1,k)
        QJ(1,l,ny,k)  = QJ(nx-1,l,2,k)
        QJ(nx,l,ny,k) = QJ(2,l,2,k)
      enddo
    enddo
    ! z lo boundary: zero-gradient on rank 0
    if (myrank == 0) then
      !$cuf kernel do(2)<<<*,*>>>
      do j = 1, ny
        do i = 1, nx
          do l = 1, 5
            QJ(i,l,j,1) = QJ(i,l,j,2)
          enddo
        enddo
      enddo
    endif
    ! z hi boundary: zero-gradient on last rank
    if (myrank == nranks - 1) then
      !$cuf kernel do(2)<<<*,*>>>
      do j = 1, ny
        do i = 1, nx
          do l = 1, 5
            QJ(i,l,j,nz) = QJ(i,l,j,nz-1)
          enddo
        enddo
      enddo
    endif
  end subroutine set_bc


  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set
```

### 4. `3D_solver/STZ/Makefile`

Copy from `3D_solver/ETGV/Makefile` verbatim (vpath and object list are identical).
Only the `TARGET` line name is cosmetic; the build target `a.out` is what `calc.sh` runs.

### 5. `3D_solver/STZ/calc.sh`

```bash
#!/bin/bash
mkdir -p data/0 data/1 data
mpiexec -n 2 ./a.out
cp mod_globals.f90 ./data
cp set.f90 ./data
```

---

## Grid Consistency Verification

With nranks=2, nz=130 (local), 2nd-order accuracy:
- `nz_int = 128` interior cells per rank, `dz1 = 1.0/(2×128) = 1/256`
- Rank 0 cell centres: `z(k) = (iz_offset + k - 1.5)*dz1` for k=2..129
  → z(2)=0.5/256, z(129)=127.5/256
- Rank 1 cell centres start at z(2)=128.5/256 (contiguous with rank 0 ✓)
- MPI exchange: `flatten_z_hi` of rank 0 packs k=nz-2*1=128 (=127.5/256);
  `reconstruct_z_lo` of rank 1 fills k=1 (ghost) with that value → correct adjacency ✓
- Shock at z=0.25=64/256: cells k≤65 are left state, k≥66 right state in rank 0;
  all rank 1 cells are right state ✓

---

## Verification Steps

1. **Build**: `cd 3D_solver/STZ && make clean && make` — must compile without error.
2. **Run**: `bash calc.sh` — 2 ranks, ~200 timesteps, T_end=0.2.
3. **Check output**: Open `data/0/Q*.vtr` and `data/1/Q*.vtr` in ParaView, combine
   the two z-slabs; the final t=0.2 snapshot should show the classic Sod profile:
   - Rarefaction fan propagating left (z decreasing)
   - Contact discontinuity and shock wave propagating right (z increasing)
4. **Compare with exact Riemann**: `ouxsbli/tests/utils/sod_exact.py` can generate
   the reference profile at t=0.2 for comparison.
5. **Switch to NS**: Change `integer(4)` → `integer(4)` for `id_visc` (keep kind=4 for NS)
   to also exercise `calc_Ev4_koff`, `calc_Fv4_koff`, `calc_Gv4_koff`.

---

## Critical Files

| File | Action |
|------|--------|
| `3D_solver/src/main.f90` | Patch 3 guards: `mod(myrank,2)==0` → also cover `kind(id_RungeKutta)==8` |
| `3D_solver/STZ/mod_globals.f90` | Create (Euler, SLAU, 2nd-order, RK4_zdec) |
| `3D_solver/STZ/set.f90` | Create (custom set_grid + set_init + set_bc for z-decomp) |
| `3D_solver/STZ/Makefile` | Copy from ETGV |
| `3D_solver/STZ/calc.sh` | Create (mkdir data/0 data/1; mpiexec -n 2) |
