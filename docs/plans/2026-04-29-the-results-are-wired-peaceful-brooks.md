# Plan: Fix WING uniform-flow bug + Add CORN (compression corner) test case

## Context

The WING curvilinear O-grid solver produces "almost uniform" flow because the shared
`src/main.f90` uses `use calc_time_dev` (Cartesian solver) and calls `RungeKutta`, not
`RungeKutta_curv` from `calc_time_dev_curv`. No WING-specific `main.f90` exists in
`3D_solver/WING/`, so the WING Makefile (which lists `main.o` depending on
`calc_time_dev_curv.mod`) cannot link `RungeKutta_curv` — causing a link failure or, if a
stale object exists, running the Cartesian solver on a domain with trivial periodic Cartesian
metrics, which produces zero flux divergence everywhere → uniform flow.

A second issue: `set_metrics_curv` (in `src/set_coordinate.f90`) stores the raw 2-D Jacobian
`J_2D = ∂x/∂ξ · ∂y/∂η − ∂x/∂η · ∂y/∂ξ` in `Jac_cpu`. The curvilinear solver convention
requires `Jac_code = 1/(J_2D · dz)` so that `QJ = Q_physical / Jac_code = Q_physical · J_2D · dz`
(cell-volume-scaled conserved variable). This convention is what makes
`dt_Szeta = dt/(Jac_code·dz) = dt·J_2D` correct, and what makes the far-field BC formula
`Q(i,1,ny,k) = rho_inf / Jacobian(i,ny)` yield `rho_inf · J_2D · dz` (correct QJ value).
The conversion `Jac_code = 1/(J_2D · dz)` must be done in `main.f90` before calling
`RungeKutta_curv`, since `set_metrics_curv` doesn't take `dz` as input.

---

## Part 1 — Fix WING: create `3D_solver/WING/main.f90`

Three changes vs `src/main.f90`:
1. `use calc_time_dev_curv` instead of `use calc_time_dev`
2. Call `set_metrics(nx, ny)` after `set_grid` to fill `Jac_cpu` (= J_2D) and metric arrays
3. Convert `Jac_cpu → Jacobian = 1/(J_2D · dz(1))` before passing to `RungeKutta_curv`
4. Call `RungeKutta_curv(...)` with `x_phys_g, y_phys_g` and all metric arrays from `set` module

**File**: `3D_solver/WING/main.f90`

```fortran
program main
  use, intrinsic :: iso_fortran_env
  use mpi
  use mod_globals, only : id_RungeKutta, id_rescale, id_recal, nx, ny, nz, Lx, Ly, Lz, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use set
  use set_coordinate, only : set_block3
  use calc_time_dev_curv
  implicit none
  integer i, j, l, m, s, mygpu, ios, errorcode
  real(8) t_start, t_end
  real(8), allocatable :: x(:), dx(:), y(:), dy(:), z(:), dz(:), Jacobian(:,:), Q(:,:,:,:)
  character(len=8) header
  character(len=40) filename
  logical is_sequential
  integer nranks, myrank, ierr, ireq, istat(MPI_STATUS_SIZE)

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  mygpu = myrank / 2

  if (mod(myrank,2) == 0) then
    call set_block3(nx, ny, nz, threads, threadsE, threadsEv, threadsF, threadsFv, &
                    threadsG, threadsGv, &
                    blocks, blocksE, blocksEv, blocksF, blocksFv, blocksG, blocksGv)
  endif

  allocate(Q(nx,5,ny,nz), x(nx), dx(nx-1), y(ny), dy(ny-1), z(nz), dz(nz-1), Jacobian(nx,ny))

  call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
  call set_metrics(nx, ny)

  ! Jac_cpu holds J_2D from set_metrics_curv.
  ! Convert to Jac_code = 1/(J_2D*dz) for pre_calc_curv: QJ = Q*J_2D*dz, dt_Szeta = dt*J_2D.
  do j = 1, ny
    do i = 1, nx
      Jacobian(i,j) = 1.d0 / (Jac_cpu(i,j) * dz(1))
    enddo
  enddo

  if (mod(myrank,2) == 0) then
    if (kind(id_recal) == 4) then
      write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
      open(10, file=filename, action="read", form="unformatted", access="sequential", &
           status="old", iostat=ios)
      if (ios /= 0) then
        print *, "Error opening file."
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
        stop
      endif
      read(10, iostat=ios) header
      close(10)
      is_sequential = (ios == 0 .and. header == 'SEQFMT01')
      if (is_sequential) then
        open(10, file=filename, action="read", form="unformatted", access="sequential", status="old")
        read(10) header
        read(10) Q
        print *, "myrank is ", myrank, "simulation has been restarted. access is sequential"
      else
        open(10, file=filename, action="read", form="unformatted", access="stream", status="old")
        rewind(10)
        read(10) Q
        print *, "myrank is ", myrank, "simulation has been restarted. access is stream"
      endif
      close(10)
    elseif (kind(id_recal) == 2) then
      write(*,*) "set initial condition"
      call set_init(myrank, nx, ny, nz, x, y, z, Q)
    else
      write(*,*) "wrong parameter was found"
    endif
  endif

  call cpu_time(t_start)
  call RungeKutta_curv(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, &
                       x_phys_g, y_phys_g, z, dz(1), &
                       Jacobian, &
                       n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
                       xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, Q)
  call cpu_time(t_end)

  if (mod(myrank,2) == 0) then
    if (t_end - t_start <= 60.d0) then
      s = int(t_end - t_start)
      print *, "calculation time:", s, " [sec]"
    else
      m = int(t_end - t_start) / 60
      s = int(t_end - t_start) - 60 * m
      print *, "calculation time:", m, " [min] ", s, " [sec]"
    endif
    do l = 1, nz
      do j = 1, ny
        do i = 1, nx
          do m = 1, 5
            Q(i,m,j,l) = Jacobian(i,j) * Q(i,m,j,l)
    enddo;enddo;enddo;enddo
    call cpu_time(t_start)
    write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
    open(10, file=filename, status="replace", action="write", form="unformatted", access="stream")
    write(10) Q
    close(10)
    call cpu_time(t_end)
    s = t_end - t_start
    print *, "output time:", s, " [sec]"
  endif

  deallocate(Q, x, dx, y, dy, z, dz, Jacobian)
  call MPI_FINALIZE(ierr)
end program main
```

---

## Part 2 — Add compression corner reflection test case: `3D_solver/CORN/`

### Purpose

Validate the curvilinear grid solver with a canonical supersonic flow problem that has an
analytical solution (oblique shock polars). The body-fitted mesh is non-orthogonal at the
corner, exercising all metric computation paths. "Reflection" = oblique shock from the
lower-wall ramp hits the flat upper wall and reflects back.

---

### 2a — Add `set_grid_c_corner` to `src/set_coordinate.f90`

Insert before `set_grid_c_wing` in the `contains` section. Follows the same coding style.

```fortran
  ! Generate body-fitted mesh for a compression corner.
  ! Lower boundary: flat plate (y=0) for xi <= x_corner, then ramp at angle theta.
  ! Upper boundary: flat wall at y = Ly.
  ! xi (i) runs from inlet to outlet (NOT periodic).
  ! eta (j) blends linearly from lower to upper wall.
  ! zeta (k) is uniform spanwise.
  subroutine set_grid_c_corner(nx, ny, nz, Lx, Ly, Lz, theta, x_corner, &
                                x_phys, y_phys, zc, dz)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz, theta, x_corner
    real(8), intent(out) :: x_phys(nx,ny), y_phys(nx,ny)
    real(8), intent(out) :: zc(nz), dz(nz-1)
    real(8) :: dxi, tj, xi(nx), yi(nx), xo(nx), yo(nx)
    integer :: i, j, k
    dxi = Lx / dble(nx-2)
    do i = 2, nx-1
      xi(i) = dxi * dble(i-2)
      if (xi(i) <= x_corner) then
        yi(i) = 0.d0
      else
        yi(i) = (xi(i) - x_corner) * tan(theta)
      endif
      xo(i) = xi(i)
      yo(i) = Ly
    enddo
    xi(1)  = xi(2)    - dxi;  yi(1)  = 0.d0
    xi(nx) = xi(nx-1) + dxi;  yi(nx) = yi(nx-1) + dxi * tan(theta)
    xo(1)  = xi(1);            yo(1)  = Ly
    xo(nx) = xi(nx);           yo(nx) = Ly
    do j = 1, ny
      tj = dble(j-1) / dble(ny-1)
      do i = 1, nx
        x_phys(i,j) = xi(i) + tj * (xo(i) - xi(i))
        y_phys(i,j) = yi(i) + tj * (yo(i) - yi(i))
      enddo
    enddo
    do k = 1, nz
      zc(k) = Lz * dble(k-1) / dble(nz-1)
    enddo
    do k = 1, nz-1
      dz(k) = zc(k+1) - zc(k)
    enddo
  end subroutine set_grid_c_corner
```

---

### 2b — `3D_solver/CORN/mod_globals.f90`

```fortran
module mod_globals
  use cudafor
  implicit none

  integer, parameter :: nx = 192, ny = 80, nz = 4

  real(8), parameter :: Lx = 2.d0        ! streamwise domain length
  real(8), parameter :: Ly = 1.d0        ! wall-normal domain height
  real(8), parameter :: Lz = 0.1d0       ! spanwise extent

  ! Compression corner geometry
  real(8), parameter :: theta    = 1.74532925199433d-1  ! 10 degrees [rad]
  real(8), parameter :: x_corner = 0.5d0                ! corner location

  ! Free-stream flow conditions (non-dimensional)
  real(8), parameter :: Ma_inf  = 2.5d0
  real(8), parameter :: gamma   = 1.4d0
  real(8), parameter :: R       = 1.d0
  real(8), parameter :: Pr      = 0.72d0

  integer, parameter :: sp      = 4

  real(8), parameter :: dt      = 5.d-5
  integer, parameter :: nstep   = 10000
  integer, parameter :: nprint  = 1000

  integer, parameter :: nt = nprint
  integer, parameter :: np = nstep / nprint

  integer, parameter :: step_offset = 0
  integer, parameter :: nt_output   = 1
  real(8), parameter :: dt_ref      = 1.d0
  integer, parameter :: rerank      = -1

  integer(2), parameter :: id_visc       = 0  ! kind=2 → Euler
  real(8),    parameter :: id_scheme     = 0  ! kind=8 → Hybrid (shock capturing)
  integer(2), parameter :: id_accuracy   = 0  ! kind=2 → 2nd order
  integer(2), parameter :: id_tvd        = 0  ! kind=2 → no TVD
  integer(2), parameter :: id_slau       = 0  ! kind=2 → SLAU1
  integer(2), parameter :: id_rescale    = 0  ! kind=2 → off
  integer(2), parameter :: id_gpumpi     = 0  ! kind=2 → CPU MPI
  integer(2), parameter :: id_RungeKutta = 0  ! kind=2 → TVD-RK3
  integer(2), parameter :: id_recal      = 0  ! kind=2 → initialize

  real(4), parameter :: threshold = 0.1_sp

  logical, parameter :: id_bc_x = .true.   ! inlet/outlet (not periodic)
  logical, parameter :: id_bc_y = .true.   ! slip walls at j=1 and j=ny
  logical, parameter :: id_bc_z = .false.  ! z-periodic

  integer, parameter :: dimension = 3

  type(dim3), parameter :: threads   = dim3(32, 8, 1)
  type(dim3), parameter :: threadsE  = dim3(32, 4, 2)
  type(dim3), parameter :: threadsFv = dim3(8,  4, 4)
  type(dim3), parameter :: threadsF  = dim3(8,  16, 2)
  type(dim3), parameter :: threadsEv = dim3(16, 4, 2)
  type(dim3), parameter :: threadsG  = dim3(8,  8, 4)
  type(dim3), parameter :: threadsGv = dim3(8,  8, 4)

  type(dim3) :: blocks, blocksE, blocksEv, blocksF, blocksFv, blocksG, blocksGv

end module mod_globals
```

---

### 2c — `3D_solver/CORN/set.f90`

```fortran
module set
  use cudafor
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, theta, x_corner, Ma_inf, gamma, R, Pr
  use set_coordinate, only : set_grid_c_corner, set_metrics_curv
  use set_bc_common
  implicit none

  real(8), allocatable, save :: x_phys_g(:,:), y_phys_g(:,:)
  real(8), allocatable, save :: n_xi_x_cpu(:,:), n_xi_y_cpu(:,:)
  real(8), allocatable, save :: n_eta_x_cpu(:,:), n_eta_y_cpu(:,:)
  real(8), allocatable, save :: xi_x_cpu(:,:), xi_y_cpu(:,:)
  real(8), allocatable, save :: eta_x_cpu(:,:), eta_y_cpu(:,:)
  real(8), allocatable, save :: Jac_cpu(:,:)

contains

  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in) :: Lx, Ly, Lz
    real(8), intent(out), allocatable :: xc(:), yc(:), zc(:)
    real(8), intent(out), allocatable :: dx(:), dy(:), dz(:)
    integer i, j
    allocate(xc(nx), yc(ny), zc(nz))
    allocate(dx(nx-1), dy(ny-1), dz(nz-1))
    allocate(x_phys_g(nx,ny), y_phys_g(nx,ny))

    call set_grid_c_corner(nx, ny, nz, Lx, Ly, Lz, theta, x_corner, &
                           x_phys_g, y_phys_g, zc, dz)

    xc = (/ (dble(i), i=1,nx) /)
    yc = (/ (dble(j), j=1,ny) /)
    dx = 1.d0
    dy = 1.d0

    if (myrank == 0) then
      print *, "CORN compression corner: nx=", nx, " ny=", ny, " nz=", nz
      print *, "  M_inf=", Ma_inf, "  theta=", theta, " rad  x_corner=", x_corner
      print *, "  Domain: Lx=", Lx, " Ly=", Ly, " Lz=", Lz
    endif
  end subroutine set_grid


  subroutine set_metrics(nx, ny)
    integer, intent(in) :: nx, ny

    allocate(n_xi_x_cpu(nx-1,ny-2), n_xi_y_cpu(nx-1,ny-2))
    allocate(n_eta_x_cpu(nx-2,ny-1), n_eta_y_cpu(nx-2,ny-1))
    allocate(xi_x_cpu(nx,ny), xi_y_cpu(nx,ny))
    allocate(eta_x_cpu(nx,ny), eta_y_cpu(nx,ny))
    allocate(Jac_cpu(nx,ny))

    call set_metrics_curv(nx, ny, x_phys_g, y_phys_g, &
                          n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
                          xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, Jac_cpu)
  end subroutine set_metrics


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,5,ny,nz)
    real(8) :: rho_inf, u_inf, v_inf, p_inf, E_inf

    rho_inf = 1.d0
    u_inf   = Ma_inf
    v_inf   = 0.d0
    p_inf   = 1.d0 / gamma
    E_inf   = p_inf / (gamma - 1.d0) + 0.5d0 * Ma_inf**2

    Q = 0.d0
    Q(:,1,:,:) = rho_inf
    Q(:,2,:,:) = rho_inf * u_inf
    Q(:,3,:,:) = rho_inf * v_inf
    Q(:,4,:,:) = 0.d0
    Q(:,5,:,:) = E_inf

    if (myrank == 0) then
      print *, "Flow: M_inf =", Ma_inf
      print *, "  rho_inf =", rho_inf, " p_inf =", p_inf, " E_inf =", E_inf
    endif
  end subroutine set_init


  !> BCs for compression corner:
  !>   i=1   (xi inlet ghost)  : supersonic Dirichlet free-stream
  !>   i=nx  (xi outlet ghost) : zero-gradient extrapolation (supersonic exit)
  !>   j=1   (lower wall ghost): Euler slip wall (ramp)
  !>   j=ny  (upper wall ghost): Euler slip wall (flat reflector)
  !>   k=1, k=nz               : z-periodic
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, eta_x, eta_y, Q)
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in),    device :: Jacobian(nx,ny)
    real(8), intent(in),    device :: eta_x(nx,ny), eta_y(nx,ny)
    real(8), intent(inout), device :: Q(nx,5,ny,nz)
    integer :: i, j, k
    real(8) :: rho_inf, u_inf, p_inf, E_inf, Jratio
    real(8) :: nxw, nyw, nmag, u_int, v_int, u_n

    rho_inf = 1.d0
    u_inf   = Ma_inf
    p_inf   = 1.d0 / gamma
    E_inf   = p_inf / (gamma - 1.d0) + 0.5d0 * Ma_inf**2

    ! (a) xi inlet ghost (i=1): supersonic free-stream Dirichlet
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do j = 1, ny
        Q(1,1,j,k) = rho_inf          / Jacobian(1,j)
        Q(1,2,j,k) = rho_inf * u_inf  / Jacobian(1,j)
        Q(1,3,j,k) = 0.d0
        Q(1,4,j,k) = 0.d0
        Q(1,5,j,k) = E_inf            / Jacobian(1,j)
      enddo
    enddo

    ! (b) xi outlet ghost (i=nx): zero-gradient (physical) extrapolation
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do j = 1, ny
        Jratio = Jacobian(nx-1,j) / Jacobian(nx,j)
        Q(nx,1,j,k) = Q(nx-1,1,j,k) * Jratio
        Q(nx,2,j,k) = Q(nx-1,2,j,k) * Jratio
        Q(nx,3,j,k) = Q(nx-1,3,j,k) * Jratio
        Q(nx,4,j,k) = Q(nx-1,4,j,k) * Jratio
        Q(nx,5,j,k) = Q(nx-1,5,j,k) * Jratio
      enddo
    enddo

    ! (c) eta lower wall (j=1): Euler slip wall, reflect via eta direction at j=1
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do i = 1, nx
        nxw   = eta_x(i,1);  nyw = eta_y(i,1)
        nmag  = sqrt(nxw*nxw + nyw*nyw)
        nxw   = nxw / nmag;  nyw = nyw / nmag
        u_int = Q(i,2,2,k) / Q(i,1,2,k)
        v_int = Q(i,3,2,k) / Q(i,1,2,k)
        u_n   = u_int*nxw + v_int*nyw
        Jratio = Jacobian(i,2) / Jacobian(i,1)
        Q(i,1,1,k) = Q(i,1,2,k) * Jratio
        Q(i,2,1,k) = (Q(i,2,2,k) - 2.d0*u_n*nxw*Q(i,1,2,k)) * Jratio
        Q(i,3,1,k) = (Q(i,3,2,k) - 2.d0*u_n*nyw*Q(i,1,2,k)) * Jratio
        Q(i,4,1,k) = Q(i,4,2,k) * Jratio
        Q(i,5,1,k) = Q(i,5,2,k) * Jratio
      enddo
    enddo

    ! (d) eta upper wall (j=ny): Euler slip wall, reflect via eta direction at j=ny
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do i = 1, nx
        nxw   = eta_x(i,ny);  nyw = eta_y(i,ny)
        nmag  = sqrt(nxw*nxw + nyw*nyw)
        nxw   = nxw / nmag;   nyw = nyw / nmag
        u_int = Q(i,2,ny-1,k) / Q(i,1,ny-1,k)
        v_int = Q(i,3,ny-1,k) / Q(i,1,ny-1,k)
        u_n   = u_int*nxw + v_int*nyw
        Jratio = Jacobian(i,ny-1) / Jacobian(i,ny)
        Q(i,1,ny,k) = Q(i,1,ny-1,k) * Jratio
        Q(i,2,ny,k) = (Q(i,2,ny-1,k) - 2.d0*u_n*nxw*Q(i,1,ny-1,k)) * Jratio
        Q(i,3,ny,k) = (Q(i,3,ny-1,k) - 2.d0*u_n*nyw*Q(i,1,ny-1,k)) * Jratio
        Q(i,4,ny,k) = Q(i,4,ny-1,k) * Jratio
        Q(i,5,ny,k) = Q(i,5,ny-1,k) * Jratio
      enddo
    enddo

    ! (e) z-periodic
    Q(:,:,:,1)  = Q(:,:,:,nz-1)
    Q(:,:,:,nz) = Q(:,:,:,2)
  end subroutine set_bc

end module set
```

---

### 2d — `3D_solver/CORN/main.f90`

Identical to `3D_solver/WING/main.f90` created in Part 1. The CORN-specific geometry
parameters (`theta`, `x_corner`) are encapsulated in CORN's `set_grid` → `set_grid_c_corner`.
No change needed in main.

---

### 2e — `3D_solver/CORN/Makefile`

Copy `3D_solver/WING/Makefile` verbatim (same OBJ list, same vpath `../src:../../src`).

---

### 2f — `3D_solver/CORN/calc.sh`

```bash
#!/bin/bash
mpirun -n 2 ./a.out
```

---

## Critical files

| File | Action |
|------|--------|
| `3D_solver/WING/main.f90` | **CREATE** — primary bug fix |
| `src/set_coordinate.f90` | **MODIFY** — insert `set_grid_c_corner` before `set_grid_c_wing` in `contains` |
| `3D_solver/CORN/mod_globals.f90` | **CREATE** |
| `3D_solver/CORN/set.f90` | **CREATE** |
| `3D_solver/CORN/main.f90` | **CREATE** (copy of WING/main.f90) |
| `3D_solver/CORN/Makefile` | **CREATE** (copy of WING/Makefile) |
| `3D_solver/CORN/calc.sh` | **CREATE** |

Files **not modified**: all `3D_solver/src/calc_*_curv.f90` kernels, WING's `set.f90` and
`mod_globals.f90`, `src/print_curv.f90`, `3D_solver/src/preprocess.f90`,
`3D_solver/src/calc_time_dev_curv.f90`.

---

## Verification

1. `cd 3D_solver/WING && make clean && make 2>&1 | head -50`
   — `main.f90` must compile against `calc_time_dev_curv.mod`, not `calc_time_dev.mod`
2. `cd 3D_solver/CORN && make clean && make 2>&1 | head -50`
3. After `set_metrics(nx, ny)` in either main, `minval(Jac_cpu) > 0` (grid orientation correct)
4. **WING run**: `bash calc.sh` → open `data/rank1/Q00001.vts` in ParaView:
   - Grid is O-grid around NACA 0012
   - Pressure shows stagnation point at LE and circulation around body (not uniform)
5. **CORN run**: `bash calc.sh` → check VTS output:
   - Oblique shock from lower-wall corner at x=0.5
   - Reflected shock off flat upper wall
   - Analytical check: for M=2.5, θ=10°, oblique shock angle β ≈ 31° (θ-β-M relation)
