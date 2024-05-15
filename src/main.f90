module mod_allocate
  implicit none
  interface alloc
    module procedure alloc2D, alloc3D
  end interface
contains
  subroutine alloc2D(Q,z)
    use mod_globals, only : nx, ny
    real(8), intent(out), allocatable :: Q(:,:,:)
    real(8), intent(out), allocatable :: z(:)
    allocate(Q(nx,ny,4),z(1))
  end subroutine alloc2D

  subroutine alloc3D(Q,z)
    use mod_globals, only : nx, ny, nz
    real(8), intent(out), allocatable :: Q(:,:,:,:)
    real(8), intent(out), allocatable :: z(:)
    allocate(Q(nx,ny,nz,5),z(nz))
  end subroutine alloc3D
end module mod_allocate

program main
  use, intrinsic :: iso_fortran_env
  use mod_allocate
  use mod_globals, only : id_RungeKutta, nx, ny, nz, Q
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  real(8) t_start, t_end
  real(8), allocatable :: x(:), xix(:), dx(:), y(:), etay(:), dy(:), z(:), Jacobian(:,:)
  real(8), allocatable :: Vin(:,:)
  allocate(x(nx),xix(nx),dx(nx),y(ny),etay(ny),dy(ny),Jacobian(nx,ny),Vin(ny,2))

  call alloc(Q,z)

  call set_grid(x,y,z,dx,dy)
  call set_xix(dx,xix)
  call set_etay(dy,etay)
  call set_Jacobian(dx,dy,Jacobian)
  call set_init(x,y,z,Q,Vin)

  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta,x,dx,xix,y,dy,etay,z,Jacobian,Q,Vin)
  call cpu_time(t_end)
  print *, "elapsed time:", t_end - t_start

  deallocate(Q,x,xix,dx,y,etay,dy,z,Jacobian,Vin)
end program main

