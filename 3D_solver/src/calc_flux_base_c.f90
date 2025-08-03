module calc_flux_base_c
  use iso_c_binding
  use calc_flux_base
  implicit none
contains
  subroutine calc_EFG_Euler_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr, QJ_ptr, E_ptr, F_ptr, G_ptr) bind(c, name="calc_EFG_Euler_c")
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), target     :: dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr
    real(8), intent(inout), target     :: QJ_ptr, E_ptr, F_ptr, G_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz, Jacobian
    real(8), dimension(:,:,:,:), device, pointer :: QJ, E, F, G
    integer(kind=2) id_visc
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(Jacobian_ptr), Jacobian, [ny])
    call c_f_pointer(c_loc(QJ_ptr), QJ, [5,nx,ny,nz])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call calc_EFG(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
  end subroutine calc_EFG_Euler_c


  subroutine calc_EFG_Euler_forcing_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr, QJ_ptr, E_ptr, F_ptr, G_ptr, &
                                      fx_ptr, fy_ptr, fz_ptr) bind(c, name="calc_EFG_Euler_forcing_c")
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), target     :: dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr
    real(8), intent(inout), target     :: QJ_ptr, E_ptr, F_ptr, G_ptr, fx_ptr, fy_ptr, fz_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz, Jacobian
    real(8), dimension(:,:,:,:), device, pointer :: QJ, E, F, G
    real(8), dimension(:,:,:), device, pointer   :: fx, fy, fz
    integer(kind=2) id_visc
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(Jacobian_ptr), Jacobian, [ny])
    call c_f_pointer(c_loc(QJ_ptr), QJ, [5,nx,ny,nz])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call c_f_pointer(c_loc(fx_ptr), fx, [nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(fy_ptr), fy, [nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(fz_ptr), fz, [nx-2,ny-2,nz-2])
    call calc_EFG(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G, fx, fy, fz)
  end subroutine calc_EFG_Euler_forcing_c

 
  subroutine calc_EFG_NS_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr, QJ_ptr, E_ptr, F_ptr, G_ptr) bind(c, name="calc_EFG_NS_c")
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), target     :: dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr
    real(8), intent(inout), target     :: QJ_ptr, E_ptr, F_ptr, G_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz, Jacobian
    real(8), dimension(:,:,:,:), device, pointer :: QJ, E, F, G
    integer(kind=4) id_visc
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(Jacobian_ptr), Jacobian, [ny])
    call c_f_pointer(c_loc(QJ_ptr), QJ, [5,nx,ny,nz])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call calc_EFG(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
  end subroutine calc_EFG_NS_c


  subroutine calc_EFG_NS_forcing_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr, QJ_ptr, E_ptr, F_ptr, G_ptr, &
                                   fx_ptr, fy_ptr, fz_ptr) bind(c, name="calc_EFG_NS_forcing_c")
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), target     :: dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr
    real(8), intent(inout), target     :: QJ_ptr, E_ptr, F_ptr, G_ptr, fx_ptr, fy_ptr, fz_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz, Jacobian
    real(8), dimension(:,:,:,:), device, pointer :: QJ, E, F, G
    real(8), dimension(:,:,:), device, pointer   :: fx, fy, fz
    integer(kind=4) id_visc
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(Jacobian_ptr), Jacobian, [ny])
    call c_f_pointer(c_loc(QJ_ptr), QJ, [5,nx,ny,nz])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call c_f_pointer(c_loc(fx_ptr), fx, [nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(fy_ptr), fy, [nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(fz_ptr), fz, [nx-2,ny-2,nz-2])
    call calc_EFG(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G, fx, fy, fz)
  end subroutine calc_EFG_NS_forcing_c


  subroutine calc_EFG_LES_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr, QJ_ptr, E_ptr, F_ptr, G_ptr) bind(c, name="calc_EFG_LES_c")
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), target     :: dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr
    real(8), intent(inout), target     :: QJ_ptr, E_ptr, F_ptr, G_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz, Jacobian
    real(8), dimension(:,:,:,:), device, pointer :: QJ, E, F, G
    integer(kind=8) id_visc
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(Jacobian_ptr), Jacobian, [ny])
    call c_f_pointer(c_loc(QJ_ptr), QJ, [5,nx,ny,nz])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call calc_EFG(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
  end subroutine calc_EFG_LES_c


  subroutine calc_EFG_LES_forcing_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr, QJ_ptr, E_ptr, F_ptr, G_ptr, &
                                     fx_ptr, fy_ptr, fz_ptr) bind(c, name="calc_EFG_LES_forcing_c")
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), target     :: dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr
    real(8), intent(inout), target     :: QJ_ptr, E_ptr, F_ptr, G_ptr, fx_ptr, fy_ptr, fz_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz, Jacobian
    real(8), dimension(:,:,:,:), device, pointer :: QJ, E, F, G
    real(8), dimension(:,:,:), device, pointer   :: fx, fy, fz
    integer(kind=8) id_visc
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(Jacobian_ptr), Jacobian, [ny])
    call c_f_pointer(c_loc(QJ_ptr), QJ, [5,nx,ny,nz])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call c_f_pointer(c_loc(fx_ptr), fx, [nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(fy_ptr), fy, [nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(fz_ptr), fz, [nx-2,ny-2,nz-2])
    call calc_EFG(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G, fx, fy, fz)
  end subroutine calc_EFG_LES_forcing_c
end module calc_flux_base_c

