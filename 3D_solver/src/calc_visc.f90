module calc_visc
  use mod_globals, only : id_visc, id_turbulence, id_av, gamma, R, Pr, Prt
  use mod_constant, only : Cp, gamma_1
  use calc_visc_base
  use calc_me4_base
  use calc_visc_visbal
  implicit none
contains
  attributes(global) subroutine calc_Ev(nx, ny, nz, dx, dy, dz, Q, E)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k
    real(8) txx, txy, txz, utxx, vtxy, wtxz, kTx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      call calc_me4_base_x(nx, ny, nz, i, j, k, dx, dy, dz, Q, txx, txy, txz, utxx, vtxy, wtxz, kTx)
    else
      call calc_derivative_x2(nx, ny, nz, i, j, k, dx, dy, dz, Q, txx, txy, txz, utxx, vtxy, wtxz, kTx)
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
  end subroutine calc_Ev
  
  attributes(global) subroutine calc_Ev_LES(nx, ny, nz, dx, dy, dz, Q, mut, qc2, E)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k
    real(8) txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      call calc_me4_base_les_x(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs)
    else
      call calc_derivative_les_x2(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs)
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx + Hsgs)
  end subroutine calc_Ev_LES
  
  attributes(global) subroutine calc_Fv(nx, ny, nz, dy, dx, dz, Q, F)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    ! 2nd-order 2
    real(8), dimension(3,2,3), device :: T323
    real(8), dimension(3,2), device   :: tmp2
    real(8), dimension(2,3), device   :: tmp3
    real(8), dimension(2), device     :: Ty, u2, v2, w2, mz, mx
    real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      call calc_me4_base_y(nx, ny, nz, i, j, k, dx, dy, dz, Q, tyx, tyy, tyz, utyx, vtyy, wtyz, kTy)
    else
      call calc_derivative_y2(nx, ny, nz, i, j, k, dx, dy, dz, Q, tyx, tyy, tyz, utyx, vtyy, wtyz, kTy)
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv
  
  attributes(global) subroutine calc_Fv_LES(nx, ny, nz, dy, dx, dz, Q, mut, qc2, F)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      call calc_me4_base_les_y(nx, ny, nz, i, j, k, dy, dx, dz, Q, mut, qc2, tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs)
    else
      call calc_derivative_les_y2(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs)
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy + Hsgs)
  end subroutine calc_Fv_LES
  
  attributes(global) subroutine calc_Gv(nx, ny, nz, dx, dy, dz, Q, G)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      call calc_me4_base_z(nx, ny, nz, i, j, k, dx, dy, dz, Q, tzx, tzy, tzz, utzx, vtzy, wtzz, kTz)
    else
      call calc_derivative_z2(nx, ny, nz, i, j, k, dx, dy, dz, Q, tzx, tzy, tzz, utzx, vtzy, wtzz, kTz)
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv

  attributes(global) subroutine calc_Gv_LES(nx, ny, nz, dx, dy, dz, Q, mut, qc2, G)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      call calc_me4_base_les_z(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs)
    else
      call calc_derivative_les_z2(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs)
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz + Hsgs)
  end subroutine calc_Gv_LES
end module calc_visc

