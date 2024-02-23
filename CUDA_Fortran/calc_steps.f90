module calc_steps
  use calc_flux
  implicit none
contains
  attributes(global) subroutine calc_step1(nx,ny,nz,ni,nj,nk,dX,dY,dZ,dt,gamma,Q,Q2,E,F,G,rho,u,v,w,p)
    integer, value :: nx, ny, nz, ni, nj, nk
    real(8), value :: dX, dY, dZ, dt
    real(8), dimension(nx,ny,nz,5,4), device :: Q
    real(8), dimension(nx,ny,nz,5,4), device :: Q2
    real(8), device :: E(nx-3,ny-4,nz-4,5,4)
    real(8), device :: F(nx-4,ny-3,nz-4,5,4)
    real(8), device :: G(nx-4,ny-4,nz-3,5,4)
    integer :: i, j, k
    k = 3 + nk * mod(int(threadIdx%x/64),8) !k = 3, nz-nk-1, nk :3 <= k <= nz-2
    j = 3 + nj * mod(int(threadIdx%x/8),8)  !j = 3, ny-nj-1, nj :3 <= j <= ny-2
    i = 3 + ni * mod(threadIdx%x,8)         !i = 3, nx-ni-1, ni :3 <= i <= nx-2
    !Q_d = Q(i-2:i+ni+1,j-2:j+nj+1,k-2:k+nk+1,:,blockIdx%x)
    Q2(i:i+ni-1,j:j+nj-1,k:k+nk-1,:,blockIdx%x) = Q(i:i+ni-1,j:j+nj-1,k:k+nk-1,:,blockIdx%x)-dt*(&
    dX*(-E(i-2:i+ni-3,j-2:j+nj-3,k-2:k+nk-3,:,blockIdx%x)+E(i-1:i+ni-2,j-2:j+nj-3,k-2:k+nk-3,:,blockIdx%x))+&
    dY*(-F(i-2:i+ni-3,j-2:j+nj-3,k-2:k+nk-3,:,blockIdx%x)+F(i-2:i+ni-3,j-1:j+nj-2,k-2:k+nk-3,:,blockIdx%x))+&
    dZ*(-G(i-2:i+ni-3,j-2:j+nj-3,k-2:k+nk-3,:,blockIdx%x)+G(i-2:i+ni-3,j-2:j+nj-3,k-1:k+nk-2,:,blockIdx%x)))
  end subroutine calc_step1
  
  attributes(global) subroutine calc_step2(nx,ny,nz,ni,nj,nk,dX,dY,dZ,dt,gamma,mu,kappa,Cp,Q,Q2,Q3,Q_d)
    integer, intent(in), value :: nx, ny, nz, ni, nj, nk
    real(8), intent(in), value :: dX, dY, dZ, dt, gamma, mu, kappa, Cp
    real(8), intent(in), dimension(nx,ny,nz,5,4) :: Q, Q2
    real(8), intent(out), dimension(nx,ny,nz,5,4) :: Q3
    integer :: i, j, k
    real(8), intent(inout) :: Q_d(ni+4,nj+4,nk+4,5)
    real(8), dimension(ni+5,nj+4,nk+4,5) :: E
    real(8), dimension(ni+4,nj+5,nk+4,5) :: F
    real(8), dimension(ni+4,nj+4,nk+5,5) :: G
    !real(8), intent(inout) :: R(ni,nj,nk,5)
    k = 3 + nk * mod(int(threadIdx%x/64),8) !k = 3, nz-nk-1, nk :3 <= k <= nz-2
    j = 3 + nj * mod(int(threadIdx%x/8),8)  !j = 3, ny-nj-1, nj :3 <= j <= ny-2
    i = 3 + ni * mod(threadIdx%x,8)         !i = 3, nx-ni-1, ni :3 <= i <= nx-2
    Q_d = Q2(i-2:i+ni+1,j-2:j+nj+1,k-2:k+nk+1,:,blockIdx%x)
    !call calc_EFG(ni+4,nj+4,nk+4,dX,dY,dZ,gamma,mu,kappa,Cp,Q_d,R)
    !R = dX * (-E(1:nx-4,:,:,:) + E(2:nx-3,:,:,:))&
    !+ dY * (-F(:,1:ny-4,:,:) + F(:,2:ny-3,:,:))&
    !+ dZ * (-G(:,:,1:nz-4,:) + G(:,:,2:nz-3,:))
    !Q3(i:i+ni-1,j:j+nj-1,k:k+nk-1,:,blockIdx%x) = 0.75d0 * Q(i:i+ni-1,j:j+nj-1,k:k+nk-1,:,blockIdx%x)&
    !+ 0.25d0 * Q2(i:i+ni-1,j:j+nj-1,k:k+nk-1,:,blockIdx%x) - 0.25d0 * dt * R
  end subroutine calc_step2
  
  attributes(global) subroutine calc_step3(nx,ny,nz,ni,nj,nk,dX,dY,dZ,dt,gamma,mu,kappa,Cp,Q3,Q,Q_d)
    integer, intent(in), value :: nx, ny, nz, ni, nj, nk
    real(8), intent(in), value :: dX, dY, dZ, dt, gamma, mu, kappa, Cp
    real(8), intent(in), dimension(nx,ny,nz,5,4) :: Q3
    real(8), intent(out), dimension(nx,ny,nz,5,4) :: Q
    integer :: i, j, k
    real(8), intent(inout) :: Q_d(ni+4,nj+4,nk+4,5)
    k = 3 + nk * mod(int(threadIdx%x/64),8) !k = 3, nz-nk-1, nk :3 <= k <= nz-2
    j = 3 + nj * mod(int(threadIdx%x/8),8)  !j = 3, ny-nj-1, nj :3 <= j <= ny-2
    i = 3 + ni * mod(threadIdx%x,8)         !i = 3, nx-ni-1, ni :3 <= i <= nx-2
    print *, i, j, k
    Q_d = Q3(i-2:i+ni+1,j-2:j+nj+1,k-2:k+nk+1,:,blockIdx%x)
    !call calc_EFG(ni+4,nj+4,nk+4,dX,dY,dZ,gamma,mu,kappa,Cp,Q_d,R)
    !Q(i:i+ni-1,j:j+nj-1,k:k+nk-1,:,blockIdx%x) = (Q(i:i+ni-1,j:j+nj-1,k:k+nk-1,:,blockIdx%x)&
    !+ 2.0d0 * Q3(i:i+ni-1,j:j+nj-1,k:k+nk-1,:,blockIdx%x) - 2.0d0 * dt * R) / 3.d0
  end subroutine calc_step3

  attributes(global) subroutine calc_quantity(nx,ny,nz,ni,nj,nk,gamma,Q,rho,u,v,w,p)
    integer, value :: nx, ny, nz, ni, nj, nk
    real(8), value :: gamma
    real(8) Q(nx,ny,nz,5,4)
    real(8), dimension(nx,ny,nz,4) :: rho, u, v, w, p
    integer :: i, j, k
    k = 3 + nk * mod(int(threadIdx%x/64),8) !k = 3, nz-nk-1, nk :3 <= k <= nz-2
    j = 3 + nj * mod(int(threadIdx%x/8),8)  !j = 3, ny-nj-1, nj :3 <= j <= ny-2
    i = 3 + ni * mod(threadIdx%x,8)         !i = 3, nx-ni-1, ni :3 <= i <= nx-2
    print *, i, j, k
    rho(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x) = Q(i:i+ni-1,j:j+nj-1,k:k+nk-1,1,blockIdx%x)
    u(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x) = Q(i:i+ni-1,j:j+nj-1,k:k+nk-1,2,blockIdx%x) / rho(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x)
    v(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x) = Q(i:i+ni-1,j:j+nj-1,k:k+nk-1,3,blockIdx%x) / rho(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x)
    w(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x) = Q(i:i+ni-1,j:j+nj-1,k:k+nk-1,4,blockIdx%x) / rho(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x)
    p(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x) = (gamma - 1.d0) * (Q(i:i+ni-1,j:j+nj-1,k:k+nk-1,5,blockIdx%x)-0.5d0*&
    rho(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x)*(u(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x)**2+v(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x)**2+w(i:i+ni-1,j:j+nj-1,k:k+nk-1,blockIdx%x)**2))
  end subroutine calc_quantity
end module calc_steps
