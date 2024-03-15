module calc_slau
  use calc_MUSCL
  implicit none
contains
  attributes(device) function energy(gamma,p,rho,u,v) result(e)
    real(8), intent(in), value :: gamma, p, rho, u , v
    real(8) :: e
    e = p / (gamma - 1.d0) + 0.5d0 * rho * (u ** 2 + v ** 2)
  end function energy

  attributes(device) function enthalpy(e,p,rho) result(h)
    real(8), intent(in), value :: e, p, rho
    real(8) :: h
    h = (e + p) / rho
  end function enthalpy

  attributes(device) function speed_of_sound(gamma,p,rho) result(c)
    real(8), intent(in), value :: gamma, p, rho
    real(8) :: c
    c = sqrt(gamma * p / rho)
  end function speed_of_sound

  attributes(device) function flux_SLAU(dim,gamma,rhol,rhor,pl,pr,Vl,Vr,Normal) result(Flux)
    integer, intent(in), value :: dim ! x:1, y:2
    real(8), intent(in), value :: gamma, rhol, rhor, pl, pr
    real(8), intent(in), dimension(2), device :: Vl, Vr
    real(8), intent(in), dimension(4), device :: Normal
    real(8) :: Flux(4)
    real(8) el, er, hl, hr, cl, cr, c
    real(8) vn, V_p, V_m, V_bar, V_bar_p, V_bar_m, M_p, M_m, M, x, g, dp, mass, beta_p, beta_m, Pressure
    real(8), dimension(4) :: phil, phir
    el = energy(gamma,pl,rhol,Vl(1),Vl(2))
    er = energy(gamma,pr,rhor,Vr(1),Vr(2))
    hl = enthalpy(el,pl,rhol)
    hr = enthalpy(er,pr,rhor)
    cl = speed_of_sound(gamma,pl,rhol)
    cr = speed_of_sound(gamma,pr,rhor)
    c = 0.5d0 * (cl + cr)
    ! What is vn?
    vn = 0.d0
    V_p = Vl(dim) - vn
    V_m = Vr(dim) - vn
    M_p = V_p / c
    M_m = V_m / c
    M = min(1.d0, sqrt(0.5d0 * (M_p ** 2 + M_m ** 2)))
    x = (1.d0 - M) ** 2
    g = -max(min(M_p, 0.d0), -1.d0) * min(max(M_m, 0.d0), 1.d0)
    V_bar = (rhol * abs(V_p) + rhor * abs(V_m)) / (rhol + rhor)
    V_bar_p = abs((1.d0 - g) * V_bar + g * V_p)
    V_bar_m = abs((1.d0 - g) * V_bar + g * V_m)
    dp = pr - pl
    mass = 0.5d0 * (rhol * (Vl(dim) + V_bar_p) + rhor * (Vr(dim) - V_bar_m) - x * dp / c)
    if (abs(M_p) < 1.d0) then
      beta_p = 0.25d0 * (2.d0 - M_p) * (M_p + 1.d0) ** 2
    else
      beta_p = 0.5d0 * (1.d0 + sign(1.d0, M_p))
    endif
    if (abs(M_m) < 1.d0) then
      beta_m = 0.25d0 * (2.d0 + M_m) * (M_m - 1.d0) ** 2
    else
      beta_m = 0.5d0 * (1.d0 + sign(1.d0, -M_m))
    endif
    Pressure = 0.5d0 * (pl + pr + (beta_p - beta_m) * (pl - pr) + (1.d0 - x) * (beta_p + beta_m - 1.d0) * (pl + pr))
    phil(:) = (/1.d0, Vl(1), Vr(2), hl/)
    phir(:) = (/1.d0, Vr(1), Vr(2), hr/)
    Flux(:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + Pressure * Normal(:)
  end function flux_SLAU

  attributes(global) subroutine calc_E(nx, ny, gamma, k, b, rho, u, v, p, E)
    use mod_globals, only : accuracy
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma, k, b
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), dimension(nx-1,ny-accuracy,4), device :: E
    integer i, j
    integer :: offset = accuracy / 2
    real(8) rhol, rhor, pl, pr
    real(8), dimension(2) :: Vl, Vr
    real(8) :: Normal(4) = (/0.d0, 1.d0, 0.d0, 0.d0/)
    real(8), dimension(4) :: d1, d2, d3, Q1, Q2, Ql, Qr
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    
    ! left edge
    if (i == 1) then
      d1(:) = 0.d0
    else
      d1(:) = (/-rho(i-1,j) + rho(i,j), -u(i-1,j) + u(i,j), -v(i-1,j) + v(i,j), -p(i-1,j) + p(i,j)/)
    endif
        
    d2(:) = (/-rho(i,j) + rho(i+1,j), -u(i,j) + u(i+1,j), -v(i,j) + v(i+1,j), -p(i,j) + p(i+1,j)/)
        
    ! right edge
    if (i == nx-1) then
      d3(:) = 0.d0
    else 
      d3(:) = (/-rho(i+1,j) + rho(i+2,j), -u(i+1,j) + u(i+2,j), -v(i+1,j) + v(i+2,j), -p(i+1,j) + p(i+2,j)/)
    endif

    Q1(:) = (/rho(i,j), u(i,j), v(i,j), p(i,j)/)
    Q2(:) = (/rho(i+1,j), u(i+1,j), v(i+1,j), p(i+1,j)/)
        
    call MUSCL(4,k,b,Q1(:),Q2(:),d1(:),d2(:),d3(:),Ql(:),Qr(:))
    rhol = Ql(1)
    rhor = Qr(1)
    Vl(:) = (/Ql(2), Ql(3)/)
    Vr(:) = (/Qr(2), Qr(3)/)
    pl = Ql(4)
    pr = Qr(4)
    E(i,j-offset,:) = flux_SLAU(1,gamma,rhol,rhor,pl,pr,Vl,Vr,Normal)
  end subroutine calc_E

  attributes(global) subroutine calc_F(nx, ny, gamma, k, b, rho, u, v, p, F)
    use mod_globals, only : accuracy
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma, k, b
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), dimension(nx-accuracy,ny-1,4), device :: F
    integer i, j
    integer :: offset = accuracy / 2
    real(8) rhol, rhor, pl, pr
    real(8), dimension(2) :: Vl, Vr
    real(8) :: Normal(4) = (/0.d0, 0.d0, 1.d0, 0.d0/)
    real(8), dimension(4) :: d1, d2, d3, Q1, Q2, Ql, Qr
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y

    ! left edge
    if (j == 1) then
      d1(:) = 0.d0
    else
      d1(:) = (/-rho(i,j-1) + rho(i,j), -u(i,j-1) + u(i,j), -v(i,j-1) + v(i,j), -p(i,j-1) + p(i,j)/)
    endif
        
    d2(:) = (/-rho(i,j) + rho(i,j+1), -u(i,j) + u(i,j+1), -v(i,j) + v(i,j+1), -p(i,j) + p(i,j+1)/)
        
    ! right edge
    if (j == ny-1) then
      d3(:) = 0.d0
    else 
      d3(:) = (/-rho(i,j+1) + rho(i,j+2), -u(i,j+1) + u(i,j+2), -v(i,j+1) + v(i,j+2), -p(i,j+1) + p(i,j+2)/)
    endif

    Q1(:) = (/rho(i,j), u(i,j), v(i,j), p(i,j)/)
    Q2(:) = (/rho(i,j+1), u(i,j+1), v(i,j+1), p(i,j+1)/)
    
    call MUSCL(4,k,b,Q1(:),Q2(:),d1(:),d2(:),d3(:),Ql(:),Qr(:))
    rhol = Ql(1)
    rhor = Qr(1)
    Vl(:) = (/Ql(2), Ql(3)/)
    Vr(:) = (/Qr(2), Qr(3)/)
    pl = Ql(4)
    pr = Qr(4)
    F(i-offset,j,:) = flux_SLAU(2,gamma,rhol,rhor,pl,pr,Vl,Vr,Normal)
  end subroutine calc_F
end module calc_slau

