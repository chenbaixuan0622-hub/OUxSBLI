  !> KEEP (Kinetic-Energy-and-Entropy-Preserving) flux, 2-point face average.
  !> Direct 1D specialization of the 2D/3D KEEP2 kernel: no transverse
  !> momentum component, and the face normal is always +1, so the `Normal`
  !> argument collapses away entirely.
  pure attributes(device) function KEEP2(rho, u, p, T) result(F)
    real(4), intent(in), dimension(2) :: rho, u, p, T
    real(4) F(3)
    F(1) = 0.25 * (rho(1) + rho(2)) * (u(1) + u(2))
    F(2) = 0.5 * (F(1) * (u(1) + u(2)) + (p(1) + p(2)))
    F(3) = F(1) * 0.5 * (T(1) + T(2)) * R_over_gamma_1 ! internal energy
    F(3) = F(3) + 0.5 * (u(1) * p(2) + u(2) * p(1))     ! pressure diffusion
    F(3) = F(3) + 0.5 * F(1) * (u(1) * u(2))            ! kinetic energy
  end function KEEP2
