  !> 2nd-order viscous flux (Sutherland's-law momentum diffusion + heat
  !> conduction), 2-point face average. Shared by both the split (calc_Ev2)
  !> and fused (calc_fused_x) kernels so the formula lives in one place.
  pure attributes(device) function VISC2(u, T, mu) result(Fv)
    real(8), intent(in), dimension(2) :: u, T, mu
    real(8) Fv(2) ! Fv(1) = txx (viscous stress), Fv(2) = viscous work + heat flux
    real(8) mudx, mux
    mudx  = 0.5d0 * (mu(1) + mu(2)) / dx
    Fv(2) = Cp_over_Pr * mudx * (-T(1) + T(2))
    mux   = mudx * (-u(1) + u(2))
    Fv(1) = two_third * 2.d0 * mux
    Fv(2) = Fv(2) + 0.5d0 * (u(1) + u(2)) * Fv(1)
  end function VISC2
