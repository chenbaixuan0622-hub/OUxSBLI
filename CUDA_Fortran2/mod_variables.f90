module mod_variables
  use cudafor
  implicit none
  real(8), dimension(3), device :: rho_d, u_d, v_d, w_d, p_d, P_keep, IE, E_total
  real(8), dimension(3), device :: RhoU, RhoUU, RhoUV, RhoUW, RhoUIE, RhoUKE, UP, RhoUU_P
  real(8), dimension(3), device :: RhoV, RhoVU, RhoVV, RhoVW, RhoVIE, RhoVKE, VP, RhoVV_P
  real(8), dimension(3), device :: RhoW, RhoWU, RhoWV, RhoWW, RhoWIE, RhoWKE, WP, RhoWW_P
end module mod_variables