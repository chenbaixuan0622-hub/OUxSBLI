module mod_globals
  implicit none
  integer, parameter :: accuracy = 2 ! only 2nd-order accuracy is available
  integer, parameter :: offset = accuracy / 2
  integer, parameter :: id_visc = 1
  !!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme ! 1  KEEP !
  !           ! 2  SLAU !
  !!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: id_scheme = 2
end module mod_globals

