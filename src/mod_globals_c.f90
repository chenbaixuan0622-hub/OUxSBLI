module mod_globals_c
  use iso_c_binding
  use mod_globals, only : nx, ny, nz, nt, np, Lx, Ly, Lz, dt, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use set_coordinate
  implicit none
contains
  function get_Nx_c() result(ans) bind(c, name="get_Nx_c")
    integer(c_int) :: ans
    ans = nx
  end function get_Nx_c
  
  function get_Ny_c() result(ans) bind(c, name="get_Ny_c")
    integer(c_int) :: ans
    ans = ny
  end function get_Ny_c
  
  function get_Nz_c() result(ans) bind(c, name="get_Nz_c")
    integer(c_int) :: ans
    ans = nz
  end function get_Nz_c

  function get_Nt_c() result(ans) bind(c, name="get_Nt_c")
    integer(c_int) :: ans
    ans = nt
  end function get_Nt_c

  function get_Np_c() result(ans) bind(c, name="get_Np_c")
    integer(c_int) :: ans
    ans = np
  end function get_Np_c

  function get_Lx_c() result(ans) bind(c, name="get_Lx_c")
    real(c_double) :: ans
    ans = Lx
  end function get_Lx_c

  function get_Ly_c() result(ans) bind(c, name="get_Ly_c")
    real(c_double) :: ans
    ans = Ly
  end function get_Ly_c

  function get_Lz_c() result(ans) bind(c, name="get_Lz_c")
    real(c_double) :: ans
    ans = Lz
  end function get_Lz_c
  
  function get_dt_c() result(ans) bind(c, name="get_dt_c")
    real(c_double) :: ans
    ans = dt
  end function get_dt_c
 
  subroutine set_block_c(nx, ny, nz) bind(c, name="set_block_c")
    integer, intent(in) :: nx, ny, nz
    call set_block(nx, ny, nz, threads, threadsE, threadsEv, threadsF, threadsFv, threadsG, threadsGv, &
                   blocks, blocksE, blocksEv, blocksF, blocksFv, blocksG, blocksGv)
  end subroutine set_block_c
end module mod_globals_c

