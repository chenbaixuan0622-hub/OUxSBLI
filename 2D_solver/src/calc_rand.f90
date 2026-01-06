module calc_rand
  use curand
  use curand_device
  implicit none
contains
  subroutine make_rand(kmax, gen, d_rand)
    integer, intent(in), value           :: kmax
    type(curandGenerator), intent(inout) :: gen
    real(4), intent(out), device         :: d_rand(kmax,kmax)
    integer :: istat, i, j, Nr
    Nr = kmax * kmax
    istat = curandGenerateUniform(gen, d_rand, Nr)
    !$cuf kernel do (2)<<<*,*>>>
    do j = 1, kmax
      do i = 1, kmax
        d_rand(i,j) = 2.d0 * d_rand(i,j) - 1.d0
    enddo;enddo
  end subroutine make_rand


  subroutine init_state(nx, ny, state)
    integer, intent(in) :: nx, ny
    type(curandStateXORWOW), intent(inout), allocatable, device :: state(:)
    integer(8), value :: seed = 12345_8
    integer i
    allocate(state(nx*ny))
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx*ny
      call curand_init(seed, int(i, kind=8), 0_8, state(i))
    enddo
  end subroutine init_state


  subroutine update_state(nx, ny, state)
    integer, intent(in) :: nx, ny
    type(curandStateXORWOW), intent(inout), device :: state(nx*ny)
    real(4) hoge
    integer i
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx*ny
      hoge = curand_uniform(state(i))
    enddo
  end subroutine update_state


  attributes(device) function Z_tilde(state) result(Zt)
    type(curandStateXORWOW), intent(in) :: state
    type(curandStateXORWOW) tmp
    real(8) Zt(3)
    tmp   = state
    Zt(1) = sqrt(2.d0) * (2.d0 * curand_uniform(tmp) - 1.d0)
    Zt(2) = sqrt(2.d0) * (2.d0 * curand_uniform(tmp) - 1.d0)
    Zt(3) = sqrt(2.d0) * (2.d0 * curand_uniform(tmp) - 1.d0)
  end function Z_tilde


  attributes(device) function Zq_x(state) result(Zq)
    type(curandStateXORWOW), intent(in) :: state
    type(curandStateXORWOW) tmp
    real(8) Zq
    tmp = state
    Zq  = 2.d0 * curand_uniform(tmp) - 1.d0
  end function Zq_x


  attributes(device) function Zq_y(state) result(Zq)
    type(curandStateXORWOW), intent(in) :: state
    type(curandStateXORWOW) tmp
    real(8) Zq
    tmp = state
    Zq  = 2.d0 * curand_uniform(tmp) - 1.d0
  end function Zq_y
end module calc_rand

