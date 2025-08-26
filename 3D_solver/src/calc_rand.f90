module calc_rand
  implicit none
contains
  attributes(device) function lcg_rand(seed) result(r)
    integer(8), intent(inout) :: seed
    real(8) :: r
    integer(8), parameter :: a = 6364136223846793005_8
    integer(8), parameter :: c = 1_8
    seed = a * seed + c
    r = dble(seed) / dble(huge(0_8) + 1_8)
  end function lcg_rand

  attributes(device) function box_muller4(seed1, seed2) result(ans)
    integer(8), intent(inout) :: seed1, seed2
    real(8), parameter :: pi = acos(-1.d0)
    real(8) u1, u2, r, t, ans(4)
    u1 = lcg_rand(seed1)
    u2 = lcg_rand(seed2)
    if (u1 <= 0.d0) u1 = 1.d-16
    r  = sqrt(-2.d0 * log(u1))
    t  = 2.d0 * pi * u2
    ans(1) = r * cos(t)
    ans(2) = r * sin(t)
    u1 = lcg_rand(seed1)
    u2 = lcg_rand(seed2)
    if (u1 <= 0.d0) u1 = 1.d-16
    r  = sqrt(-2.d0 * log(u1))
    t  = 2.d0 * pi * u2
    ans(3) = r * cos(t)
    ans(4) = r * sin(t)
  end function box_muller4

  subroutine init_seed(nx, ny, nz, seed)
    integer, intent(in)             :: nx, ny, nz
    integer(8), intent(out), device :: seed(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          seed(i,j,k) = nx * ny * int(k-1, kind=8) + nx * int(j-1, kind=8) + int(i, kind=8)
    enddo;enddo;enddo
  end subroutine
end module calc_rand

