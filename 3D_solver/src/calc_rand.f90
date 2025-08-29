module calc_rand
  implicit none
contains
  attributes(device) function lcg_rand(seed) result(r)
    integer(8), intent(in) :: seed
    real(8) :: s, r
    integer(8), parameter :: a = 1103515245_8
    integer(8), parameter :: c = 12345_8
    integer(8), parameter :: m = 2147483647_8
    s = mod(a * seed + c, m)
    r = dble(s) / dble(m)
  end function lcg_rand

  function lcg_rand_cpu(seed) result(r)
    integer(8), intent(in) :: seed
    real(8) :: s, r
    integer(8), parameter :: a = 1103515245_8
    integer(8), parameter :: c = 12345_8
    integer(8), parameter :: m = 2147483647_8
    s = mod(a * seed + c, m)
    r = dble(s) / dble(m)
  end function lcg_rand_cpu
  
  attributes(device) function box_muller(seed1, seed2) result(ans)
    integer(8), intent(in) :: seed1, seed2
    real(8), parameter :: pi = acos(-1.d0)
    real(8) u1, u2, r, t, ans(2)
    u1 = lcg_rand(seed1)
    u2 = lcg_rand(seed2)
    if (u1 <= 0.d0) u1 = 1.d-16
    r  = sqrt(-2.d0 * log(u1))
    t  = 2.d0 * pi * u2
    ans(1) = r * cos(t)
    ans(2) = r * sin(t)
  end function box_muller

  function box_muller_cpu(seed1, seed2) result(ans)
    integer(8), intent(in) :: seed1, seed2
    real(8), parameter :: pi = acos(-1.d0)
    real(8) u1, u2, r, t, ans(2)
    u1 = lcg_rand_cpu(seed1)
    u2 = lcg_rand_cpu(seed2)
    if (u1 <= 0.d0) u1 = 1.d-16
    r  = sqrt(-2.d0 * log(u1))
    t  = 2.d0 * pi * u2
    ans(1) = r * cos(t)
    ans(2) = r * sin(t)
  end function box_muller_cpu

  attributes(device) function Z_tilde(seed) result(Zt)
    integer(8), intent(in), value :: seed
    integer(8) s1, s2, s3, s4
    real(8) Z(6), Zt(6)
    s1 = seed
    s2 = seed + 1_8
    s3 = seed + 2_8
    s4 = seed + 3_8
    Z(1:2) = box_muller(s1, s2)
    Z(3:4) = box_muller(s2, s3)
    Z(5:6) = box_muller(s3, s4)
    Zt(:)  = sqrt(2.d0) * Z(:)
  end function Z_tilde

  attributes(device) function Zq_x(seed) result(Zq)
    integer(8), intent(in), value :: seed
    integer(8) s1, s2
    real(8) Zq, tmp(2)
    s1  = seed
    s2  = seed + 1_8
    tmp = box_muller(s1, s2)
    Zq  = tmp(1)
  end function Zq_x

  attributes(device) function Zq_y(seed) result(Zq)
    integer(8), intent(in), value :: seed
    integer(8) s1, s2
    real(8) Zq, tmp(2)
    s1  = seed + 1_8
    s2  = seed + 2_8
    tmp = box_muller(s1, s2)
    Zq  = tmp(1)
  end function Zq_y

  attributes(device) function Zq_z(seed) result(Zq)
    integer(8), intent(in), value :: seed
    integer(8) s1, s2
    real(8) Zq, tmp(2)
    s1  = seed + 2_8
    s2  = seed + 3_8
    tmp = box_muller(s1, s2)
    Zq  = tmp(1)
  end function Zq_z

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
  end subroutine init_seed
  
  subroutine init_seed_cpu(nx, ny, nz, seed)
    integer, intent(in)     :: nx, ny, nz
    integer(8), intent(out) :: seed(nx,ny,nz)
    integer i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          seed(i,j,k) = nx * ny * int(k-1, kind=8) + nx * int(j-1, kind=8) + int(i, kind=8)
    enddo;enddo;enddo
  end subroutine init_seed_cpu
 
  subroutine update_seed(nx, ny, nz, seed)
    integer, intent(in), value        :: nx, ny, nz
    integer(8), intent(inout), device :: seed(nx,ny,nz)
    integer(8), parameter :: a = 1103515245_8
    integer(8), parameter :: c = 12345_8
    integer(8), parameter :: m = 2147483647_8
    integer i, j, k
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          seed(i,j,k) = mod(a * seed(i,j,k) + c, m)
    enddo;enddo;enddo
  end subroutine update_seed
  
  subroutine calc_rand_uniform(nx, ny, nz, seed, rand)
    integer, intent(in)       :: nx, ny, nz
    integer(8), intent(inout) :: seed(nx,ny,nz)
    real(8), intent(inout)    :: rand(nx,ny,nz)
    integer i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rand(i,j,k) = lcg_rand_cpu(seed(i,j,k))
    enddo;enddo;enddo
  end subroutine calc_rand_uniform

  subroutine calc_rand_normal(nx, ny, nz, seed, rand)
    integer, intent(in)       :: nx, ny, nz
    integer(8), intent(inout) :: seed(nx,ny,nz)
    real(8), intent(inout)    :: rand(nx,ny,nz)
    real(8) tmp(2)
    integer i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          tmp = box_muller_cpu(seed(i,j,k), seed(i,j,k))
          rand(i,j,k) = tmp(1)
    enddo;enddo;enddo
  end subroutine calc_rand_normal
end module calc_rand


!program test_rand
!  use cudafor
!  use calc_rand
!  implicit none
!  integer, parameter :: nx = 32, ny = 32, nz = 32
!  integer(8), allocatable :: seed(:,:,:)
!  real(8), allocatable    :: rand(:,:,:)
!  real(8) mean, std
!  allocate(seed(nx,ny,nz), rand(nx,ny,nz))
!  call init_seed_cpu(nx, ny, nz, seed)
!  call calc_rand_normal(nx, ny, nz, seed, rand)
!  !call calc_rand_uniform(nx, ny, nz, seed, rand)
!  mean = sum(rand) / size(rand)
!  std  = sqrt(sum(rand**2) / size(rand) - mean**2)
!  print *, "mean=", mean, " std=", std
!  deallocate(seed, rand)
!end program test_rand

