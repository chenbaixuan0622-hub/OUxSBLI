module mod_statistics
  use mod_pdf
  implicit none
contains
  subroutine normal_distribution(n,mean,stddev,normal_samples)
    integer, intent(in)   :: n
    real(8), intent(in)   :: mean
    real(8), intent(in)   :: stddev
    real(8), intent(out)  :: normal_samples(n)
    real(8) :: u1, u2, z0, z1
    real(8) :: pi = 4.d0 * atan(1.d0)
    integer :: i

    ! seed te random number generator
    call random_seed()

    do i = 1, n, 2
      ! Generate two random numbers between 0 and 1
      call random_number(u1)
      call random_number(u2)

      ! Apply Boc-Muller transform
      z0 = sqrt(-2.d0 * log(u1)) * cos(2.d0 * pi * u2)
      z1 = sqrt(-2.d0 * log(u1)) * sin(2.d0 * pi * u2)

      ! Store the generated normal random numbers
      normal_samples(i) = mean + stddev * z0
      if (i + 1 <= n) then
        normal_samples(i+1) = mean + stddev * z1
      endif
    enddo
  end subroutine normal_distribution

  subroutine calc_KLD(p,q,KLD)
    type(pdf), intent(in) :: p, q
    real(8), intent(out)  :: KLD
    integer i
    ! check
    if (p%s /= q%s) then
      write(*,*) "invalid input"
      stop
    elseif (p%w /= q%w) then
      write(*,*) "invalid input"
      stop
    endif
    ! set init value KLD = 0
    KLD = 0.d0
    do i = 1, p%s
      if (p%d(i) /= 0.d0 .and. q%d(i) /= 0.d0) then
        KLD = KLD + p%d(i) * log(p%d(i) / q%d(i)) * p%w
      endif
    enddo
  end subroutine calc_KLD
end module mod_statistics

