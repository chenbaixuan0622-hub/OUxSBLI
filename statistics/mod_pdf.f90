module mod_pdf
  implicit none
  type pdf
    integer               :: s    ! size
    real(8), allocatable  :: d(:) ! density
    real(8), allocatable  :: c(:) ! bin center
    real(8)               :: w    ! width
  end type pdf
contains
  subroutine setPDF(n, data, p, min, max)
    integer, intent(in)           :: n
    real(8), intent(in)           :: data(:)
    type(pdf), intent(out)        :: p
    real(8), intent(in), optional :: min, max
    real(8) :: min_val, max_val
    real(8), allocatable :: bin_edges(:)
    integer :: i, bin, num_data

    num_data = size(data)
    if (num_data == 0) then
      print *, "Error: Input data is empty."
      return
    end if

    ! Determine the minimum and maximum values of the data
    if (present(min) .and. present(max)) then
      min_val = min
      max_val = max
    else
      min_val = minval(data)
      max_val = maxval(data)
    endif

    p%s = n
    allocate(p%d(n), p%c(n))
    allocate(bin_edges(n + 1))

    ! Determine the bin width and bin edges
    p%w = (max_val - min_val) / dble(n)
    do i = 0, n
      bin_edges(i + 1) = min_val + dble(i) * p%w
    end do

    ! Initialize pdf
    p%d(:) = 0.d0

    ! Count the number of data points in each bin
    do i = 1, num_data
      bin = int((data(i) - min_val) / p%w) + 1 
      if (bin == n + 1) bin = n  ! Correct for data points exactly equal to max_val
      if (1 <= bin .and. bin <= n) then
        p%d(bin) = p%d(bin) + 1.d0
      end if
    end do

    ! Normalize the histogram to form a probability density function
    p%d(:) = p%d(:) / (dble(num_data) * p%w)

    ! Calculate the bin centers
    do i = 1, n
      p%c(i) = (bin_edges(i) + bin_edges(i + 1)) / 2.d0
    end do
    deallocate(bin_edges)
  end subroutine setPDF
end module mod_pdf

