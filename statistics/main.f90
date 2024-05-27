program main
  use mod_pdf
  use mod_statistics
  implicit none
  integer, parameter    :: n = 100000
  integer, parameter    :: n_pdf = 1000
  real(8)               :: mean1, mean2, stddev1, stddev2, min_val, max_val, KLD
  real(8), dimension(n) :: normal_samples1, normal_samples2
  type(pdf)             :: p, q
  integer i

  mean1 = 1.d0
  mean2 = 0.d0
  stddev1 = 5.d0
  stddev2 = 1.d0

  call normal_distribution(n,mean1,stddev1,normal_samples1)
  call normal_distribution(n,mean2,stddev2,normal_samples2)

  min_val = min(minval(normal_samples1), minval(normal_samples2))
  max_val = max(maxval(normal_samples1), maxval(normal_samples2))

  call setPDF(n_pdf,normal_samples1,p,min_val,max_val)
  call setPDF(n_pdf,normal_samples2,q,min_val,max_val)

  call calc_KLD(p,q,KLD)
  print *, "Kullback Leibler divergence", KLD

  open(10,file='data/pdf1.d', status='replace', action='write')
  do i = 1, p%s
    write(10,*) p%c(i), p%d(i)
  enddo
  close(10)

  open(10,file='data/pdf2.d', status='replace', action='write')
  do i = 1, q%s
    write(10,*) q%c(i), q%d(i)
  enddo
  close(10)
end program main

