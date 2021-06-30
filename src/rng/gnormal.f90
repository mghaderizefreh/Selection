subroutine gnormal(mean, cov, dim, N, output, seed)
  ! this subroutine produces normal distribution using boxMuller transformation


  ! Inputs:
  !     `dim`    : dimension of the distribution (integer)
  !     `mean`   : mean of the distribution (1D real(KINDR) array of shape `dim`)
  !     `cov`    : covariance matrix (2D real(KINDR) array of shape `dim` x `dim`)
  !     `N`   : number of samples (integer)
  !     `output` : (2D real(KINDR) array of shape `dim` x `N`) output
  !     `nseed`  : seed for random stream (integer)
  ! matrix `cov`. `N` is the length of the array and `output` is a real(KINDR)
  ! array of shape `dim` x `output`.
  !
  ! written by Masoud Ghaderi Zefreh
  ! first revision : 1 December 2020
  
  use constants

  implicit none
  integer, intent(in) :: dim, N
  real(KINDR), dimension(1:dim), intent(in) :: mean
  real(KINDR), dimension(1:dim, 1:dim), intent(in) :: cov
  real(KINDR), dimension(1:N, 1:dim), intent(inout) :: output
  integer, intent(in), dimension(:), optional :: seed

  integer :: dim2, i, j
  real(KINDR) :: sd, val
  real(KINDR), dimension(dim, dim) :: covcopy
  real(KINDR), dimension(:,:), allocatable :: uniform

  external :: dpotrf, dgemm, normdev

  covcopy = cov

  if (present(seed)) then
     call random_seed(put = seed)
  end if

  if (dim == 1) then
     sd = sqrt(cov(1,1))
     do i = 1, N
        call normdev(val)
        output(i,1) = mean(1) + sd * val
     end do
     return
  else

     dim2 = dim
     if (mod(dim, 2) .eq. 1) dim2 = dim2 + 1

     ! making uniform distribution
     call alloc2D(uniform, N, dim2, "uniform", "gnormal")
     call random_number(uniform)

     ! making independent pairs of normally distributed data
     do i = 1, dim, 2
        j = i + 1
        output(:,i) = sqrt(-2 * log(uniform(:,i))) * cos(2 * PI * uniform(:,j))
     end do
     do j = 2, dim, 2
        i = j - 1
        output(:,j) = sqrt(-2 * log(uniform(:,i))) * sin(2 * PI * uniform(:,j))
     end do

     if (dim .ne. dim2)then
        deallocate(uniform)
        call alloc2D(uniform, N,dim, "uniform", "gnormal")
     end if
     uniform = output

     ! cholesky factorisation
     call dpotrf('L',dim,covcopy,dim,i)
     if (i.gt.0) then
        write(STDERR, *) ' Cholesky factorisation failed!'
        write(STDERR, *) ' Covariance matrix is not positive definite'
        write(STDERR, *) ' Exiting...'
        stop 2
     elseif (i.lt.0) then
        write(STDERR, *) "some other shit happened"
        stop 2
     end if

     ! the cholesky factorisation routine ?potrf(with 'L') returns the lower
     ! part of the cholesky. The upper part must be zero.
     do i = 1, dim - 1
        do j = i + 1, dim
           covcopy(i,j) = ZERO
        end do
     end do

     ! the output is normal(nxd) = unifrm(nxd) * chol.Transpose
     call dgemm('N','T', N, dim, dim, ONE, uniform, N, covcopy, dim, &
          ZERO, output, N)

     forall( i = 1:N)
        output(i,1:DIM) = mean(1:DIM) + output(i,1:DIM)
     end forall
  end if
  deallocate(uniform)
end subroutine gnormal

subroutine normdev(val)
! making a random varaible form standard normal
! copy of ricardo's code (goto replaced by while)
! Routine 'gasdev' ,Numerical Recipes, page 203
  use constants
  implicit none
  save
  real(KINDR)    ::val,v1,v2,fac, r = 2.
  real(KINDR)    ::gset =0.0
  integer ::iset =0

  if (iset .eq. 0) then
     do while(r.ge.1..or.r.eq.0.)
        call random_number(v1)
        call random_number(v2)
        v1 = 2. * v1 - 1.
        v2 = 2. * v2 - 1.
        r = v1 ** 2 + v2 ** 2
     end do
     fac = sqrt(-2. * log(r) / r)
     gset= v1 * fac
     val = v2 * fac
     iset= 1
  else
     val = gset
     iset = 0
  end if
end subroutine normdev
