subroutine leastSquare(verbose, nobs, nfix, id, env, y, effects)
  use constants
  use quickSort
  implicit none
  !! ================ variable definitions  ================ !!
  logical, intent(in) :: verbose
  integer, intent(in) :: nobs, nfix
  integer, dimension(nobs), intent(in) :: id
  real(KINDR), dimension(nobs), intent(in) :: y, env
  real(KINDR), dimension(nfix), intent(out) :: effects

  character :: uplo, trans
  integer :: info
  integer :: i, j

  real(KINDR) :: val1, val2

  integer, dimension(:), allocatable, save :: tempInd
  real(KINDR), dimension(:,:), allocatable, save :: x , XtX
  real(KINDR), dimension(:), allocatable, save :: Xty
  real(KINDR), dimension(:), allocatable, save :: temp, ipiv
  
  real(KINDR), dimension(:), allocatable, save :: work
  integer :: lwork
  if (.not.allocated(tempInd)) allocate(tempInd(nobs), temp(nobs))

  ! counting distinct values in the environment (nfix)
  call sortrx(nobs, env, tempInd)
  i = 1
  j = 1
  temp(j) = env(tempInd(i))
  do i = 2, nobs
     if (env(tempInd(i)) .ne. env(tempInd(i-1))) then
        j = j + 1
        temp(j) = env(tempInd(i))
     end if
  end do
  if (nfix .ne. j) then 
     write(STDERR, '(a)') "Error:"
     write(STDERR, *) " more/less levels than number of fixed effects"
     write(STDERR, *) " nfix:", nfix, "; counted:", j
     stop 2
  end if
  if (verbose) write(STDOUT, *) " number of fixed effects are ", nfix

  ! now allocation
  if (.not.allocated(x)) then
     allocate(x(nobs, nfix), XtX(nfix, nfix), xty(nfix))
  end if
  
  X(1:nobs, 1:nfix) = ZERO
  ! todo: use sparse matrices: x is very sparse
  do i = 1, nobs
     do j = 1, nfix
        if (temp(j) == env(i)) then
           X(i, j) = ONE
           exit
        end if
     end do
  end do
  
  ! todo: again sparse matrices must be used here
  ! or methods specific to this: look at Nill's book
  uplo = 'u'
  trans = 't'
  val1 = ONE
  val2 = ZERO
  call dsyrk(uplo, trans, nfix, nobs, val1, X, nobs, val2, XtX, nfix)
  if (verbose) write(STDOUT, *) " dsyrk done building x'x"
  do i = 1, nfix
     do j = 1, (i - 1)
        XtX(i, j) = XtX(j, i)
     end do
  end do

  ! todo: again sparse matrices must be used here
  i = 1
  trans = 't'
  val1 = ONE
  val2 = ZERO
  call dgemv(trans, nobs, nfix, val1, X, nobs, y, i, val2, xty, i)
  if (verbose) write(STDOUT, *) " dgemv done with x'y"

  ! solving the system
  uplo = 'u'
  trans = 't'
  val1 = ONE
  val2 = ZERO
  i = 1
  lwork = (nfix * nfix + nfix) * int(nfix  /2)
  if (.not.allocated(work)) allocate(work(lwork), ipiv(nfix))
  call dsysv(uplo, nfix, i, XtX, nfix, ipiv, xty, nfix, work, lwork, info)
  if (info .ne. 0) then
     write(STDERR, '(a)') "Error:"
     if (info .lt. 0) then
        write(STDERR, *) " illegal value at ", -info
     else
        write(STDERR, *) " matrix is not positive definite"
     end if
     stop 2
  else
     if (verbose) then
        write(STDOUT, *) " dsysv solved the system"
        write(STDOUT, *) " work(1), lwork:", work(1), lwork
     end if
  end if
  
  effects(1:nfix) = xty(1:nfix)
end subroutine leastSquare
