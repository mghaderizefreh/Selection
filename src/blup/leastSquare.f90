subroutine leastSquare(verbose, nobs, nfix, id, env, y, effects, tempInd, temp, info)
  use constants, only : KINDR, alloc1D, alloc2D, STDOUT, STDERR, ONE, ZERO
  use quickSort, only : sortrx
  implicit none
  !! ================ variable definitions  ================ !!
  logical, intent(in) :: verbose
  integer, intent(in) :: nobs, nfix
  integer, dimension(nobs), intent(in) :: id
  real(KINDR), dimension(1:nobs), intent(in) :: y, env
  real(KINDR), dimension(1:nfix), intent(inout) :: effects
  integer, dimension(1:nobs), intent(inout) :: tempInd
  real(KINDR), dimension(1:nobs), intent(inout) :: temp
  character :: uplo, trans
  integer, intent(inout) :: info
  integer :: i, j

  real(KINDR) :: val1, val2
  real(KINDR), dimension(:,:), allocatable :: X
  real(KINDR), dimension(1:nfix,1:nfix) :: XtX
  real(KINDR), dimension(1:nfix) :: ipiv
  real(KINDR), dimension(1) :: junk
  real(KINDR), dimension(:), allocatable :: work
  integer :: lwork

  external :: DSYRK, DGEMV, DSYSV

  ! may need id later
  if (id(1) > 0) then
  end if
  
  call alloc2D(X,nobs, nfix, "x", "leastSquare")
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
     write(STDERR, '(a)') "Error in leastSquare:"
     write(STDERR, *) " more/less levels than number of fixed effects"
     write(STDERR, *) " nfix:", nfix, "; counted:", j
     stop 2
  end if
  if (verbose) write(STDOUT, *) " number of fixed effects are ", nfix

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
  do i = 1, nfix ! todo: do I really need this loop?
     do j = 1, (i - 1)
        XtX(i, j) = XtX(j, i)
     end do
  end do

  ! todo: again sparse matrices must be used here
  i = 1
  trans = 't'
  val1 = ONE
  val2 = ZERO
  call dgemv(trans, nobs, nfix, val1, X, nobs, y, i, val2, effects, i)
  if (verbose) write(STDOUT, *) " dgemv done with x'y"

  ! solving the system
  uplo = 'u'
  trans = 't'
  val1 = ONE
  val2 = ZERO
  i = 1
  ! to get proper value for lwork first it is passed as -1
  lwork = -1
  call dsysv(uplo, nfix, i, XtX, nfix, ipiv, effects, nfix, junk, lwork, info)
  if (info /= 0) then
     write(STDERR, '(A)') "Error in leastSquare:"
     write(STDERR, *) "Could not get optimum lwork size. Exiting"
     stop 2
  end if
  lwork = int(junk(1))
  call alloc1D(work, lwork, "work", "leastSquare")
  call dsysv(uplo, nfix, i, XtX, nfix, ipiv, effects, nfix, work, lwork, info)
  if (info .ne. 0) then
     write(STDERR, '(a)') "Error in leastSquare:"
     if (info .lt. 0) then
        write(STDERR, *) " illegal value at ", -info
     else
        write(STDERR, *) " matrix is not positive definite. Exiting"
     end if
     stop 2
  else
     if (verbose) then
        write(STDOUT, *) " dsysv solved the system"
        write(STDOUT, *) " work(1), lwork:", work(1), lwork
     end if
  end if

  deallocate(X)
  deallocate(work)
end subroutine leastSquare
