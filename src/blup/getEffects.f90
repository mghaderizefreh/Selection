!This subroutine calculates the fixed and random effects based on the given variances (theta). 
! In writing the random effects First comes slope effects; then intercept effects; then individual slope residual
!written by Masoud Ghaderi
subroutine getEffects(nobs, maxid, nfix, nvar, nran, theta, Gmatrix, Vhat,&
     Py, y, X, id, fixeff, raneff, verbose)
  use constants, only : KINDR, Jarr, STDOUT, alloc1D, ONE, ZERO
  implicit none
  logical :: verbose
  integer, intent(in) :: nobs, nfix, nvar, maxid, nran
  real(KINDR), dimension(1:(maxid*(maxid+1))), intent(in) :: Gmatrix
  integer, dimension(1:nobs), intent(in) :: id
  real(KINDR), dimension(1:(nvar+1)), intent(in) :: theta
  real(KINDR), dimension(1:nobs), intent(in) :: Py, y
  real(KINDR), dimension(1:nfix,1:nobs), intent(in) :: Vhat
  real(KINDR), dimension(1:nobs,1:nfix), intent(in) :: X
  real(KINDR), dimension(1:nfix), intent(out) :: fixeff
  type (Jarr), dimension(1:nran), intent(inout) :: raneff

  real(KINDR), dimension(:), allocatable :: temp
  integer :: i, j
  real(KINDR) :: val1
  type (Jarr), dimension(:), allocatable :: theZPy

  external :: dgemm
  if (verbose) write(STDOUT, *) " Inside getEffects"

  ! fixed effects
  call dgemm('n', 'n', nfix, 1, nobs, ONE, Vhat, nfix, y, nobs, ZERO, fixeff, nfix)
  if (verbose) write(STDOUT, *) "  Fixed effects estimated"

  allocate(theZPy(nran))
  if (nran == 1) then
     call alloc1D(theZPy(1)%array, maxid, "theZPy(1)%array", "getEffects")
     theZPy(1)%array(1:maxid) = ZERO
     raneff(1)%array(1:maxid) = ZERO
     do i = 1, nobs
        j = id(i)
        theZPy(1)%array(j) = theZPy(1)%array(j) + Py(i)
     end do
     theZPy(1)%array(1:maxid) = theZPy(1)%array(1:maxid) * theta(1)
     do i = 1, maxid
        do j = 1, i
           val1 = Gmatrix((i * (i - 1)) / 2 + j)
           raneff(1)%array(i) = raneff(1)%array(i) + theZPy(1)%array(j) * val1
           if (i .ne. j) then
              raneff(1)%array(j) = raneff(1)%array(j) + theZPy(1)%array(i) * val1
           end if
        end do
     end do
     return
  end if
  ! allocation
  call alloc1D(theZPy(1)%array,maxid, "theZPy(1)%array", "getEffects")!slope effect (genetic)
  call alloc1D(theZPy(2)%array,maxid, "theZPy(2)%array", "getEffects")!int effect (genetic)
  call alloc1D(theZPy(3)%array,nobs, "theZPy(3)%array", "getEffects") !env slope eff (diag.)

  ! initialisation
  do i = 1, 3
     theZPy(i)%array(:) = ZERO
     raneff(i)%array(:) = ZERO
  end do

  if (verbose) write(STDOUT, *) "  initialisation done for random effects"
  call alloc1D(temp, maxid, "temp", "getEffects")

  ! random effects
  do i = 1, nobs
     j = id(i)
     val1 = Py(i) * X(i,1)
     theZPy(1)%array(j) = theZPy(1)%array(j) + val1
     theZPy(2)%array(j) = theZPy(2)%array(j) + Py(i)
     raneff(3)%array(i) = X(i,1) * Py(i) * theta(3)
  end do

  if (nvar == 4) then
     temp(1:maxid) = theZPy(1)%array(1:maxid)
     theZPy(1)%array(1:maxid) = theZPy(1)%array(1:maxid) * theta(1) + &
        theZPy(2)%array(1:maxid) * theta(4)
     theZPy(2)%array(1:maxid) = theZPy(2)%array(1:maxid) * theta(2) + &
        temp(1:maxid)            * theta(4)
  else
     theZPy(1)%array(1:maxid) = theZPy(1)%array(1:maxid) * theta(1)
     theZPy(2)%array(1:maxid) = theZPy(2)%array(1:maxid) * theta(2)
  end if

  do i = 1, maxid
     do j = 1, i
        val1 = Gmatrix((i * (i - 1)) / 2 + j)
        raneff(1)%array(i) = raneff(1)%array(i) + theZPy(1)%array(j) * val1
        raneff(2)%array(i) = raneff(2)%array(i) + theZPy(2)%array(j) * val1
        if (i .ne. j) then
           raneff(1)%array(j) = raneff(1)%array(j) + theZPy(1)%array(i) * val1
           raneff(2)%array(j) = raneff(2)%array(j) + theZPy(2)%array(i) * val1
        end if
     end do
  end do
end subroutine getEffects
