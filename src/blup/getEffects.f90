!This subroutine calculates the fixed and random effects based on the given variances (theta). 
! In writing the random effects First comes slope effects; then intercept effects; then individual slope residual
!written by Masoud Ghaderi
subroutine getEffects(nobs, maxid, nfix, nvar, nran, theta, Gmatrix, Vhat,&
     Py, y, X, id, fixeff, raneff, verbose)
  use constants
  use global_module
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
  type (doublePre_Array), dimension(1:nran), intent(out) :: raneff

  real(KINDR), dimension(:), allocatable :: temp
  integer :: i, j
  real(KINDR) :: val1, s1, s2
  type (doublePre_Array), dimension(1:nran) :: theZPy

  external :: dgemm
  if (verbose) write(STDOUT, *) " Inside getEffects"
  ! fixed effects
  call dgemm('n', 'n', nfix, 1, nobs, ONE, Vhat, nfix, y, nobs, ZERO, fixeff, nfix)
  if (verbose) write(STDOUT, *) "  Fixed effects estimated"
  
  if (nran == 1) then
     allocate(theZPy(1)%level(maxid))
     theZPy(1)%level(1:maxid) = ZERO
     raneff(1)%level(1:maxid) = ZERO
     do i = 1, nobs
        j = id(i)
        theZPy(1)%level(j) = theZPy(1)%level(j) + Py(i)
     end do
     theZPy(1)%level(1:maxid) = theZPy(1)%level(1:maxid) * theta(1)
     do i = 1, maxid
        do j = 1, i
           val1 = Gmatrix((i * (i - 1)) / 2 + j)
           raneff(1)%level(i) = raneff(1)%level(i) + theZPy(1)%level(j) * val1
           if (i .ne. j) then
              raneff(1)%level(j) = raneff(1)%level(j) + theZPy(1)%level(i) * val1
           end if
        end do
     end do
     return
  end if
  ! allocation
  allocate(theZPy(1)%level(maxid)) ! slope effect (genetic)
  allocate(theZPy(2)%level(maxid)) ! intercept effect (genetic)
  allocate(theZPy(3)%level(nobs))   ! environment slope effect (diagonal)

  ! initialisation
  do i = 1, 3
     theZPy(i)%level(:) = ZERO
     raneff(i)%level(:) = ZERO
  end do

  if (verbose) write(STDOUT, *) "  initialisation done for random effects"
  allocate(temp(maxid))

  ! random effects
  do i = 1, nobs
     j = id(i)
     val1 = Py(i) * X(i,1)
     theZPy(1)%level(j) = theZPy(1)%level(j) + val1
     theZPy(2)%level(j) = theZPy(2)%level(j) + Py(i)
     raneff(3)%level(i) = X(i,1) * Py(i) * theta(3)
  end do
  
  if (nvar == 4) then
     s1 = theta(4) / theta(1)
     s2 = theta(4) / theta(2)
     temp(1:maxid) = theZPy(1)%level(1:maxid)
     theZPy(1)%level(1:maxid) = theZPy(1)%level(1:maxid) * theta(1) + theZPy(2)%level(1:maxid) * theta(4)
     theZPy(2)%level(1:maxid) = theZPy(2)%level(1:maxid) * theta(2) + temp(1:maxid)            * theta(4)
  else
     theZPy(1)%level(1:maxid) = theZPy(1)%level(1:maxid) * theta(1)
     theZPy(2)%level(1:maxid) = theZPy(2)%level(1:maxid) * theta(2)
  end if

  do i = 1, maxid
     do j = 1, i
        val1 = Gmatrix((i * (i - 1)) / 2 + j)
        raneff(1)%level(i) = raneff(1)%level(i) + theZPy(1)%level(j) * val1
        raneff(2)%level(i) = raneff(2)%level(i) + theZPy(2)%level(j) * val1
        if (i .ne. j) then
           raneff(1)%level(j) = raneff(1)%level(j) + theZPy(1)%level(i) * val1
           raneff(2)%level(j) = raneff(2)%level(j) + theZPy(2)%level(i) * val1
        end if
     end do
  end do
  
end subroutine getEffects
