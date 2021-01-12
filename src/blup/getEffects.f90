!This subroutine calculates the fixed and random effects based on the given variances (theta). 
!It then creates (or overwrites) three files with these information.
! In writing the random effects First comes slope effects; then intercept effects; then individual slope residual
!written by Masoud Ghaderi
!last modified 25 March 2020
subroutine getEffects(nobs, maxid, nfix, nvar, theta, Gmatrix, Vhat, Py, y, X,&
     id, fixeff, raneff, verbose)
  use constants
  use global_module
  implicit none
  logical                                                             :: verbose
  integer, intent(in)                                                 :: nobs, nfix, nvar, maxid
  double precision, dimension(:), intent(in)                          :: Gmatrix
  integer, dimension(:), intent(in)                                   :: id
  double precision, dimension(:), intent(in)                          :: theta, Py, y
  double precision, dimension(:,:), intent(in)                        :: Vhat, X
  double precision, dimension(:), intent(out)                         :: fixeff
  type (doublePre_Array), dimension(:), intent(out)                   :: raneff

  double precision, dimension(:), allocatable                         :: temp
  integer                                                             :: i, j
    double precision                                                    :: val1, s1, s2
  type (doublePre_Array), dimension(:), allocatable, target           :: theZPy
  ! allocation
  allocate(theZPy(3))
  allocate(theZPy(1)%level(maxid)) ! slope effect (genetic)
  allocate(theZPy(2)%level(maxid)) ! intercept effect (genetic)
  allocate(theZPy(3)%level(nobs))   ! environment slope effect (diagonal)

  ! initialisation
  do i = 1, 3
     theZPy(i)%level(:) = 0.d0
     raneff(i)%level(:) = 0.d0
  end do

  if (verbose) write(stdout, *) "inside getEffects; initialisation done"
  allocate(temp(maxid))

  ! fixed effects
  call dgemm('n', 'n', nfix, 1, nobs, 1.d0, Vhat, nfix, y, nobs, 0.d0, fixeff, i)

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
