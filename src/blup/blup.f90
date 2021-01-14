subroutine blup(id, X, y, nfix, nobs, maxid, Gmatrix, nvar, theta, &
     fixEffects, ranEffects, verbose, EmIterations, maxIters)

  use constants
  use global_module
  implicit none
  !! ================ variable definitions  ================ !!
  logical, intent(in)                            :: verbose
  integer, intent(in)                            :: maxid, nvar, nobs, nfix
  integer, dimension(:), intent(in)              :: id ! real id of animals
  double precision, dimension(:), intent(in)     :: y ! phenotypes
  double precision, dimension(:,:), intent(in)   :: X ! incid. matrix
  double precision, dimension(:), intent(in)     :: Gmatrix
  double precision, dimension(:),intent(inout)   :: theta
  integer, intent(in), optional                  :: EmIterations, maxIters

  double precision, dimension(:), intent(out)    :: fixEffects
  type(doublePre_Array),dimension(:),intent(out) :: ranEffects

  type(doublePre_Array),dimension(:),allocatable :: theZGZ
  double precision, dimension(:), allocatable    :: Py, P, V, work
  integer                                        :: ifail, i, j
  double precision                               :: val1, val2
  double precision, dimension(:,:), allocatable  :: Vhat
  integer, dimension(:), allocatable             :: ipiv
  double precision, external                     :: dnrm2, ddot, dasum
  !! ================ No defintion after this line ================ !!
  allocate(Py(nobs), Vhat(nfix, nobs))
  I = nobs * (nobs + 1) / 2
  allocate(P(I),V(I))
  I = nobs * nobs
  allocate(work(I),ipiv(nobs))

  if (present(EmIterations)) I = EmIterations
  if (present(maxIters)) I = maxIters

  if (nvar == 3) then
     if (any ( theta < 0 )) then
        write(stderr, *) "n_var and initial guess not consistent"
        stop 2
     end if
     write(stdout, '(2x,a22)') "no correlation assumed"
  elseif (nvar > 3) then
     write(stdout, '(2x,a30)') "correlation taken into account"
  end if

  allocate(theZGZ(nvar))
  i = nobs * (nobs + 1) / 2
  do j = 1, nvar
     if (j .eq. 3) then
        allocate(theZGZ(j)%level(nobs))
     else
        allocate(theZGZ(j)%level(i))
     end if
  end do

  if (nvar .eq. 3) then
     call getMatricesUncorrelated(verbose, nobs, X, Gmatrix, id, &
          theZGZ(1)%level, theZGZ(2)%level, theZGZ(3)%level)
  elseif (nvar .eq. 4) then
     call getMatricesCorrelated(verbose, nobs, X, Gmatrix, id, &
          theZGZ(1)%level, theZGZ(2)%level, theZGZ(3)%level, &
          theZGZ(4)%level)
  else
     call getMatrices(verbose, nobs, X, Gmatrix, id, theZGZ(1)%level)
  end if

  call calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose)
  if (verbose) write(stdout, *) " V is calculated"

  call detInv(nobs, V, val1, ipiv, work, verbose)
  if (verbose) write(stdout, *) " V is replaced by its inverse"

  call calculateP(nobs, nfix, V, X, P, val2, Vhat, verbose)
  if (verbose) write(stdout, *) " P is calcuated"

  call dspmv('u', nobs, 1.d0, P, y, 1, 0.d0, Py, 1)
  if (verbose) write(stdout, *) "  DSPMV finished calculating Py (=P * y)"

  do i = 1, size(ranEffects)
     ranEffects(i)%level(:) = 0.d0
  end do

  call getEffects(nobs, maxid, nfix, nvar, theta, Gmatrix, Vhat, Py, y, X,&
       id, fixEffects, ranEffects, verbose)

end subroutine blup
