subroutine blup(id, X, y, nfix, nobs, maxid, nelement, Gmatrix, nvar, nran,&
     theta, fixEffects, ranEffects, verbose, ipiv, Py, P, V, Vhat, temp)

  use constants
  use global_module
  implicit none
  !! ================ variable definitions  ================ !!
  logical, intent(in) :: verbose
  integer, intent(in) :: maxid, nvar, nobs, nfix, nran, nelement
  integer, dimension(1:nobs), intent(in) :: id ! real(KINDR) id of animals
  real(KINDR), dimension(1:nobs), intent(in) :: y ! phenotypes
  real(KINDR), dimension(1:nobs,1:nfix), intent(in) :: X ! incid. matrix
  real(KINDR), dimension(1:(maxid*(maxid+1)/2)), intent(in) :: Gmatrix
  real(KINDR), dimension(1:(nvar+1)),intent(inout) :: theta
  integer, dimension(1:nobs), intent(inout) :: ipiv
  real(KINDR), dimension(1:nobs), intent(inout) :: Py
  real(KINDR), dimension(1:nelement), intent(inout) :: P
  real(KINDR), dimension(1:nelement), intent(inout) :: V
  real(KINDR), dimension(1:nfix,1:nobs), intent(inout) :: Vhat
  real(KINDR), dimension(1:nobs,1:nfix), intent(inout) :: temp

  real(KINDR), dimension(nfix), intent(out) :: fixEffects
  type(doublePre_Array), dimension(1:nran), intent(inout) :: ranEffects

  type(doublePre_Array), dimension(1:nvar) :: theZGZ
  real(KINDR), dimension(:), allocatable :: work
  ! nelements = nobs * (nobs + 1) / 2
  integer :: ifail, i, j
  real(KINDR) :: val1, val2
  external :: dspmv
  !! ================ No defintion after this line ================ !!
  i = nobs * nobs
  allocate(work(I))

  if (nvar == 3) then
     if (any ( theta < 0 )) then
        write(STDERR, *) "n_var and initial guess not consistent"
        stop 2
     end if
     if (verbose) write(STDOUT, '(2x,a22)') "no correlation assumed"
  elseif (nvar > 3) then
     if (verbose) write(STDOUT, '(2x,a30)') "correlation taken into account"
  end if

  i = nobs * (nobs + 1) / 2
  do j = 1, nvar
     if (j .eq. 3) then
        allocate(theZGZ(j)%level(nobs))
     else
        allocate(theZGZ(j)%level(i))
     end if
  end do

  if (nvar .eq. 3) then
     call getMatricesUncorrelated(verbose, nobs, nfix, maxid, X, Gmatrix, id, &
          theZGZ(1)%level, theZGZ(2)%level, theZGZ(3)%level)
  elseif (nvar .eq. 4) then
     call getMatricesCorrelated(verbose, nobs, nfix, maxid, X, Gmatrix, id, &
          theZGZ(1)%level, theZGZ(2)%level, theZGZ(3)%level, &
          theZGZ(4)%level)
  else
     call getMatrices(verbose, nobs, nfix, maxid, X, Gmatrix, id, theZGZ(1)%level)
  end if
  if (verbose) write(STDOUT, *) " ZGZ created"

  call calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose)
  if (verbose) write(STDOUT, *) " V is calculated"

  call detInv(nobs, V, val1, ipiv, work, verbose)
  if (verbose) write(STDOUT, *) " V is replaced by its inverse"

  ! val = det(x'*vinv*x)
  call calculateP(nobs, nfix, V, X, P, val2, Vhat, work, temp, verbose)
  if (verbose) write(STDOUT, *) " P is calcuated"

  call dspmv('u', nobs, 1.d0, P, y, 1, 0.d0, Py, 1)
  if (verbose) write(STDOUT, *) "  DSPMV finished calculating Py (=P * y)"

  call getEffects(nobs, maxid, nfix, nvar, nran, theta, Gmatrix, Vhat,&
       Py, y, X, id, fixEffects, ranEffects, verbose)
  if (verbose) write(STDOUT, *) " Effects are estimated"

  deallocate(work)
end subroutine blup
