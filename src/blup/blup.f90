subroutine blup(id, X, y, nfix, nobs, maxid, nelement, Gmatrix, nvar, nran,&
     theta, fixEffects, ranEffects, verbose, ipiv, Py, P, V, Vhat, temp, info)

  use constants, only: KINDR, Jarr, STDERR, STDOUT, alloc1D, ONE, ZERO
  use global_module, only : detInv
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
  integer, intent(inout) :: info

  real(KINDR), dimension(nfix), intent(out) :: fixEffects
  type(Jarr), dimension(1:nran), intent(inout) :: ranEffects

  type(Jarr), dimension(1:nvar) :: theZGZ
  real(KINDR), dimension(:), allocatable :: work
  ! nelements = nobs * (nobs + 1) / 2
  integer :: i, j
  real(KINDR) :: val1, val0
  external :: dspmv
  external :: getMatrices, getMatricesUncorrelated, getMatricesCorrelated
  external :: calculateP, calculateV, getEffects
  !! ================ No defintion after this line ================ !!
  i = nobs * nobs
  call alloc1D(work, I, "work", "blup")

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
        call alloc1D(theZGZ(j)%array, nobs, "theZGZ(3)%array", "blup")
     else
        call alloc1D(theZGZ(j)%array, i, "theZGZ(j)%array", "blup")
     end if
  end do

  if (nvar .eq. 3) then
     call getMatricesUncorrelated(verbose, nobs, nfix, maxid, X, Gmatrix, id, &
          theZGZ(1)%array, theZGZ(2)%array, theZGZ(3)%array)
  elseif (nvar .eq. 4) then
     call getMatricesCorrelated(verbose, nobs, nfix, maxid, X, Gmatrix, id, &
          theZGZ(1)%array, theZGZ(2)%array, theZGZ(3)%array, theZGZ(4)%array)
  else
     call getMatrices(verbose, nobs, nfix, maxid, X, Gmatrix, id, theZGZ(1)%array)
  end if
  if (verbose) write(STDOUT, *) " ZGZ created"

  call calculateV(nobs, nvar, theta, theZGZ, V, verbose)
  if (verbose) write(STDOUT, *) " V is calculated"

  call detInv(nobs, V, val1, ipiv, work, verbose, info)
  if (info /= 0) then
     write(STDOUT, 211)
     return
  end if
  if (verbose) write(STDOUT, *) " V is replaced by its inverse"

  call calculateP(nobs, nfix, V, X, P, val0, Vhat, work, temp, verbose, info)
  if (info /= 0) then
     write(STDOUT, 211)
     return
  end if
  if (verbose) write(STDOUT, *) " P is calcuated"

  val0 = ZERO
  val1 = ONE
  i = 1
  call dspmv('u', nobs, val1, P, y, i, val0, Py, i)
  if (verbose) write(STDOUT, *) "  DSPMV finished calculating Py (=P * y)"

  call getEffects(nobs, maxid, nfix, nvar, nran, theta, Gmatrix, Vhat,&
       Py, y, X, id, fixEffects, ranEffects, verbose)
  if (verbose) write(STDOUT, *) " Effects are estimated"

211 format("   warning: error in blup")

  info = 0 ! successful completion
  deallocate(work)
end subroutine blup
