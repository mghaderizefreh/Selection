subroutine calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose) 
  use constants
  use global_module
  implicit none
  logical, intent(in) :: verbose
  integer, intent(in) :: nobs, nvar
  real(KINDR), dimension(1:(nvar+1)), intent(in) :: theta
  type (doublePre_array), dimension(1:nvar), intent(in) :: theZGZ
  integer, intent(out) :: ifail
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(out) :: V

  integer, parameter :: k = 3  ! the index of the diagonal matrix ZsZs
  integer :: i, ipos, irow, isize
  real(KINDR) :: val1

  if (verbose) write(STDOUT,*) "  In the subroutine CalculateV"
  ifail = 1
  isize = (nobs + 1) * nobs/2
  v(1 : isize) = 0.d0
  val1 = theta(nvar + 1)
  ipos = 0

  if (nvar .ge. k) then
     do irow = 1, nobs
        ipos    = ipos + irow   ! the position of this diagonal
        V(ipos) = val1 + theta(k) * theZGZ(k)%level(irow)
     end do
  else
     ipos = 0
     do irow = 1, nobs
        ipos = ipos + irow
        V(ipos) = val1
     end do
  end if

  do i = 1, nvar
     if (i .ne. k) then
        v(1 : isize) = v(1 : isize) + theZGZ(i)%level(1 : isize) * theta(i)
     end if
  end do

  ifail=0
  if (verbose) write(STDOUT,*) "  calculateV returned succesfully"
end subroutine calculateV

