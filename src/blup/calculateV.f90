subroutine calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose) 
  use constants
  use global_module
  implicit none
  logical, intent(in)                              :: verbose
  integer, intent(in)                              :: nobs, nvar
  double precision, dimension(:), intent(in)       :: theta
  type (doublePre_array), dimension(:), intent(in) :: theZGZ
  integer, intent(out)                             :: ifail
  double precision, dimension(:), intent(out)      :: V

  integer, parameter                               :: k = 3  ! the index of the diagonal matrix ZsZs
  integer                                          :: i, ipos, irow, isize
  double precision                                 :: Val1

  if (verbose) write(stdout,*) "  In the subroutine CalculateV"
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
  if (verbose) write(stdout,*) "  calculateV returned succesfully"
end subroutine calculateV

