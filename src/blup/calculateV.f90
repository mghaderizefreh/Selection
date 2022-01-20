subroutine calculateV(nobs, nvar, theta, theZGZ, V, verbose)
  use constants, only : KINDR, Jarr, STDOUT, ZERO
  implicit none
  logical, intent(in) :: verbose
  integer, intent(in) :: nobs, nvar
  real(KINDR), dimension(1:(nvar+1)), intent(in) :: theta
  type(Jarr), dimension(1:nvar), intent(in) :: theZGZ
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(out) :: V !main output

  integer :: i, ipos, irow, isize
  real(KINDR) :: val1

  associate ( k => 3)! the index of the diagonal matrix ZsZs
    if (verbose) write(STDOUT,*) "  In the subroutine CalculateV"
    isize = (nobs + 1) * nobs/2
    V(1 : isize) = ZERO
    val1 = theta(nvar + 1)
    ipos = 0

    if (nvar .ge. k) then
       do irow = 1, nobs
          ipos    = ipos + irow   ! the position of this diagonal
          V(ipos) = val1 + theta(k) * theZGZ(k)%array(irow)
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
          V(1 : isize) = V(1 : isize) + theZGZ(i)%array(1 : isize) * theta(i)
       end if
    end do
  end associate
  if (verbose) write(STDOUT,*) "  calculateV returned succesfully"
end subroutine calculateV

