subroutine calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose) 
  use global_module
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs, nvar
  double precision, dimension(:), intent(in)                          :: theta
  type (doublePre_array), dimension(:), intent(in)                    :: theZGZ
  integer, intent(out)                                                :: ifail
  double precision, dimension(:), intent(out)                         :: V

  integer                                                             :: i, k, ipos, irow, isize
  double precision                                                    :: Val1

  if (verbose) write(6,*) "  In the subroutine CalculateV"
  ifail = 1
  isize = (nobs + 1) * nobs/2
  v(1 : isize) = 0.d0
  val1 = theta(nvar + 1)
  ipos = 0
  k = 3 ! the index of the diagonal matrix ZsZs
  do irow = 1, nobs
     ipos    = ipos + irow   ! the position of this diagonal
     V(ipos) = val1 + theta(k) * theZGZ(k)%level(irow)
  end do

  do i = 1, nvar
     if (i .ne. k) then
        v(1 : isize) = v(1 : isize) + theZGZ(i)%level(1 : isize) * theta(i)
     end if
  end do

  ifail=0
  if (verbose) write(6,*) "  calculateV returned succesfully"
end subroutine calculateV

