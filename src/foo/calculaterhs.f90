subroutine calculaterhs(nobs, nvar, theZGZ, P, Py, rhs, work, verbose)
  use global_module
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs, nvar
  type (doublePre_array), dimension(:), intent(in)                    :: theZGZ
  double precision, dimension(:), intent(in)                          :: P, Py
  double precision, dimension(:), intent(out)                         :: rhs
  double precision, external                                          :: ddot 
  double precision, dimension(:), intent(inout)                       :: work
  integer                                                             :: i, k

  if (verbose) write(6, *) "  In the subroutine calculateRHS"

  rhs(nvar + 1) = -0.5d0 *(traceA(P, nobs) - ddot(nobs, Py, 1, Py, 1)) 
  k = 3
  do i = 1, nvar
     if (i .ne. k) then
        call dspmv('u', nobs, 1.d0, theZGZ(i)%level, Py, 1, 0.d0, work(1:nobs), 1)
        rhs(i) = -0.5d0 * (traceAxB(P, theZGZ(i)%level, nobs) - ddot(nobs, Py, 1, work, 1))
     else
        work(1 : nobs) = theZGZ(i)%level * Py
        rhs(i) = -0.5d0 * (traceAxBdiag(P, theZGZ(i)%level, nobs) - ddot(nobs, Py, 1, work, 1))
     end if
  end do

  if (verbose) write(6, *) "  calculateRHS returned successfully"
end subroutine calculaterhs

