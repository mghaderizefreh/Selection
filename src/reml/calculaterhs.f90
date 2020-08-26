subroutine calculaterhs(nobs, nvar, theZGZ, P, Py, rhs, f, verbose)
  use constants
  use global_module
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs, nvar
  type (doublePre_array), dimension(:), intent(in)                    :: theZGZ
  double precision, dimension(:), intent(in)                          :: P, Py
  double precision, dimension(:), intent(out)                         :: rhs
  type (ArrOfArr), dimension(:), intent(inout)                        :: f
  double precision, external                                          :: ddot 
  integer                                                             :: i, k

  if (verbose) write(stdout, *) "  In the subroutine calculateRHS"
  
  f(nvar+1)%array(1:nobs) = Py(1:nobs)
  k = 3 ! index of the diagonal matrix ZsZs
  do i = 1, nvar 
     if (i .ne. k) then
        call dspmv('u', nobs, 1.d0, theZGZ(i)%level, Py, 1, 0.d0, f(i)%array, 1)
     else
        f(i)%array(1:nobs) = theZGZ(i)%level(1:nobs) * Py(1:nobs)
     end if
     if (verbose) write(stdout, '(2x,a30, i1, a1)') "DSPMV finished calculating f(", i, ")"
  end do

  k = 3
  do i = 1, nvar
     if (i .ne. k) then
        rhs(i) = -0.5d0 * (traceAxB(P, theZGZ(i)%level, nobs) - ddot(nobs, Py, 1, f(i)%array, 1))
     else
        rhs(i) = -0.5d0 * (traceAxBdiag(P, theZGZ(i)%level, nobs) - ddot(nobs, Py, 1, f(i)%array, 1))
     end if
  end do
  rhs(nvar + 1) = -0.5d0 *(traceA(P, nobs) - ddot(nobs, Py, 1, f(i)%array, 1)) 

  if (verbose) write(stdout, *) "  calculateRHS returned successfully"
end subroutine calculaterhs

