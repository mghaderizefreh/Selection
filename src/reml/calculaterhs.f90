subroutine calculaterhs(nobs, nvar, theZGZ, P, Py, rhs, f, verbose)
  use constants
  use global_module
  implicit none
  logical, intent(in)                                 :: verbose
  integer, intent(in)                                 :: nobs, nvar
  type (doublePre_array), dimension(:), intent(in)    :: theZGZ
  real(KINDR), dimension(:), intent(in)          :: P, Py
  real(KINDR), dimension(:), intent(out)         :: rhs
  type (JArrD), dimension(:), intent(inout)           :: f
  real(KINDR), external                          :: ddot 
  integer                                             :: i
  integer, parameter                                  :: k=3!ind of diag mat ZsZs

  external                                            :: dspmv

  if (verbose) write(STDOUT, *) "  In the subroutine calculateRHS"
  
  f(nvar+1)%array(1:nobs) = Py(1:nobs)
  do i = 1, nvar 
     if (i .ne. k) then
        call dspmv('u', nobs, 1.d0, theZGZ(i)%level, Py, 1, 0.d0, f(i)%array, 1)
     else
        f(i)%array(1:nobs) = theZGZ(i)%level(1:nobs) * Py(1:nobs)
     end if
     if (verbose) write(STDOUT, '(2x,a30, i1, a1)') "DSPMV finished calculating f(", i, ")"
  end do

  do i = 1, nvar
     if (i .ne. k) then
        rhs(i) = -0.5d0 * (traceAxB(P, theZGZ(i)%level, nobs) - ddot(nobs, Py, 1, f(i)%array, 1))
     else
        rhs(i) = -0.5d0 * (traceAxBdiag(P, theZGZ(i)%level, nobs) - ddot(nobs, Py, 1, f(i)%array, 1))
     end if
  end do
  rhs(nvar + 1) = -0.5d0 *(traceA(P, nobs) - ddot(nobs, Py, 1, f(nvar+1)%array, 1)) 

  if (verbose) write(STDOUT, *) "  calculateRHS returned successfully"
end subroutine calculaterhs

