subroutine calculaterhs(nobs, nvar, theZGZ, P, Py, rhs, f, verbose)
  use constants, only: KINDR, JArr, STDOUT, ZERO, ONE, HALF
  use global_module, only: traceA, traceAxB, traceAxBdiag
  implicit none
  logical, intent(in) :: verbose
  integer, intent(in) :: nobs, nvar
  type (Jarr), dimension(1:nvar), intent(in) :: theZGZ
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(in) :: P
  real(KINDR), dimension(1:nobs), intent(in) :: Py
  real(KINDR), dimension(1:(nvar+1)), intent(out) :: rhs
  type (JArr), dimension(1:(nvar+1)), intent(inout) :: f
  real(KINDR), external :: ddot 
  integer :: i
  external :: dspmv
  associate(k => 3)!ind of diag mat ZsZs

  if (verbose) write(STDOUT, *) "  In the subroutine calculateRHS"
  
  f(nvar+1)%array(1:nobs) = Py(1:nobs)
  do i = 1, nvar 
     if (i .ne. k) then
        call dspmv('u', nobs, ONE, theZGZ(i)%array, Py, 1, ZERO, f(i)%array, 1)
     else
        f(i)%array(1:nobs) = theZGZ(i)%array(1:nobs) * Py(1:nobs)
     end if
     if (verbose) write(STDOUT, '(2x,a30, i1, a1)') "DSPMV finished calculating f(", i, ")"
  end do

  do i = 1, nvar
     if (i .ne. k) then
        rhs(i) = -HALF * (traceAxB(P, theZGZ(i)%array, nobs) - ddot(nobs, Py, 1, f(i)%array, 1))
     else
        rhs(i) = -HALF * (traceAxBdiag(P, theZGZ(i)%array, nobs) - ddot(nobs, Py, 1, f(i)%array, 1))
     end if
  end do
  end associate
  rhs(nvar + 1) = -HALF *(traceA(P, nobs) - ddot(nobs, Py, 1, f(nvar+1)%array, 1)) 

  if (verbose) write(STDOUT, *) "  calculateRHS returned successfully"
end subroutine calculaterhs

