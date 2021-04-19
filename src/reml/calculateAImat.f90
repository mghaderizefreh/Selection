subroutine calculateAImatrix(nobs, nvar, P, AI, f, verbose)
  use constants
  use global_module
  implicit none
  logical, intent(in) :: verbose
  integer, intent(in) :: nobs, nvar
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(in) :: P
  real(KINDR), dimension(1:((nvar+1)*(nvar+2)/2)), intent(out) :: AI
  type (JArr), dimension(1:(nvar+1)), intent(in) :: f
  real(KINDR), external :: ddot
  integer :: i, k, j
  real(KINDR), dimension(:), allocatable :: temp

  external :: dspmv

  if (verbose) write(STDOUT, *) "  In the subroutine calculateAImatrix"
  allocate(temp(nobs))
  k = 1
  do i = 1, (nvar + 1)
     call dspmv('u', nobs, 1.d0, P, f(i)%array, 1, 0.d0, temp, 1)

     if (verbose) write(STDOUT, '(2x,a32, i1, a1)') "DSPMV finished calculating P*f(", i, ")"
     do j = 1, i
        AI(k) = 0.5d0 * ddot(nobs, temp, 1, f(j)%array, 1)
        k = k + 1 
     end do
  end do
  deallocate(temp)
  if (verbose) write(STDOUT, *) "  calculateAImatrix returend successfully"
end subroutine calculateAImatrix

