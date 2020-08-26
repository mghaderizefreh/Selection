subroutine calculateAImatrix(nobs, nvar, P, AI, f, verbose)
  use constants
  use global_module
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs, nvar
  double precision, dimension(:), intent(in)                          :: P
  double precision, dimension(:), intent(out)                         :: AI
  type (ArrOfArr), dimension(:), intent(in)                           :: f
  double precision, external                                          :: ddot
  integer                                                             :: i, k, j
  double precision, dimension(:), allocatable, save                   :: temp

  if (verbose) write(stdout, *) "  In the subroutine calculateAImatrix"
  if (.not. allocated(temp))  allocate(temp(nobs))

  k = 1
  do i = 1, (nvar + 1)
     call dspmv('u', nobs, 1.d0, P, f(i)%array, 1, 0.d0, temp, 1)

     if (verbose) write(stdout, '(2x,a32, i1, a1)') "DSPMV finished calculating P*f(", i, ")"
     do j = 1, i
        AI(k) = 0.5d0 * ddot(nobs, temp, 1, f(j)%array, 1)
        k = k + 1 
     end do
  end do

  if (verbose) write(stdout, *) "  calculateAImatrix returend successfully"
end subroutine calculateAImatrix

