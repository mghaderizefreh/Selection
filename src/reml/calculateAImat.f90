subroutine calculateAImatrix(nobs, nvar, theZGZ, P, Py, AI, verbose)
  use constants
  use global_module
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs, nvar
  type (doublePre_array), dimension(:), intent(in)                    :: theZGZ
  double precision, dimension(:), intent(in)                          :: P, Py
  double precision, dimension(:), intent(out)                         :: AI
  double precision, external                                          :: ddot
  integer                                                             :: i, k, j
  type ArrOfArr
     double precision, dimension(:), allocatable                      :: array
  end type ArrOfArr
  type (ArrOfArr), dimension(:), allocatable, save                    :: f
  double precision, dimension(:), allocatable, save                   :: temp

  if (verbose) write(stdout, *) "  In the subroutine calculateAImatrix"
  if (.not. allocated(f)) allocate(f(nvar+1))
  if (.not. allocated(f(nvar+1)%array)) allocate(f(nvar+1)%array(nobs))
  if (.not. allocated(temp))  allocate(temp(nobs))

  f(nvar+1)%array(1:nobs) = Py(1:nobs)
  k = 3 ! index of the diagonal matrix ZsZs
  do i = 1, nvar 
     if (.not. allocated(f(i)%array)) allocate(f(i)%array(nobs))
     if (i .ne. k) then
        call dspmv('u', nobs, 1.d0, theZGZ(i)%level, Py, 1, 0.d0, f(i)%array, 1)
     else
        f(i)%array(1:nobs) = theZGZ(i)%level(1:nobs) * Py(1:nobs)
     end if
     if (verbose) write(stdout, '(2x,a30, i1, a1)') "DSPMV finished calculating f(", i, ")"
  end do

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

