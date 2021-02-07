subroutine SimulatePhenotype(verbose, nAnim, nComp, indiv, TBV, vars, means, &
     nobs, locations, ids, phen, proc, X)

  use constants
  use rng_module
  implicit none

  logical, intent(in) :: verbose
  integer, intent(in) :: nAnim, nComp
  integer, dimension(nanim), intent(in) :: indiv
  real(KINDR), dimension(nanim, ncomp), intent(in) :: TBV
  type(variances), intent(in) :: vars
  real(KINDR), dimension(ncomp), intent(in) :: means
  integer, intent(in) :: nobs
  real(KINDR), dimension(:,:), intent(in) :: locations
  integer, dimension(nobs), intent(out) :: ids
  real(KINDR), dimension(nobs), intent(out) :: phen
  character(len = *) :: proc
  real(KINDR), dimension(:,:), allocatable, intent(out) :: X

  integer :: i, j, k, len
  real(KINDR), dimension(:,:), allocatable, save :: temp2
  real(KINDR), dimension(:), allocatable, save :: tempr ! to hold mean of zero
  real(KINDR), allocatable, dimension(:,:), save :: E!, A, PE

  external :: covariate

  if (verbose) write(STDOUT, *) " simulating phenotypes"

  if (.not.allocated(temp2)) then
     allocate(temp2(nComp, nComp), tempr(nComp))
     ! covariance structure for E (covariances are zero)
     temp2(1:nComp, 1:nComp) = ZERO
     do i = 1, nComp
        temp2(i,i) = vars%E(i)
        tempr(1:nComp) = ZERO
     end do
  end if

  j = size(locations)

  if (j == 1) then
     i = 1
     len = nAnim
     allocate(X(len, i))
     ids = indiv
  elseif ((j < nAnim) .and. (size(locations, 1) .eq. 1)) then
     i = 2
     len = j * nAnim
     allocate(X(len, i))
     do k = 1, j
        ids(k:len:j) = indiv(1:nAnim)
        X(k:len:j, 1) = locations(1, k) !for mu_slope (locations are repeated)
     end do
  elseif (size(locations, 1) == nAnim) then
     i = 2
     len = j
     allocate(X(len, i))
     j = size(locations, 2)
     do k = 1, j
        ids(k:len:j) = indiv(1:nAnim) ! ids are repeated
        X(k:len:j, 1) = locations(1:nAnim, k) ! different locations
     end do
  else
     write(STDERR, *) "Error: The format of 'locations' is wrong"
     write(STDERR, *) "Use either:"
     write(STDERR, *) " - [size:    1 x 1   ] all animals in the same location"
     write(STDERR, *) " - [size:    1 x nLox] all have phenotypes in nlox places"
     write(STDERR, *) " - [size:nAnim x nLox] phenotypes at different places"
     stop 2
  end if
  if (.not.allocated(E)) allocate(E(len,nComp))
  X(1:len, i) = ONE ! for mu_int
  call gnormal(tempr, temp2, nComp, len, E)
   ! todo: implementation for PE if really it is required
  select case (proc(1:3))
  case ("COV", "COv", "CoV", "Cov", "cOV", "cOv", "coV", "cov")
     call covariate(nComp, len, nAnim, TBV, E, phen, locations, &
          size(locations,1), size(locations,2), means)
  case default
     write(STDERR, *) "error:"
     write(STDERR, *) "case '", proc, "' not implemented"
     stop 2
  end select

  if (verbose) write(STDOUT, *) " end of SimulatePhenotype subroutine"
end subroutine SimulatePhenotype
