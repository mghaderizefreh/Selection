subroutine SimulatePhenotype(verbose, nAnim, nComp, nFix, nLox, nran,&
     indiv, TBV, vars, means, nobs, locations, ids, phen, proc, X)
! note that nfix and ncomp are different. nComp refers to number of components
! in simulation: for intercept and slope ncomp = 2. However, nfix depends on the
! type of analysis: if a single trait is desired nfix = 1, if random regression
! is to be conducted nfix = 2
  use constants
  use rng_module
  implicit none

  logical, intent(in) :: verbose
  integer, intent(in) :: nAnim, nComp, nLox, nFix, nran
  integer, dimension(nanim), intent(in) :: indiv
  real(KINDR), dimension(nAnim, nComp), intent(in) :: TBV
  type(variances), intent(in) :: vars
  real(KINDR), dimension(nComp), intent(in) :: means
  integer, intent(in) :: nobs
  real(KINDR), dimension(nAnim,nLox), intent(in) :: locations
  integer, dimension(nObs), intent(out) :: ids
  real(KINDR), dimension(nObs), intent(out) :: phen
  character(len = *) :: proc
  real(KINDR), dimension(nobs,nFix), intent(out) :: X

  integer :: i, k
  real(KINDR), dimension(:,:), allocatable, save :: temp2
  real(KINDR), dimension(:), allocatable, save :: tempr ! to hold mean of zero
  real(KINDR), allocatable, dimension(:,:), save :: E!, A, PE

  external :: covariate

  X(1:nobs, 1:nfix) = ZERO
  if (.not.allocated(temp2)) then
     allocate(temp2(nComp, nComp), tempr(nComp))
     ! covariance structure for E (covariances are zero)
     temp2(1:nComp, 1:nComp) = ZERO
     do i = 1, nComp
        temp2(i,i) = vars%E(i)
     end do
     tempr(1:nComp) = ZERO
  end if
  if (nRan .eq. 3) then
     do k = 1, nlox
        X(k:nobs:nLox, 1) = locations(1:nAnim, k) ! for mu_slo if nfix == 2
     end do
  end if
  do k = 1, nlox
     ids(k:nobs:nLox) = indiv(1:nAnim) ! ids are repeated
  end do
  X(1:nobs, nFix) = ONE ! for mu_int

  if (.not.allocated(E)) allocate(E(nobs,nComp))
  call gnormal(tempr, temp2, nComp, nobs, E)
   ! todo: implementation for PE if really it is required

  select case (proc(1:3))
  case ("COV", "COv", "CoV", "Cov", "cOV", "cOv", "coV", "cov")
     if (verbose) write(STDOUT, *) "Using locations as covariates"
     ! location is used instead of X because X is to be used in analysis
     ! and varies in dimension depending on the type of analysis
     call covariate(nComp, nObs, nAnim, nLox, TBV, E, phen, locations, means)
     if (verbose) write(STDOUT, *) " number of locations used:", nlox
  case default
     write(STDERR, '(a)') "Error:"
     write(STDERR, *) "case '", proc, "' not implemented"
     stop 2
  end select

end subroutine SimulatePhenotype


!!!! ============================================================ !!!!
subroutine allocateInd(nAnim, nlox, nfarm, allocation, farmBounds, &
     farmRange, farmInd, locations)
  use constants
  implicit none
  integer, intent(in) :: nAnim, nLox, nFarm
  integer, intent(in) :: allocation
  real(KINDR), dimension(1:nfarm, 2), intent(in) :: farmBounds
  real(KINDR) :: farmRange
  integer, dimension(1:(nAnim*nLox)), intent(out) :: farmInd
  real(KINDR), dimension(nAnim, nLox), intent(out) :: locations

  integer :: nobs 
  integer :: i, j, k
  real(KINDR) :: val
  real(KINDR), dimension(nlox) :: temp
  nobs = nlox * nAnim
  select case (allocation)
  case(1) ! random
     ! each individual is assigned to (only) one random farm
     do i = 1, nAnim
        call random_number(val)
        do while(val.eq.ONE)
           call random_number(val)
        end do
        j = int(val * nfarm) + 1
        k = nlox * (i - 1) + 1
        farmInd(k:(nlox*i)) = j
        call random_number(temp)
        ! shifting and scaling to match interval
        temp(1:nlox) = temp(1:nlox) * farmRange + farmBounds(j, 1)
        locations(i, 1:nlox) = temp(1:nlox)
        ! random_number is in [0,1] so no need to shift by 0 and divide by 1
     end do
  case(2:)
  end select
  
end subroutine allocateInd

!!!! ============================================================ !!!!
subroutine defineFarms(interval, nfarm, diameter, farms)
  use constants
  use quickSort
  implicit none
  integer, intent(in) :: nfarm
  real(KINDR), dimension(2), intent(in) :: interval
  real(KINDR), intent(in) :: diameter
  real(KINDR), dimension(nfarm, 2), intent(out) :: farms

  real(KINDR), dimension(nfarm) :: centres
  integer, dimension(nfarm) :: ind
  
  call random_number(centres)
  call sortrx(nfarm, centres, ind)
  centres(1:nfarm) = centres(ind(1:nfarm))
  centres(1:nfarm) = (centres(1:nfarm) - centres(1)) / &
       (centres(nfarm) - centres(1)) * &
       (interval(2) - interval(1) - diameter) +&
       interval(1) + diameter / 2
  farms(1:nfarm, 1) = centres(1:nfarm) - diameter / 2
  farms(1:nfarm, 2) = centres(1:nfarm) + diameter / 2
  ! arithmetic floating point messes up with first and last boundary
  ! although experiment showed this is not required most of the times
  farms(1      , 1) = interval(1)
  farms(nfarm  , 2) = interval(2)
    
end subroutine defineFarms
