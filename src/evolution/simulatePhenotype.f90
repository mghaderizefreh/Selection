subroutine SimulatePhenotype(verbose, nAnim, nComp, nFix, nLox, nran,&
     indiv, TBV, vars, means, nobs, locations, ids, phen, proc, X)
! note that nfix and ncomp are different. nComp refers to number of components
! in simulation: for intercept and slope ncomp = 2. However, nfix depends on the
! type of analysis: if a single trait is desired nfix = 1, if random regression
! is to be conducted nfix = 2
  use constants, only: KINDR, variances, alloc2D, ZERO, ONE, STDERR, STDOUT
  use rng_module, only: gnormal
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
  real(KINDR), dimension(1:ncomp,1:ncomp) :: temp2
  real(KINDR), dimension(1:ncomp) :: tempr ! to hold mean of zero
  real(KINDR), allocatable, dimension(:,:) :: E!, A, PE

  external :: covariate

  X(1:nobs, 1:nfix) = ZERO
  
  ! covariance structure for E (covariances are zero)
  temp2(1:nComp, 1:nComp) = ZERO
  do i = 1, nComp
     temp2(i,i) = vars%E(i)
  end do
  tempr(1:nComp) = ZERO
  
  if (nRan .eq. 3) then
     do k = 1, nlox
        X(k:nobs:nLox, 1) = locations(1:nAnim, k) ! for mu_slo if nfix == 2
     end do
  end if
  do k = 1, nlox
     ids(k:nobs:nLox) = indiv(1:nAnim) ! ids are repeated
  end do
  X(1:nobs, nFix) = ONE ! for mu_int

  call alloc2D(E, nobs,nComp, "e", "simulatePhenotype")
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

  deallocate(E)

end subroutine SimulatePhenotype


!!!! ============================================================ !!!!
subroutine allocateInd(nAnim, nlox, nobs, nfarm, allocation, farmBounds,&
     farmRange, farmInd, locations, pedigree, nm, male)
  use constants, only: KINDR, ONE, alloc1I, STDERR
  use rng_module, only: choice
  implicit none
  integer, intent(in) :: nAnim, nLox, nobs, nFarm
  integer, intent(in) :: allocation
  real(KINDR), dimension(1:nfarm, 2), intent(in) :: farmBounds
  real(KINDR) :: farmRange
  integer, dimension(1:nobs), intent(out) :: farmInd
  real(KINDR), dimension(nAnim, nLox), intent(out) :: locations
  integer, dimension(1:nAnim, 3), optional :: pedigree
  integer, optional, intent(in) :: nm
  integer, optional, intent(in) :: male(:)

  integer, dimension(:), allocatable :: temp1, temp2
  integer :: spf! each farm has spf sires
  integer :: isire
  integer :: i, j, k
  real(KINDR) :: val
  real(KINDR), dimension(nlox) :: temp

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
  case(2) ! some sort of clustering:
     !one farm/sire & spf sires/farm & farms and sires random
     ! i.e., nf farms (as indices 1...nf) are randomised and each farm in this
     ! index array (temp2) will contain spf sires. The mating (pedigree) is 
     ! already decided.
     ! requires pedigree
     if (.not.present(pedigree)) then
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) " pedigree is required for allocation = 2"
        stop 2
     elseif (.not.present(nm)) then
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) " number of sires is required for allocation = 2"
        stop 2
     elseif (.not.present(male)) then
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) " male array is required for allocation = 2"
        stop 2
     elseif (nm < nfarm) then
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) " number of farms is more than number of sires"
        stop 2
     end if
     spf = int(nm / nfarm)
     if (spf * nfarm .ne. nm) then
        write(STDERR,'(a)') "Error:"
        write(STDERR, *) "This is something that should not happen"
        write(STDERR, *) "The number of sires is not divisible to nfarms"
        stop 2
     end if

     call alloc1I(temp1, nm, "temp1", "allocateInd")
     call alloc1I(temp2, nm, "temp2", "allocateInd")
     do i = 1, spf
        j = (i-1)*nfarm + 1
        k = i * nfarm
        temp1(j:k) = (/(isire, isire = 1, nfarm)/)
     end do
     call choice(temp1, nm, nm, nm, temp2, nm) ! temp2 contains farms
     ! could do smarter way, but would require 2 more arrays
     !
     ! the method below requires the pedigree to be sorted, as this is the
     ! case for my simulations, I don't need to make this general.
     ! for the first individual, take the sire
     i = 1
     isire = pedigree(i,2)
     ! find its index in male array
     do j = 1, nm
        if (isire == male(j)) exit
     end do
     ! assign the farm to that individual
     k = temp2(j)
     farmind(((i-1)*nlox+1):(i*nlox)) = k
     call random_number(temp)
     locations(i, 1:nlox) = temp(1:nlox) * farmRange + farmBounds(k, 1)
     do i = 2, nAnim
        !for all other individual, get the sire
        isire = pedigree(i,2)
        ! if its another offspring of the same sire, just assign the same
        ! again, assuming pedigree is sorted, this check is easy
        if (isire == pedigree(i-1, 2)) then
           k = temp2(j)
           farmind(((i-1)*nlox + 1):(i*nlox)) = k
           call random_number(temp)
           locations(i, 1:nlox) = temp(1:nlox) * farmRange + farmBounds(k, 1)
           cycle ! and go to the next individual
        end if
        ! otherwise, find the index and do the rest
        do j = 1, nm
           if (isire == male(j)) exit
        end do
        k = temp2(j)
        farmind(((i-1)*nlox + 1):(i*nlox)) = k
        call random_number(temp)
        locations(i, 1:nlox) = temp(1:nlox) * farmRange + farmBounds(k, 1)
     end do
     deallocate(temp1)
     deallocate(temp2)
  case(3:) ! later
  end select
  
end subroutine allocateInd

!!!! ============================================================ !!!!
subroutine defineFarms(interval, nfarm, diameter, farms)
  use constants, only: KINDR
  use quickSort, only: sortrx
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
