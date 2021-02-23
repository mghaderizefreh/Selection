subroutine GetMutRecArray(verbose, maxVal, chrL, mutationRate, nLoci, &
     chiasmaCumP, chiasmacumP0, totalChiasma, maxchiasma,&
     mutationCumP, mutationCumP0, totalMutation, maxmutations)

  use constants
  use rng_module
  implicit none

  logical, intent(in) :: verbose
  integer, intent(in) :: maxVal, nloci
  real(KINDR), intent(in) :: chrl, mutationRate
  real(KINDR), allocatable, dimension(:), intent(out) :: chiasmaCumP, mutationCumP
  real(KINDR), intent(out) :: chiasmacumP0, mutationCumP0
  integer, allocatable, dimension(:), intent(out) :: totalChiasma, totalMutation
  integer, intent(out) :: maxchiasma, maxmutations

  real(KINDR) :: v1, v2, v3


  allocate(chiasmacumP(maxval), mutationcumP(maxval))
  allocate(totalchiasma(0:maxval), totalmutation(0:maxval))

  totalmutation(0:maxval) = 0
  totalchiasma(0:maxval) = 0

  maxchiasma = 0
  call poissonProb(chrL, maxchiasma, chiasmacumP0)
  v1 = chiasmacumP0
  if (verbose) write(STDOUT, *)' # recom', maxchiasma, chiasmacumP0

  do
     maxchiasma = maxchiasma + 1
     call poissonProb(chrL,maxchiasma,v2)
     v1 = v1 + v2
     if (v1 > ONE) v1 = ONE
     chiasmacumP(maxchiasma) = v1
     if (verbose) write(STDOUT,*) ' # recom', maxchiasma, v2, v1
     if (abs(ONE - v1) <= 0.0000001_KINDR) then
        chiasmacumP(maxchiasma) = ONE
        exit
     end if
  end do

  ! calculating the cumulative distribution for number of mutation events
  maxmutations = 0
  v3 = mutationRate * nloci
  if (verbose) write(STDOUT, *) " average number of mutations", v3
  call poissonProb(v3, maxmutations, mutationcumP0)
  v1 = mutationcumP0
  if (verbose) write(STDOUT, *) " # mutations", maxmutations, mutationcumP0, v1
  do
     maxmutations = maxmutations + 1
     call poissonProb(v3, maxmutations, v2)
     v1 = v1 + v2
     if (v1 > ONE) v1 = ONE
     mutationcumP(maxmutations) = v1
     if (verbose) write(STDOUT, *) " # mutations", maxmutations, v2, v1
     if (abs(ONE - v1) <=  0.0000001_KINDR) then
        mutationcumP(maxmutations) = ONE
        exit
     end if
  end do

end subroutine GetMutRecArray
