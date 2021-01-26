program run1
  use constants
  use global_module
  use rng_module
  use evolution_module
  use blup_module
  use reml_module

  implicit none

  logical, parameter :: verbose = .true.
  integer, parameter :: nChr = 1
  integer :: genestart, nanim, istore
  integer :: nloci, nblock, maxloci, maxblock, ifail
  character(len=256) :: startFile, filename1, filename2
!  character(len=60) :: formato
  type(chromosome), dimension(:), allocatable, target :: genome1, genome2
  type(chromosome), dimension(:), pointer :: Parentgenome, Offgenome!, thisGenome
  integer, dimension(:), allocatable :: seed
  integer, dimension(:), pointer :: indiv, male, female
  logical, dimension(:), pointer :: sex !true = male, false = female

!  integer :: iunfix, iunvar, iunran
!  character(len = 30) :: fixEffFile, ranEfffile, varFile
  
  logical :: random
  real(KINDR), dimension(:,:), allocatable :: X !incidence matrix
!  real(KINDR), dimension(:,:), allocatable :: frequency
  integer :: ngen, withpos
  real(KINDR) :: chrL = ONE, mutationRate = 1.e-6_KINDR 


  !ncomp is the number of traits by a QTL (slope, intercept --> 2, otherwise 1)
  integer, parameter :: nSNP = 1000, nQTL = 500, nComp = 2
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array) :: QTLlist

  real(KINDR), allocatable, dimension(:,:) :: TBV, locations, CTE
  real(KINDR), allocatable, dimension(:) :: Amat, means, phenotypes, theta, fixeff
  real(KINDR), allocatable, dimension(:) :: chiasmaCumP
  real(KINDR), pointer, dimension(:) :: mutationCumP
  real(KINDR), pointer :: mutationCumP0
  type(doublePre_Array), dimension(:), allocatable :: raneff
  integer, allocatable, dimension(:) :: ids
  integer, allocatable, dimension(:,:) :: pedigree 
  integer :: ivar, iscaled, imiss, addDom, nobs, maxid, nfix, nvar
  real(KINDR) :: v1, v2, v3, chiasmaCumP0
  type(variance) :: vars
  integer :: i, j, iChr
  integer, parameter :: n_m = 10, n_fpm = 30, n_opf = 3
  integer :: maxchiasma, maxmutations

  startfile = "inicio.dat"
  call istart(seed, startfile, ifail)
  if (ifail /= 0) then
     write(STDERR, '(a)') "reading/setting seed faild"
     stop 2
  end if
  call random_seed(put = seed)

  ! ==============================================
  ! simulating first individuals from genepool
  ! ==============================================
  genestart = 3
  nanim = 900 ! TODO: this is an input
  istore = 1
  withpos = 1
  filename1 = "pedigree.ch"
  filename2 = "SNPpositions.t"
  allocate(genome1(nChr))
  allocate(genome2(nChr))
  ! read the files to initialise the genotype of the first generation
  call initialiseGenotypes(verbose, nchr, nanim, genestart, nloci, nblock, &
       istore, genome1, maxloci, maxblock, ifail, filename1)
  ! copy common attributes to genome2
  do iChr = 1, nChr
     genome2(iChr)%nloci = genome1(iChr)%nloci
     genome2(iChr)%chrL = genome1(iChr)%ChrL
     genome2(iChr)%positions => genome1(iChr)%positions
  end do
  ! ==============================================
  ! calculating the cumulative distribution for number of recombination events
  ! ==============================================
  maxchiasma = 0
  i = 500 ! upto 500 recombinations
  allocate(chiasmacumP(i), mutationcumP(i))
!  allocate(totalchiasma(0:i), totalmutation(0:i))

!  ich=0
!  imu=0
!  imuM=0
!  ichM=0
!  totalmutation(:)=0
!  totalchiasma(:)=0

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
     if (abs(ONE - v1) <= 0.0000001) then
        chiasmacumP(maxchiasma) = ONE
        exit
     end if
  end do
  
  allocate(mutationcumP0)  !becuase it was ocnverted topointer
  ! calculating the cumulative distribution for number of mutation events
  maxmutations = ZERO
  v3 = mutationRate * nloci
  if (verbose) write(STDOUT, *) " average number of mutations", v3
  call poissonProb(v3, maxmutations, mutationcumP0)
  v1 = mutationcumP0
  write(STDOUT, *) " # mutations", maxmutations, mutationcumP0, v1
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
  ! ================================================================
  ! Getting QTL and SNP list
  ! ================================================================
  ! initialising genetic variance and correlations
  allocate(vars%A(nComp), vars%E(nComp), vars%PE(nComp), vars%corr(nComp,nComp))
  vars%A = (/ 1.0_KINDR, 3.0_KINDR /)
  vars%E = (/ 9.0_KINDR, 7.0_KINDR /)
  vars%PE = ZERO
  vars%corr(1,:) = (/ ONE, 0.4_KINDR /)
  vars%corr(2,:) = (/ 0.4_KINDR, ONE /)

  call getQTLandSNP(verbose, nChr, nQTL, nSNP, nComp, .true., genome1, QTLlist, &
       SNPlist, vars%corr)
  ! print*, 'reading snplist'
  !  allocate(SNPlist(nChr, k))
  !  open(1, file = 'SNPlist.txt')
  !  do iChr = 1, nChr
  !     do j = 1, 1000
  !        read(1, *) i, k
  !        SNPlist(iChr, j) = k
  !     end do
  !  end do

  ! ================================================================
  ! Simulating TBV
  ! ================================================================
  ! initialising first generation individuals
  allocate(indiv(nanim), sex(nanim), CTE(nComp, 2))
  indiv= (/( i, i = 1, nanim )/)
  sex(1:nanim) = indiv(1:nanim) < 301
  call SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, genome1, QTLlist, SNPlist, &
       TBV, verbose)

  ! this shall run only for the first generation
  CTE(:, 1) = (/( sum(TBV(1:nAnim, i)) / nAnim , i = 1, nComp )/)
  do i = 1, nComp
     v1 = sum(TBV(1:nAnim, i)) / nAnim
     CTE(i, 1) = sum((TBV(1:nAnim, i) - v1) ** 2) / (nAnim - 1)
     QTLlist%values(1:nChr, 1:nQTL, i) = QTLlist%values(1:nChr, 1:nQTL, i) * &
          sqrt(vars%A(i)/ CTE(i, 1))
     TBV(1:nAnim, i) = TBV(1:nAnim, i) * sqrt(vars%A(i)/ CTE(i, 1))
     CTE(i, 2) = sum(TBV(1:nAnim, i)) / nAnim
     TBV(1:nAnim, i) = TBV(1:nAnim, i) - CTE(i, 2)
  end do

  ! ================================================================
  ! Simulating Phenotypes
  ! ================================================================
  ! initialising first generation individuals
  ! getting locations for phenotyped individuals
  allocate(locations(nAnim, 2))
  !  locations(1,1) = 1.0_kindr ! at x = 1
  call random_number(locations)

  ! setting mean values for phenotypes
  allocate(means(2))
  means = (/ -2.0, 11.0 /)
  call SimulatePhenotype(nAnim, nComp, indiv, TBV, vars, means, locations, ids,&
       phenotypes, "cov", v1, X)

  !  open(1, file = 'phentest')
  !  do i = 1, size(phenotypes,1)
  !     write(1, *) ids(i), phenotypes(i)
  !  end do
  !  close(1)

  ! ================================================================
  ! Genomic Evaluation
  ! ================================================================
  ! making Gmatrix
  iscaled = 1 !(0:no, 1:yes)
  ivar  = 1!(0:sample, 1:2pq, 2:2p'q')
  imiss = 0!(0:mean, 1:ignore)
  addDom = 1!(1:additive, 2:dominance)
  call getGmatrix(nanim, nChr, nSNP, indiv, genome1, SNPlist, iscaled, ivar, &
       imiss, addDom, Amat, verbose)

  nobs = size(phenotypes, 1)
  maxid = maxval(ids)
  nfix = 2
  nvar = 4
  ! setting initial guess for variances
  allocate(theta(nvar + 1))
  theta(1:2) = vars%A(1:2)
  theta(3) = vars%E(1)
  theta(5) = vars%E(2)
  theta(4) = vars%corr(1,2) * sqrt(theta(1) * theta(2))
  ! or use "call blup" with same argument list
  call blup(ids, X, phenotypes, nfix, nobs, maxid, Amat, nvar, theta,&
       fixEff, ranEff, verbose, 6, 13)
  if (verbose) write(STDOUT, *) "Genomic evaluation done"
  !====================================
  ! selection
  !====================================
  allocate(male(10), female(300))
  ! assuming ids are sorted and uniq(ids) = indiv
  call selectByIntercept(nanim, indiv, sex, 10, 30, male, female, ranEff)

  allocate(pedigree(nanim, 3))
  pedigree = makePedigree(10, 30, 3, male, female, nanim)

  open(1, file = 'newped')
  do i = 1, 900
     write(1, '(2(i4,x),i4)') pedigree(i, :)
  end do
  close(1)

  i = n_m + n_fpm * n_m
  j = n_m * n_fpm
  do iChr = 1, nChr
     allocate(genome2(iChr)%genotypes(i,2,size(genome1(iChr)%genotypes,3)))
     genome2(iChr)%genotypes(1:n_m,:,:) = &
          genome1(iChr)%genotypes(male(1:n_m),:,:)
     genome2(iChr)%genotypes((n_m+1):i,:,:) = &
          genome1(iChr)%genotypes(female(1:j),:,:)
  end do
  ParentGenome => genome2
  Offgenome    => genome1
  !  allocate(frequency(maxloci,2))
  !  write(STDOUT, '(a)') " reading positions"
  !  do ichr=1,nchr
  !     genome1(ichr)%chrL = chrL
  !     allocate(genome1(ichr)%positions(genome1(ichr)%nloci))
  !     write(filename1, '(a,i3.3)') trim(filename2),ichr
  !     open(newunit = iun, file = filename1, status = 'old')
  !     read(iun, *) ! skip first line
  !     do i = 1, genome1(ichr)%nloci
  !        read(iun, *) k, rand
  !        if (rand < 0.d0 .or. rand > genome1(ichr)%chrL) then
  !           write(STDERR, *) " error, wrong position", rand
  !           write(STDERR, *) " should be between 0 and ", genome1(ichr)%chrL
  !           stop 2
  !        end if
  !        genome1(ichr)%positions(i) = rand
  !     end do
  !  end do
  !
  !  ! now initialise the offpring genome
  !  nanimNext = 2000 ! number of individuals for the next generation
  !  do ichr = 1, nchr
  !     genome2(ichr)%nLoci     =  genome1(ichr)%nLoci
  !     genome2(ichr)%nblock    =  genome1(ichr)%nblock
  !     !this array does not need to be initalised,. only allocated
  !     allocate(genome2(ichr)%genotypes(nanimNext, 2, genome2(ichr)%nblock))   
  !     genome2(ichr)%chrL      =  genome1(ichr)%chrL
  !     genome2(ichr)%positions => genome1(ichr)%positions
  !  end do
  !
  !  Parentgenome => genome1
  !  Offgenome    => genome2
  !  
  !  !eliminationType = -1 (keep all snps) 
  write(6, *) 'finished'
  !  
  !
end program run1
