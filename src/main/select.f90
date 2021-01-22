program select
  use constants
  use global_module
  use rng_module
  use evolution_module
  use blup_module
  use reml_module

  implicit none

  logical, parameter :: verbose = .true.
  integer, parameter :: nChr = 1
  integer :: genestart, nanim, istore, nanimNext
  integer :: nloci, nblock, maxloci, maxblock, ifail
  character(len=256) :: startFile, filename1, filename2
  character(len=30) :: formato
  type(chromosome), dimension(:), allocatable, target :: genome1, genome2
  type(chromosome), dimension(:), pointer :: Parentgenome, Offgenome, thisGenome
  integer, dimension(:), allocatable :: seed
  integer, dimension(:), pointer :: indiv

  logical :: random
  real(KINDR), dimension(:,:), allocatable :: frequency
  integer :: ngen, withpos
  real(KINDR) :: chrL, mutationRate

  !ncomp is the number of traits by a QTL (slope, intercept --> 2, otherwise 1)
  integer, parameter :: nSNP = 1000, nQTL = 500, nComp = 2
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array) :: QTLlist
  real(KINDR), dimension(:,:), allocatable :: covMat

  real(KINDR), allocatable, dimension(:,:) :: TBV, locations, CTE, phenotypes
  real(KINDR), allocatable, dimension(:) :: Amat, means
  integer, allocatable, dimension(:) :: ids
  integer :: ivar, iscaled, imiss, addDom
  real(KINDR) :: v1, v2
  type(variance) :: vars
  real(KINDR) :: rand
  integer :: i, j, k, ichr, iun, iun2

  startfile = "inicio.dat"
  call istart(seed, startfile, ifail)
  if (ifail /= 0) then
     write(STDERR, '(a)') "reading/setting seed faild"
     stop 2
  end if
  call random_seed(put = seed)

  ngen = 1
  chrL = 1.d0 !TODO: is this an input?
  mutationRate = 1.e-6 !TODO: is this an input?

  ! First individuals are simulated
  genestart = 3
  nanim = 900 ! TODO: this is an input
  istore = 1
  withpos = 1
  filename1 = "pedigree.ch"
  filename2 = "SNPpositions.t"
  allocate(genome1(nChr))
  allocate(genome2(nChr))
  Parentgenome => genome1
  Offgenome    => genome2

  allocate(vars%A(nComp), vars%E(nComp), vars%PE(nComp), vars%corr(nComp,nComp))
  vars%A = (/ 1.0_KINDR, 3.0_KINDR /)
  vars%E = (/ 9.0_KINDR, 7.0_KINDR /)
  vars%PE = ZERO
  vars%corr(1,:) = (/ ONE, 0.4_KINDR /)
  vars%corr(2,:) = (/ 0.4_KINDR, ONE /)
  ! read the files to initialise the genotype of the first generation
  call initialiseGenotypes(verbose, nchr, nanim, genestart, nloci, nblock, istore, genome1,&
       maxloci, maxblock, ifail, filename1)

  random = .true.
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

  allocate(indiv(nanim), CTE(nComp, 2))
  indiv= (/( i, i = 1, nanim )/)
  call SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, genome1, QTLlist, SNPlist, TBV, verbose)

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

  allocate(locations(1,1))
  locations(1,1) = 1_kindr
  allocate(means(2))
  means = (/ -2.0, 11.0 /)

  call SimulatePhenotype(nAnim, nComp, indiv, TBV, vars, means, locations, ids, phenotypes, "cov", v1)

!  open(1, file = 'phentest')
!  do i = 1, size(phenotypes,1)
!     write(1, *) ids(i), phenotypes(i,:)
!  end do
!  close(1)
  ! genomic evaluation
  ! making Gmatrix
  iscaled = 1 !(0:no, 1:yes)
  ivar  = 1!(0:sample, 1:2pq, 2:2p'q')
  imiss = 0!(0:mean, 1:ignore)
  addDom = 1!(1:additive, 2:dominance)
  call getGmatrix(nanim, nChr, nSNP, indiv, genome1, SNPlist, iscaled, ivar, &
       imiss, addDom, Amat, verbose)

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
end program select
