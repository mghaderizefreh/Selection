program select
  use constants
  use global_module
  use rng_module
  use evolution_module
  use blup_module
  use reml_module

  implicit none

  logical, parameter :: verbose = .true.
  integer, parameter :: nChr = 26
  integer :: genestart, nanim, istore, nanimNext
  integer :: nloci, nblock, maxloci, maxblock, ifail
  character(len=256) :: startFile, filename1, filename2
  character(len=30) :: formato
  type(chromosome), dimension(:), allocatable, target :: genome1, genome2
  type(chromosome), dimension(:), pointer :: Parentgenome, Offgenome, thisGenome
  integer, dimension(:), allocatable :: seed
  integer, dimension(:), pointer :: indiv

  double precision, dimension(:,:), allocatable :: frequency
  integer :: ngen, withpos
  double precision :: chrL, mutationRate

  !ncomp is the number of traits by a QTL (slope, intercept --> 2, otherwise 1)
  integer, parameter :: nSNP = 1000, nQTL = 500, nComp = 2
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array) :: QTLlist
  real, dimension(:,:), allocatable :: covMat

  double precision, allocatable, dimension(:,:) :: TBV, locations
  double precision, allocatable, dimension(:) :: phenotypes, Amat
  integer, allocatable, dimension(:) :: ids
  integer :: ivar, iscaled, imiss, addDom
  double precision :: v1, v2
  real :: rand
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

  ! read the files to initialise the genotype of the first generation
  call initialiseGenotypes(nchr, nanim, genestart, nloci, nblock, istore, genome1,&
       maxloci, maxblock, ifail, filename1)

  allocate(covMat(ncomp,ncomp))
  covMat(1,:) = (/ 1.0, 0.3/)
  covMat(2,:) = (/ 0.3, 1.0/)
  j = 500
  k = 1000
  !  call getQTLandSNP(nChr, j, k, ncomp, .true., genome1, QTLlist, &
  !       SNPlist, covMat)
  print*, 'reading snplist'
  allocate(SNPlist(nChr, k))
  open(1, file = 'SNPlist.txt')
  do iChr = 1, nChr
     do j = 1, 1000
        read(1, *) i, k
        SNPlist(iChr, j) = k
     end do
  end do
  

  allocate(indiv(nanim))
  indiv= (/( i, i = 1, nanim )/)
!  call SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, genome1, QTLlist, SNPlist, TBV)

!  allocate(locations(1,1))
!  call SimulatePhenotype(nAnim, nComp, indiv, TBV, locations, ids, phenotypes)

  ! genomic evaluation
  ! making Gmatrixo
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
