program select
  use constants
  use global_module
  use rng_module
  use evolution_module
  use blup_module
  use reml_module

  implicit none

  integer, parameter :: nChr = 1
  integer :: genestart, nanim, istore, nanimNext
  integer :: nloci, nblock, maxloci, maxblock, ifail
  character(len=256) :: startFile, filename1, filename2
  character(len=30) :: formato
  type(chromosome), dimension(:), allocatable, target :: genome1, genome2
  type(chromosome), dimension(:), pointer :: Parentgenome, Offgenome, thisGenome
  integer, dimension(:), allocatable :: seed

  double precision, dimension(:,:), allocatable :: frequency
  integer :: ngen, withpos
  double precision :: chrL, mutationRate

  !ncomp is the number of traits by a QTL (slope, intercept --> 2, otherwise 1)
  integer, parameter :: nSNP = 1000, nQTL = 500, nComp = 1 
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array), dimension(nChr) :: QTLlist
  
  real :: rand
  integer :: i, j, k, ichr, iun, iun2

  startfile = "inicio.dat"
  call istart(seed, startfile, ifail)
  if (ifail /= 0) then
     write(stderr, '(a)') "reading/setting seed faild"
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
  filename1 = "dumpedigree.ch"
  filename2 = "SNPpositions.t"
  allocate(genome1(nChr))
  allocate(genome2(nChr))
  Parentgenome => genome1
  Offgenome    => genome2
  
  ! read the files to initialise the genotype of the first generation
  call initialiseGenotypes(nchr, nanim, genestart, nloci, nblock, istore, genome1,&
       maxloci, maxblock, ifail, filename1)

  ! TODO : make list creator a subroutine returing QTL and SNP list
  ! inputs: nChr, nQTL, nComp, randomOrNot
  allocate(SNPList(nChr, nSNP))
  ! otherwise, name of files (QTLlist and SNPlist and nComp)
  open(newUnit = iun, file = "QTLlist.txt")
  open(newUnit = iun2, file = "SNPlist.txt")
  
  write(formato, '(a,i1,a)') "(i3,3x,i6,3x,",nComp,"f15.7)"
  do iChr = 1, nChr
     QTLlist(iChr)%nQTL = nQTL
     QTLlist(iChr)%nComp = 1
     allocate(QTLlist(iChr)%indices(QTLlist(iChr)%nQTL), &
          QTLlist(iChr)%values(QTLlist(iChr)%nQTL,QTLlist(iChr)%nComp))
     do i = 1, QTLlist(iChr)%nQTL
        read(iun, formato) k, QTLlist(iChr)%indices(i), &
             (QTLlist(iChr)%values(i,j), j = 1, ncomp)
     end do
     do i = 1, nSNP
        read(iun2, '(i3,3x,i6)') k, SNPlist(ichr, i)
     end do
  end do

  allocate(frequency(maxloci,2))
  write(stdout, '(a)') " reading positions"
  do ichr=1,nchr
     genome1(ichr)%chrL = chrL
     allocate(genome1(ichr)%positions(genome1(ichr)%nloci))
     write(filename1, '(a,i3.3)') trim(filename2),ichr
     open(newunit = iun, file = filename1, status = 'old')
     read(iun, *) ! skip first line
     do i = 1, genome1(ichr)%nloci
        read(iun, *) k, rand
        if (rand < 0.d0 .or. rand > genome1(ichr)%chrL) then
           write(stderr, *) " error, wrong position", rand
           write(stderr, *) " should be between 0 and ", genome1(ichr)%chrL
           stop 2
        end if
        genome1(ichr)%positions(i) = rand
     end do
  end do

  ! now initialise the offpring genome
  nanimNext = 2000 ! number of individuals for the next generation
  do ichr = 1, nchr
     genome2(ichr)%nLoci     =  genome1(ichr)%nLoci
     genome2(ichr)%nblock    =  genome1(ichr)%nblock
     !this array does not need to be initalised,. only allocated
     allocate(genome2(ichr)%genotypes(nanimNext, 2, genome2(ichr)%nblock))   
     genome2(ichr)%chrL      =  genome1(ichr)%chrL
     genome2(ichr)%positions => genome1(ichr)%positions
  end do

  Parentgenome => genome1
  Offgenome    => genome2
  
  !eliminationType = -1 (keep all snps) 
  write(6, *) 'finished'
  

end program select
