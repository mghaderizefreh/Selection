!     Last change:  R     6 Oct 2014   11:27 am
program generations
!  use sampling
!  use quickSort
  use rng_module
  implicit none

  integer                                  :: istore  ! how SNP are stored
  integer                                  :: withpos

  character(len = 11) :: seedfile
  integer, dimension(:), allocatable :: nseed
  character(len = 256) :: pedfile, filename, filename1, filename2

  integer :: nchr, ichr
  type(chromosome), dimension(:), allocatable, target  :: genome1, genome2
  type(chromosome), dimension(:),              pointer :: Parentgenome, Offgenome, thisGenome

  integer, dimension(:,:,:), pointer :: parentGen
  integer, dimension(:,:,:), pointer :: offGen      ! if the genotyps of offspring are store in a different array
  integer, dimension(:,:,:), pointer :: this        ! a working pointer to sawp the arrays
  real, dimension(:), allocatable, target :: positions   ! if the loci are not equidistance, then the position (it does not matter if Morgan or cM).


  integer :: nloci,nblock    !nblock and nloci are the same if storage is as integer
  integer :: nremained, newnblock, originalnblock, originalnloci, maxblock, maxloci 

  integer, dimension(:,:), allocatable :: pedigree, finalPedigree
  integer, dimension(:), allocatable :: sires,dams
  integer :: npairs
  integer :: nsires,ndams, nanim, ngen, noff, nfinaloff

  double precision :: chrL, mutationRate

  integer :: maxchiasma, maxmutations
  integer :: nchiasma
  double precision :: chiasmacumP0                           ! probability of no recombination
  double precision, dimension(:), allocatable :: chiasmacumP ! cumulative probability of i recombination

  double precision, pointer :: mutationcumP0              ! probability of no mutation
  double precision, dimension(:), pointer :: mutationcumP ! cumulative probability of i mutations

  double precision, dimension(:), allocatable, target :: AllChrmutationcumP0 ! probability of no mutation
  double precision, dimension(:,:), allocatable, target :: AllChrmutationcumP ! cumulative probability of i mutations

  double precision, dimension(:,:), allocatable :: frequency  ! frequency of both alleles
  integer, dimension(:), allocatable :: eliminated ! flag if SNP will be eliminated
  integer, dimension(:), allocatable :: oldnumber 

  double precision, dimension(:), allocatable :: LD  ! LD between consecutive loci
  double precision :: averLD

  double precision, dimension(:), allocatable :: MAF  ! frequency of both alleles
  integer, dimension(:), allocatable :: MAFsorted  !
  integer :: eliminationType, nkept, genstart
  double precision :: minMAF

  double precision :: val1, val2, val3
  integer :: i, j, k, igen, mutotal, chtotal, a1, a2, iloci, iun, iline, header
  integer :: parent, id, igam, iblock
  integer :: iparent

  integer :: icheck, ifail
  real :: rand

  character(len=1) :: theStatus, expectStatus 

  write(seedfile,'(a)') "inicio.dat"
  call istart(seed, seedfile, i)
  if (i /= 0) then
     write(STDERR, '(a)') "reading/setting seed faild"
     stop 2
  end if
  call random_seed(put = seed)

  !=========================================================
  ! passing parameters
  !=========================================================

  call askInteger(npairs, ' number of mating pairs (n is twice the mating pairs)')
  nsires = npairs
  ndams = npairs
  nanim = npairs*2

  ngen = 1 ! number of generations

  !TODO : this is an input possibly
  nchr = 1 ! it should be 26 or input

  ! chromosome length is 1
  chrL = 1.d0 ! TODO: does this need to be changeable

  ! mutation rate: 1e-6
  ! TODO: does this need to be changeable
  mutationRate = 1.e-6


  genstart = 3 ! meaning the start is from previous run
  !todo : remove this variable


  ! todo: these need to be in form of arrays/objects not files and input to subroutine
  write(*,'(a)')' prefix for genotype file name (genotype for each chr will be in file named <prefix+XXX>)'
  read(*,*) filename1
  write ( *, '(a)' ) ' prefix for position file name (position for each chr will be in file named <prefix+XXX>)'
  read ( *, * ) filename2
  istore=1    !only this option for this one
  withpos=1  ! only this option for this one

  ALLOCATE ( genome1( nChr ) )
  ALLOCATE ( genome2( nChr ) )
  Parentgenome => genome1
  Offgenome    => genome2

  call initialiseGenotypes(nchr,nanim,genstart,nloci, nblock, istore,genome1,nseed,maxloci, maxblock, ifail, filename1)
  write(*,*)'maxloci maxblock',maxloci, maxblock

  ALLOCATE(frequency(maxloci,2))
  if(genstart ==1 .or. genstart==2) then
     ALLOCATE(positions(maxloci))
     frequency(:,:)=0.d0
     val1 = chrL/DBLE(nloci-1)
     positions(1)=0.d0
     do i=2,nloci
        positions(i)=DBLE(i-1)*val1
     end do
     do ichr = 1, nchr
        genome1(ichr)%chrL = chrL
        genome1( ichr )%positions => positions  !all chrmosomes have the same (same equidistant number SNP at the same position)
     enddo
  elseif(genstart==3) then
     write(*,*)' reaing positions'
     do ichr=1,nchr
        genome1(ichr)%chrL = chrL
        allocate( genome1(ichr)%positions( genome1(ichr)%nloci ) )
        WRITE(filename,'(a,i3.3)')trim(filename2),ichr
        OPEN(newunit=iun,FILE=filename,STATUS='old')
        read(iun,*) 
        do i=1,genome1(ichr)%nloci
           read(iun,*)k, rand
           if(rand < 0.d0 .or. rand > genome1(ichr)%chrL )then
              write(*,*)' error, wrong position', rand
              write(*,*)' should be between 0 and ', genome1(ichr)%chrL
              stop
           endif
           genome1(ichr)%positions(i)=rand
        enddo
     enddo
  endif


  ! now initialise the offpring genome
  do ichr=1,nchr
     genome2(ichr)%nLoci     =  genome1(ichr)%nLoci
     genome2(ichr)%nblock    =  genome1(ichr)%nblock
     ALLOCATE ( genome2( ichr )%genotypes( nanim, 2, genome2(ichr)%nblock ) )   !this array does not need to be initalised,. only allocated
     genome2(ichr)%chrL      =  genome1(ichr)%chrL
     genome2(ichr)%positions => genome1(ichr)%positions
  enddo


  Parentgenome => genome1
  Offgenome    => genome2





  icheck=100
  WRITE(*,*)' icheck '
  READ( *,*)icheck
  IF(icheck <=0) icheck =INT(ngen/10)
  IF(icheck <=0) icheck=1
  WRITE(*,*)icheck

  !-------------------------------
  ! output for last generation 
  !-------------------------------
  WRITE(*,*)
  WRITE(*,*)' SNP Output for last generation'
  write(*,*)'-1= kept all SNP fixed or segragating'
  WRITE(*,*)' 0= remove non segregating SNP'
  WRITE(*,*)' 1= all below a minimum MAF'
  WRITE(*,*)' 2= keep a number of segregating SNP with top MAF. (or less if not more segregating)'
  WRITE(*,*)' 3= keep a number of segregating SNP chosen at random from a given MAF theshold (or less if not more segregating)'

  READ( *,*)eliminationType
  IF(eliminationType < -1 .OR. eliminationType > 3) eliminationType=0 ! default
  write(*,*)eliminationType

  minMAF=0.d0
  IF(eliminationType == 1)then
     WRITE(*,*)' minimum MAF for leaving the SNP (between 0:0.5)'
     READ(*,*)minMAF
     WRITE(*,*)minMAF
     IF(minMAF > 0.5d0 .OR. minMAF <0.0d0) then
        WRITE(*,*) ' wrong threshold for min MAF ', minMAF
        stop
     endif
  ELSEIF(eliminationType==2)then
     WRITE(*,*)' number of SNP to be kept'
     READ(*,*)nkept
     WRITE(*,*)nkept

  ELSEIF(eliminationType==3)then
     WRITE(*,*)' number of SNP to be kept'
     READ(*,*)nkept
     WRITE(*,*)nkept



     WRITE(*,*)' minimum MAF for leaving the SNP (between 0:0.5)'
     READ(*,*)minMAF
     WRITE(*,*)minMAF
     IF(minMAF > 0.5d0 .OR. minMAF <0.0d0) then
        WRITE(*,*) ' wrong threshold for min MAF ', minMAF
        stop
     endif


  endif

  call askInteger ( nfinaloff, ' Number of offsprings for the last generation:')


  !============================================
  ! Allocating arrays and setting starting values
  !============================================


  ! the pedigree for the historical population
  !--------------------------------------------

  nsires=npairs
  ndams =npairs
  nanim =npairs*2

  ALLOCATE(pedigree( nanim,2))
  allocate(finalPedigree(npairs*nfinalOff,2))
  !pedigree for next generation
  i=0
  a1=-1
  a2=0
  id=0
  do k=1,npairs
     a1=a1+2  !father   !odd are males parents
     a2=a2+2  !mother   !even are females parents
     id=id+1
     pedigree(id,1)=a1
     pedigree(id,2)=a2

     id=id+1
     pedigree(id,1)=a1
     pedigree(id,2)=a2
  end do

  ! same loop as above but now for finalPedigree
  do k = 1 , npairs
     do i = 0, nfinalOff - 1
        a1 = (k-1) * nfinaloff + i + 1
        finalpedigree(a1,1) = 2 * k - 1
        finalpedigree(a1,2) = 2 * k
     end do
  end do

  WRITE(*,*)' total parents ',a1,a2



  !-----------------------------------------------------------------------------
  ! the cumulative probabilities for number of recombination and mutation events
  !-----------------------------------------------------------------------------

  ALLOCATE(chiasmacumP( 100) )     !upto 100 max
  ALLOCATE(mutationcumP(100) )     !upto 100 max
  ALLOCATE(totalchiasma(0:100))
  ALLOCATE(totalmutation(0:100))
  ich=0
  imu=0
  imuM=0
  ichM=0
  totalmutation(:)=0
  totalchiasma(:)=0

  ! calculating the cumulative distribution for number of recombination events
  maxchiasma=0
  call poissonProb(chrL,maxchiasma,chiasmacumP0)
  VAL1=chiasmacumP0
  WRITE(*,*)' # recom',maxchiasma,chiasmacumP0,val1
  do
     maxchiasma=maxchiasma+1
     call poissonProb(chrL,maxchiasma,val2)
     val1=val1+val2
     IF(val1 > 1.0) val1 =1.d0
     chiasmacumP(maxchiasma) = val1

     WRITE(*,*)' # recom',maxchiasma,val2,val1
     IF(ABS(1.d0-val1) <=  0.0000001) then
        chiasmacumP(maxchiasma) = 1.d0
        exit
     endif
  end do
  WRITE(*,*)

  allocate(mutationcumP0)  !becuase it was ocnverted topointer
  ! calculating the cumulative distribution for number of mutation events
  maxmutations=0
  val3=mutationRate*nloci
  WRITE(*,*)' average number of mutations',val3
  call poissonProb(val3,maxmutations,mutationcumP0)
  VAL1=mutationcumP0
  WRITE(*,*)' # mutations',maxmutations,mutationcumP0,val1
  do
     maxmutations=maxmutations+1
     call poissonProb(val3,maxmutations,val2)
     val1=val1+val2
     IF(val1 > 1.0) val1 =1.d0
     mutationcumP(maxmutations) = val1

     WRITE(*,*)' # mutations',maxmutations,val2,val1
     IF(ABS(1.d0-val1) <=  0.0000001) then
        mutationcumP(maxmutations) = 1.d0
        exit
     endif
  end do
  WRITE(*,*)

  !=================================================================
  ! starting sampling the historical generations
  !=================================================================
  mutotal=0
  chtotal=0
  do igen =1, ngen
     IF(MOD(igen,icheck) ==1) WRITE(*,*)'igen',igen,ich,imu

     do ichr =1, nChr


        Parentgen     => ParentGenome(ichr)%genotypes
        Offgen        =>    OffGenome(ichr)%genotypes
	nloci=ParentGenome(ichr)%nloci
	nblock=ParentGenome(ichr)%nblock
        do id=1,nanim
           do igam =1,2
              Parent=pedigree(id,igam)

              IF(withpos ==1) then
                 i=1 ! no loci with same pos  (if there
                 call sampleGamete(parent,id,igam,parentGen,ParentGenome(ichr)%nloci,ParentGenome(ichr)%nblock,maxchiasma,chiasmaCumP0,chiasmaCumP,istore,nseed, &
                      offGenInp=offGen, positions=ParentGenome(ichr)%positions,samepos=i)
              else
                 call sampleGamete(parent,id,igam,parentGen,ParentGenome(ichr)%nloci,ParentGenome(ichr)%nblock,maxchiasma,chiasmaCumP0,chiasmaCumP,istore,nseed, &
                      offGenInp=offGen)
              endif
              call sampleMutation(id,igam,offgen,ParentGenome(ichr)%nloci,ParentGenome(ichr)%nblock,maxmutations,mutationCumP0,mutationCumP,istore,nseed)

              !===========================================
              totalmutation(imu)=totalmutation(imu)+1
              totalchiasma(ich)=totalchiasma(ich)+1
              mutotal=mutotal+imu
              chtotal=chtotal+ich
              IF(imu > imuM) imuM=imu
              IF(ich > ichM) ichM=ich
              !===========================================

           end do
        end do
     end do
     IF(MOD(igen,icheck) ==1)WRITE(*,*)'total chiasma mutation',chtotal,mutotal,ichM,imuM

     ! now offspring become parents and new offspring wuill be saved in other array
     thisGenome  =>  ParentGenome
     ParentGenome  =>     OffGenome
     OffGenome  =>    thisGenome
     ! swapping pedigree (performa sampling without replacement)
     DO igam=1,2
        call randomiseInt(pedigree(:,igam),nanim,nseed)
     end do

  end do






  OffGenome => ParentGenome

  !=================================================================
  ! end sampling historcal generations
  !=================================================================


  WRITE(*,*)

  !============================================
  ! saving information from historical generations
  !============================================

  OPEN(1,FILE="events.txt", STATUS='unknown')

  originalnloci=nloci
  originalnblock=nblock
  ALLOCATE(LD(maxloci))
  ALLOCATE(eliminated(maxloci))
  ALLOCATE(oldnumber( maxloci))
  IF(eliminationType > 0) then
     ALLOCATE(MAF(      maxloci))
     ALLOCATE(MAFsorted(maxloci))
  endif



  do ichr=1,nChr
     WRITE(*,*)
     WRITE(*,*)' processing Chr ', ichr
     parentGen=> ParentGenome(ichr)%genotypes
     OffGen   =>    OffGenome(ichr)%genotypes
     nblock=OffGenome(ichr)%nblock
     nloci=OffGenome(ichr)%nloci
     WRITE(*,*)' ichr OffGenome(ichr)%nblock, originalnblock ', ichr,nblock, originalnblock
     !---------------------------------------------
     WRITE(*,*)" calculating frequencies",nloci, nblock,size(frequency, dim=1)
     call CalcFrequency(nanim,OffGen,nloci,nblock,istore,frequency,k)
     WRITE(*,*) ' segregating loci',k

     !---------------------------------------------
     WRITE(*,*) ' eliminating SNP which are fixed'


     if(eliminationType == -1) then
        do i = 1, nloci
           eliminated(i)=0
        enddo
        j=nloci
     else
        j=0
        IF(eliminationType==0) then
           do i=1,nloci
              eliminated(i)=1

              IF(frequency(i,1) > 0.d0 .and. frequency(i,2) > 0.d0) then
                 eliminated(i)=0
                 j=j+1
              endif
           enddo

        ELSEIF(eliminationType==1) then
           do i=1,nloci
              eliminated(i)=1

              IF(MAF(i) > 0.d0 .and. MAF(i) >= minMAF ) then
                 eliminated(i)=0
                 j=j+1
              endif
           enddo
        ELSEIF(eliminationType==2) then
           call  SORTDX(nloci,MAF,MAFsorted)
           eliminated(:)=1
           j=0
           k=nloci
           do WHILE(j < nkept)
              i=MAFsorted(k)
              IF(MAF(i) == 0.d0 .OR. MAF(i) < minMAF) exit
              !if it reaches here then the SNP is chosen
              eliminated(i)=0
              j=j+1
              k=k-1
           enddo
        ELSEIF(eliminationType==3) then
           write(*,*)'e liination option 3'
           eliminated(:)=1
           j=0
           do i=1,nloci
              eliminated(i)=1

              IF(MAF(i) > 0.d0 .and. MAF(i) >= minMAF ) then
                 j=j+1
                 MAFsorted(j)=i
              endif
           enddo
           call randomiseInt(MAFsorted,j,nseed)
           k=min(j,nkept)
           j=0
           do i=1,k
              eliminated( MAFsorted(i) )=0
              j=j+1
           enddo
        endif
     endif
     WRITE(*,*)' number of loci to eliminate =',nloci-j
     call eliminateLoci(nanim,OffGen,nloci,nblock,istore,eliminated,oldnumber,nremained,newnblock)
     WRITE(*,*)' number of loci started and blocks started',nloci,nblock
     WRITE(*,*)' number of loci left and blocks',nremained,newnblock
     OffGenome(ichr)%nloci=nremained
     OffGenome(ichr)%nblock=newnblock



     !---------------------------------------------
     nloci  = nremained
     nblock = newnblock

     WRITE(*,*) ' Calculating LD with segregating loci'

     call calculateLD(nanim,OffGen,OffGenome(ichr)%nloci,OffGenome(ichr)%nblock,istore,LD,averLD)
     WRITE(*,*) ' average LD ',averLD

     WRITE(*,*)' calculating frequencies with segregating loci'
     call CalcFrequency(nanim,OffGen,OffGenome(ichr)%nloci,OffGenome(ichr)%nblock,istore,frequency,k)
     WRITE(*,*) ' segregating loci',k

     WRITE(filename,'(a14,i3.3)')"SNPsummary.txt",ichr
     OPEN(3,FILE=filename,STATUS='UNKNOWN')
     WRITE(3,*)' n oldposi freq ld'
     do i=1,OffGenome(ichr)%nloci
        WRITE(3,1000) i,oldnumber(i),frequency(i,1),LD(i)
     end do
     close(3)
1000 FORMAT(i7,i7,f8.5,f8.5)

     WRITE ( filename, '(a14,i3.3)' ) "SNPpositions.txt", ichr
     OPEN ( 3, FILE = filename, STATUS = 'UNKNOWN' )
     WRITE ( 3, * ) ' n position oldposi freq ld'
     do i = 1, OffGenome(ichr)%nloci
        WRITE ( 3, 1001 ) i, OffGenome(ichr)%positions( oldnumber( i ) ), oldnumber( i ), frequency( i, 1 ), LD( i )
     end do
     close ( 3 )
1001 FORMAT ( i7, f15.7, i7, f8.5, f8.5 )

     !============================================
     WRITE(*,*) ' writting the SNP'

     iline=1
     iun=2
     header=1
     WRITE(filename,'(a17,i3.3)')"genotypestest.txt",ichr
     OPEN(iun,FILE=filename, STATUS="UNKNOWN")

     call writeLoci(nanim,OffGen,nloci,nblock,istore,iun,iline,header)
     CLOSE(iun)

     WRITE(*,*) ' saving  the SNP',nloci,nblock,nanim
     WRITE(filename,'(a16,i3.3)')"genotypesInt.txt",ichr
     OPEN(iun,FILE=filename, STATUS="UNKNOWN")

     call saveLoci(nanim,OffGen,nloci,nblock,iun)
     CLOSE(iun)

  enddo
  !==============================================

  ! Here, we let sires and dams mate again to produce `nfinaloff` offsprings
  write(*,*) 'Last generation mating for',nfinaloff, 'offsprings'
  genStart = 3
  nAnim = nPairs * nFinalOff
  filename = filename(1:16)
  write(*,*) 
  call initialiseGenotypes(nChr,nAnim,genStart,nLoci, nBlock, iStore,genome1,nseed,maxloci, maxblock, ifail, filename)
  do ichr = 1,nchr

     parentGen => parentGenome(ichr)%genotypes
     offGen => offGenome(ichr)%genotypes
     nloci = parentGenome(ichr)%nloci
     nblock = parentGenome(ichr)%nblock
     do id = 1, npairs * nfinalOff
        do igam = 1, 2
           parent = finalpedigree(id, igam)
           if(withpos == 1) then
              call sampleGamete(parent,id,igam,parentGen,ParentGenome(ichr)%nloci,ParentGenome(ichr)%nblock,maxchiasma,chiasmaCumP0,chiasmaCumP,istore,nseed, &
                   offGenInp=offGen, positions=ParentGenome(ichr)%positions,samepos=i)
           else
              call sampleGamete(parent,id,igam,parentGen,ParentGenome(ichr)%nloci,ParentGenome(ichr)%nblock,maxchiasma,chiasmaCumP0,chiasmaCumP,istore,nseed, &
                   offGenInp=offGen)
           end if
        end do
     end do
  end do

  write(*,*) 'Success!'


  call ifinal(seed,seedfile)

end program generations



