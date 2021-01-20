!     Last change:  R     6 Oct 2014   11:27 am

!module sampling
! use precision_types
! USE allf77
!implicit none
!
!  integer , DIMENSION(:)    , allocatable     :: totalchiasma  ! cumulative probability of i recombination
!  integer , DIMENSION(:)    , allocatable     :: totalmutation  ! cumulative probability of i recombination
!  INTEGER :: ich,imu
!  INTEGER :: ichm,imum
!
!contains
!include 'SNP_map.f90'
!include 'poissonProb.f90'
!INCLUDE 'sampleGamete.f90'
!include 'randomise.f90'
!include 'printingDateTime.f90'
!
!end module
!

program generations

  USE sampling
  use quickSort

  implicit none

  INTEGER                                  :: istore  ! how SNP are stored
  INTEGER                                  :: withpos

  CHARACTER(LEN=11) :: seedfile
  INTEGER                                  :: nseed
  CHARACTER ( LEN = 256 ) :: pedfile,filename, filename1, filename2

  integer :: nchr,ichr
  TYPE(chromosome), DIMENSION(:), ALLOCATABLE, TARGET  :: genome1,genome2
  TYPE(chromosome), DIMENSION(:),              POINTER :: Parentgenome,Offgenome, thisGenome


  INTEGER          , DIMENSION(:,:,:), pointer      :: parentGen
  INTEGER          , DIMENSION(:,:,:), pointer      :: offGen      ! if the genotyps of offspring are store in a different array
  INTEGER          , DIMENSION(:,:,:), pointer      :: this        ! a working pointer to sawp the arrays
  REAL             , DIMENSION(:    ), allocatable, target  :: positions   ! if the loci are not equidistance, then the position (it does not matter if Morgan or cM).


  INTEGER   :: nloci,nblock    !nblock and nloci are the same if storage is as integer
  INTEGER   :: nremained,newnblock,originalnblock, originalnloci,maxblock, maxloci 

  INTEGER , DIMENSION(:,:), allocatable :: pedigree, finalPedigree
  INTEGER , DIMENSION(:  ), allocatable :: sires,dams
  INTEGER   :: npairs
  INTEGER   :: nsires,ndams,nanim,ngen,noff,nfinaloff


  DOUBLE PRECISION :: chrL, mutationRate

  INTEGER          :: maxchiasma, maxmutations
  INTEGER          :: nchiasma
  DOUBLE precision :: chiasmacumP0                                       ! probability of no recombination
  DOUBLE precision , DIMENSION(:)    , allocatable     :: chiasmacumP    ! cumulative probability of i recombination

  !  DOUBLE PRECISION :: mutationcumP0                                      ! probability of no mutation
  !  DOUBLE precision , DIMENSION(:)    , allocatable     :: mutationcumP   ! cumulative probability of i mutations

  DOUBLE PRECISION , pointer:: mutationcumP0                                      ! probability of no mutation
  DOUBLE precision , DIMENSION(:)    , pointer     :: mutationcumP   ! cumulative probability of i mutations

  DOUBLE precision , DIMENSION(:)    , allocatable, target     :: AllChrmutationcumP0  ! probability of no mutation
  DOUBLE precision , DIMENSION(:,:)  , allocatable, target     :: AllChrmutationcumP   ! cumulative probability of i mutations


  DOUBLE precision , DIMENSION(:,:)  , allocatable     :: frequency  ! frequency of both alleles
  INTEGER          , DIMENSION(:)    , allocatable     :: eliminated ! flag if SNP will be eliminated
  INTEGER          , DIMENSION(:)    , allocatable     :: oldnumber  !

  DOUBLE precision , DIMENSION(:)    , allocatable     :: LD  ! LD between consecutive loci
  DOUBLE PRECISION :: averLD

  DOUBLE precision , DIMENSION(:)    , allocatable     :: MAF  ! frequency of both alleles
  INTEGER          , DIMENSION(:)    , allocatable     :: MAFsorted  !
  INTEGER :: eliminationType, nkept, genstart
  DOUBLE PRECISION :: minMAF

  DOUBLE PRECISION :: val1,val2,val3
  INTEGER :: i,j,k,igen,mutotal,chtotal,a1,a2,iloci,iun,iline,header
  INTEGER :: parent, id,igam,iblock
  INTEGER :: iparent

  INTEGER ,DIMENSION(8):: clock_elements,end_clock

  integer :: icheck,ifail
  REAL :: rand

  character(len=1) :: theStatus, expectStatus 



  iun=6  ! the
  j=1
  call printingDateTime(iun,j,clock_elements)

  ! in case the seed file  does not exist, then starting seed calculated as function of clock time
  nseed= clock_elements(1)*10000+clock_elements(2)*100 + clock_elements(3) + &      !the date
       clock_elements(5)*10000+clock_elements(6)*100 + clock_elements(7) + &      !the time (h/m/s)
       clock_elements(4)+clock_elements(8)                                        !time millisecond

  seedfile="inicio.dat"
  CALL  istart(nseed,seedfile)

  !call askFilename ( pedfile, "pedigree file", theStatus, expectStatus )

  !=========================================================
  ! passing parameters
  !=========================================================

  call askInteger ( npairs, ' number of mating pairs (n is twice the mating pairs)' )
  nsires=npairs
  ndams =npairs
  nanim =npairs*2

  call askInteger ( ngen, ' number of generations (final gen = n + 1)' )

  call askInteger ( nchr, ' number of chromosomes' )

  WRITE ( *, * ) ' chromosome length (all equal)'
  READ ( *, * ) chrL
  WRITE ( *, * ) chrL

  WRITE ( *, * ) ' mutation rate '
  READ ( *, * ) mutationRate
  WRITE ( *, * ) mutationRate


  write ( *, '(a)' )' Start of initial generation'
  write ( *, '(a)' )'     1. all SNP fixed '
  write ( *, '(a)' )'     2. all SNP with freq 0.5 in HWE and without LD (LD will be created as generation progress)'
  write ( *, '(a)' )'     3. SNP genotype given by user (from previous run)'
  write ( *, '(a)' )' options 1 and 2 asume equal number of equidistant SNP in each chromosome'
  read(*,*) genstart
  if(genstart < 0 .or. genstart > 3) genstart=1  ! default is 1
  write(*,*) genstart




  if(genstart == 3)then
     write(*,'(a)')' prefix for genotype file name (genotype for each chr will be in file named <prefix+XXX>)'
     read(*,*)filename1
     write ( *, '(a)' ) ' prefix for position file name (position for each chr will be in file named <prefix+XXX>)'
     read ( *, * ) filename2
     istore=1    !only this option for this one
     withpos=1  ! only this option for this one
  else

     call askInteger (nloci, ' number of markers')
     !     WRITE ( *, * ) ' number of markers'
     !     READ ( *, * ) nloci
     !     WRITE ( *, * ) nloci

     istore=1
     WRITE(*,*)' storage (binary =1; integer =0)'
     READ( *,*) istore
     IF(istore .NE. 0)istore=1
     WRITE(*,*) istore

     withpos=0
     WRITE(*,*)' using position variable (Y=1,N=0)'
     read ( *,*)withpos
     IF(withpos /= 1) withpos=0
     WRITE(*,*) withpos

  endif
  ALLOCATE ( genome1( nChr ) )
  ALLOCATE ( genome2( nChr ) )
  Parentgenome => genome1
  Offgenome    => genome2

  write(*,*)' get into initiailise genotypes',genstart
  call initialiseGenotypes(nchr,nanim,genstart,nloci, nblock, istore,genome1,nseed,maxloci, maxblock, ifail, filename1)
  write(*,*)'maxloci maxblock',maxloci, maxblock

  !maxloci=0
  !maxblock=0

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
!!!!nkept=nloci
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
     ! IF(nkept > nloci .OR. nkept <0) then
     !   WRITE(*,*) ' wrong threshold for number of SNP to be kept ', nkept
     !   stop
     ! endif

  ELSEIF(eliminationType==3)then
     WRITE(*,*)' number of SNP to be kept'
     READ(*,*)nkept
     WRITE(*,*)nkept
     ! IF(nkept > nloci .OR. nkept <0) then
     !   WRITE(*,*) ' wrong threshold for number of SNP to be kept ', nkept
     !   stop
     ! endif
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



  !-----------------------------------------------------------------------------
  ! the array for the genoypes of parenst and offsprings
  !-----------------------------------------------------------------------------
  !ALLOCATE ( genome1( nChr ) )
  !ALLOCATE ( genome2( nChr ) )
  !Parentgenome => genome1
  !Offgenome    => genome2
  !if( genstart < 3) then
  !  IF(istore == 1) then
  !    i=32
  !    nblock= nloci/i
  !    IF(MOD(nloci,i) > 0) nblock=nblock+1
  !  else
  !    nblock=nloci
  !  endif
  !  maxblock=nblock
  !  maxloci=nloci
  !  do ichr = 1, nchr
  !!    genome1( ichr )%nloci  = nloci
  !    genome2( ichr )%nloci  = nloci
  !    genome1( ichr )%nblock = nblock
  !    genome2( ichr )%nblock = nblock
  !    genome1( ichr )%chrL = chrL
  !    genome2( ichr )%chrL = chrL
  !    ALLOCATE ( genome1( ichr )%genotypes( nanim, 2, nblock ) )
  !    ALLOCATE ( genome2( ichr )%genotypes( nanim, 2, nblock ) )
  !  end do
  !!ALLOCATE(parentGen(nanim,2,nblock) )
  !!ALLOCATE(offGen(   nanim,2,nblock) )
  !
  !!  Parentgenome => genome1
  !  Offgenome    => genome2
  !
  !  ALLOCATE(frequency(nloci,2))
  !  ALLOCATE(positions(nloci))
  !  frequency(:,:)=0.d0
  !  val1 = chrL/DBLE(nloci-1)
  !  positions(1)=0.d0
  !  do i=2,nloci
  !    positions(i)=DBLE(i-1)*val1
  !  end do
  !  do ichr = 1, nchr
  !    genome1( ichr )%positions => positions  !all chrmosomes have the same (same equidistant number SNP at the same position)
  !    genome2( ichr )%positions => positions  ! so all chr can share the same array of positions
  !  enddo
  !
  !  WRITE(*,*)' nblock ',nblock
  !  k = 0
  !  IF ( istore == 1 ) then
  !    do j = 31, 0, - 1
  !      k = ibclr( k, j )
  !    end do
  !    WRITE ( *, 10001 ) ' starting hapl ', k, k
  !  10001 FORMAT ( a15, i20, ' ( ', b32.32, ' ) ' )
  !  else
  !    k = 1  !genotype are 1
  !  end if
  !  if(genstart==1) then
  !    do ichr =1, nchr
  !      genome1(ichr)%genotypes(:,:,:)=k !all animalss are homozygous for all genes
  !      genome2(ichr)%genotypes(:,:,:)=k
  !    end do
  !  else
  !    ! if the initial genotype are segregating at freq 0.5 and in LE across SNPs
  !
  !    ! when genotypes are store in bits
  !    if(istore==1) then
  !      do ichr = 1, nchr
  !        do iblock = 1, nblock
  !          do igam = 1, 2
  !            do id = 1, nanim
  !              a1=k
  !              do j = 31, 0, - 1
  !                call ran1 ( nseed, rand )
  !                if ( rand > 0.5 ) a1 = ibset( a1, j )
  !              end do
  !              genome1( ichr )%genotypes( id, igam, iblock ) = a1
  !            enddo
  !          enddo
  !        enddo
  !      enddo
  !    else
  !      !  when genotypes are store as integer
  !      a1=1
  !      do ichr = 1, nchr
  !        do iblock = 1, nblock
  !          do igam = 1, 2
  !            do id = 1, nanim
  !              call ran1 ( nseed, rand )
  !              if ( rand > 0.5 ) a1 = 3-a1  !change the allele (it is ok as frequency is the same of both alleles)
  !              genome1( ichr )%genotypes( id, igam, iblock ) = a1
  !            end do
  !          end do
  !        end do
  !      end do
  !    endif
  !  endif
  !
  !else
  !  maxblock=0
  !  maxloci=0
  !  do ichr=1,nChr
  !    WRITE(filename,'(a,i3.3)')trim(filename1),ichr
  !    OPEN(newunit=iun,FILE=filename,STATUS='old')
  !    read(iun,*) i, nloci,nblock
  !	if(i <= nanim) then
  !	  write(*,*)' error. Genotype file does not have enough individuals', trim(filename)!
  !	  stop
  ! 	endif
  !	if(nloci > maxloci) maxloci=nloci
  !	if(nblock > maxblock) maxblock=nblock
  !
  !    genome1( ichr )%nloci  = nloci
  !    genome2( ichr )%nloci  = nloci
  !    genome1( ichr )%nblock = nblock
  !    genome2( ichr )%nblock = nblock
  !    genome1( ichr )%chrL = chrL
  !    genome2( ichr )%chrL = chrL
  !    ALLOCATE ( genome1( ichr )%genotypes( nanim, 2, nblock ) )
  !    ALLOCATE ( genome2( ichr )%genotypes( nanim, 2, nblock ) )
  !    rewind(iun)
  !    call openLoci ( nanim, genome1( ichr )%genotypes, nloci, nblock, iun )!
  !	close(iun)
  !    ALLOCATE ( genome1( ichr )%positions( nLoci ) )
  !    ALLOCATE ( genome2( ichr )%positions( nloci ) )
  !    genome2( ichr )%positions => genome1( ichr )%positions
  !	WRITE(filename,'(a,i3.3)')trim(filename2),ichr
  !    OPEN(newunit=iun,FILE=filename,STATUS='old')
  !    read(iun,*) !skiping the header'
  !    do i=1,nloci
  !      read(iun,*) j, rand
  !	  genome1(ichr)%positions(i)=rand
  !    enddo
  !  enddo
  !endif
  !

  !=================================================================
  ! starting sampling the historical generations
  !=================================================================
  mutotal=0
  chtotal=0
  do igen =1, ngen
     IF(MOD(igen,icheck) ==1) WRITE(*,*)'igen',igen,ich,imu

     do ichr =1, nChr
        !!    IF ( MOD( igen, icheck ) == 1 ) WRITE ( *, * ) 'igen chr', igen, ichr
        Parentgen     => ParentGenome(ichr)%genotypes
        Offgen        =>    OffGenome(ichr)%genotypes
	nloci=ParentGenome(ichr)%nloci
	nblock=ParentGenome(ichr)%nblock
        do id=1,nanim
           do igam =1,2
              Parent=pedigree(id,igam)
              !WRITE(*,*)' id, igam parent',id,igam,parent
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
     !  this      => ParentGen
     !  ParentGen => offGen
     !  offGen    => this
     thisGenome  =>  ParentGenome
     ParentGenome  =>     OffGenome
     OffGenome  =>    thisGenome
     ! swapping pedigree (performa sampling without replacement)
     DO igam=1,2
        call randomiseInt(pedigree(:,igam),nanim,nseed)
     end do
     !  IF ( MOD( igen, icheck ) == 1 ) WRITE ( *, * ) 'end of loop',igen

  end do






  !offGen => ParentGen
  OffGenome => ParentGenome

  !=================================================================
  ! end sampling historcal generations
  !=================================================================


  WRITE(*,*)

  !============================================
  ! saving information from historical generations
  !============================================

  OPEN(1,FILE="events.txt", STATUS='unknown')

  !---------------------------------------------
  ! number of recombinations
  !---------------------------------------------
  !k=SUM(totalchiasma)    ! number of meiosis
  !j=totalchiasma(0)      !number meiosis without recombination
  ! WRITE(1,*)" # number_recomb number_meiosis probab_recombin_numb cummul_meiosis total_meiosis"
  ! WRITE(1,*)0,totalchiasma(0),chiasmaCumP0,j,k
  ! do i=1,maxchiasma
  !   j=j+totalchiasma(i)
  !   WRITE(1,*)i,totalchiasma(i), chiasmaCumP(i),j,k
  ! end do
  ! WRITE(1,*)
  ! WRITE(1,*)

  !---------------------------------------------
  ! number of mutations
  !---------------------------------------------
  !WRITE(1,*)" # number_mutation number_meiosis probab_mutation_numb cummul_meiosis total_meiosis"
  !k=SUM(totalmutation)
  !j=totalmutation(0)
  !    WRITE(1,*)0,totalmutation(0),mutationCumP0,j,k
  !do i=1,maxmutations
  !  j=j+totalmutation(i)
  !    WRITE(1,*)i,totalmutation(i),mutationCumP(i),j,k
  !end do





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

     !    WRITE(1,*)
     !    WRITE(1,*)
     !    WRITE(1,*)" #allele frequencies"
     !    k=0    !number of segregating snp
     !    igen=0 !number of snp with MAF > 0.05
     !    do i=1,OffGenome(ichr)%nloci
     !      j=0
     !      IF(frequency(i,1) > 0.d0 .and. frequency(i,2) > 0.d0) j=1
     !      IF(frequency(i,1) >= 0.05d0 .and. frequency(i,2) >= 0.05d0) igen=igen+1
     !      k=k+j
     !      IF(eliminationType > 0) MAF(i)= MIN(frequency(i,1), frequency(i,2))
     !      WRITE(1,*)i,frequency(i,1), frequency(i,2),j,k,igen
     !    end do
     !    WRITE(*,*) ' segregating loci (total / MAF >= 0.05) ',k,igen


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
     !   WRITE(1,*)
     !   WRITE(1,*)
     !   WRITE(1,*)" #LD for segregaing SNP"
     !   k=0
     !   do i=1,OffGenome(ichr)%nloci
     !     j=0
     !     WRITE(1,*)i,LD(i)
     !   end do

     WRITE(*,*)' calculating frequencies with segregating loci'
     call CalcFrequency(nanim,OffGen,OffGenome(ichr)%nloci,OffGenome(ichr)%nblock,istore,frequency,k)
     WRITE(*,*) ' segregating loci',k

     !    OPEN(3,FILE="SNPsummary.txt",STATUS='UNKNOWN')
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
     !    OPEN(iun,FILE="genotypestest.txt", STATUS="UNKNOWN")
     call writeLoci(nanim,OffGen,nloci,nblock,istore,iun,iline,header)
     CLOSE(iun)

     WRITE(*,*) ' saving  the SNP',nloci,nblock,nanim
     WRITE(filename,'(a16,i3.3)')"genotypesInt.txt",ichr
     OPEN(iun,FILE=filename, STATUS="UNKNOWN")
     !    OPEN(iun,FILE="genotypesInt.txt",STATUS="unknown")
     call saveLoci(nanim,OffGen,nloci,nblock,iun)
     CLOSE(iun)

  enddo
  !==============================================



  !WRITE(*,*) 'expanding'
  !ParentGen => OffGen
  !DEALLOCATE(OffGen)
  !DEALLOCATE(pedigree)
  !nanim = npairs * noff


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


  call ifinal(nseed,seedfile)
  iun=6  ! the
  j=1
  call printingDateTime(iun,j,end_clock)
  call printingElapseTime(iun,clock_elements,end_clock)

  stop


end program generations



