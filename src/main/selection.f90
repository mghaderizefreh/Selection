program selection
  use constants
  use math
  use global_module
  use rng_module
  use evolution_module
  use blup_module
  use reml_module

  implicit none

  logical :: verbose
  integer :: nanim, nloci, nblock, maxloci, maxblock, ifail, nChr
  character(len=100) :: startFile, filename1, filename2, inputfile
  character(len=60) :: formato
  type(chromosome), dimension(:), allocatable, target :: genome1, genome2
  type(chromosome), dimension(:), pointer :: Parentgenome, Offgenome, thisGenome
  integer, dimension(:), allocatable :: seed
  integer, dimension(:), pointer :: indiv, male, female
  logical, dimension(:), pointer :: sex !true = male, false = female
  integer, dimension(:,:,:), pointer :: parentGen, offGen!, this

  logical :: doreml
  real(KINDR), dimension(:,:), allocatable :: X !incidence matrix
  real(KINDR) :: chrL, mutationRate

  !ncomp is the number of traits by a QTL (slope, intercept --> 2, otherwise 1)
  integer :: nSNP, nQTL, nComp
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array) :: QTLlist

  real(KINDR), allocatable, dimension(:,:) :: TBV, locations, CTE
  real(KINDR), allocatable, dimension(:) :: Amat, means, phenotypes, theta, fixeff
  real(KINDR), allocatable, dimension(:) :: chiasmaCumP
  real(KINDR), allocatable, dimension(:) :: mutationCumP
  real(KINDR) :: mutationCumP0, chiasmaCumP0
  type(doublePre_Array), dimension(:), allocatable :: raneff
  integer, allocatable, dimension(:) :: ids
  integer, allocatable, dimension(:) :: totalChiasma, totalMutation
  integer, allocatable, dimension(:,:) :: pedigree
  integer :: ivar, iscaled, nobs, maxid, nfix, nvar
  type(variances) :: vars
  integer :: n_m, n_fpm, n_opf, ngen
  integer :: selectionType, analysisType
  integer :: maxchiasma, maxmutations!, chTotal, muTotal
  ! counters
  integer :: i, j, k, iChr, iGam, id, igen!, iCh, iMu

  character(len=30) :: status, estatus

!  ich = 0
!  imu = 0
!  mutotal = 0
!  chtotal = 0
  startfile = "inicio.dat"
  call istart(seed, startfile, ifail)
  if (ifail /= 0) then
     write(STDERR, '(a)') "reading/setting seed faild"
     stop 2
  end if
  call random_seed(put = seed)

  status = "old"
  call askFileName(inputfile, "input file: ", status, estatus)
  if (status(1:1) .eq. "x") then
     write(STDERR, *) "error in openning input file ", trim(inputfile)
     stop 2
  end if

  call readInput(inputfile, verbose, nanim, nchr, filename1, filename2,&
       chrL, mutationRate, ncomp, vars, nQTL, nSNP, locations, &
       selectionType, nobs, means, analysisType, theta, n_m, n_fpm, &
       n_opf, ngen, doreml, nfix, nvar)

  allocate(indiv(nanim), sex(nanim), CTE(nComp, 2))
  indiv= (/( i, i = 1, nanim )/)
  sex(1:nanim) = .false.
  sex(1:nanim:3) = .true.
  allocate(pedigree(nanim, 3))
  allocate(male(10), female(300))
  allocate(ids(nobs), phenotypes(nobs))
  allocate(genome1(nChr))
  allocate(genome2(nChr))
  maxid = nanim ! maxval(indiv)
  ! ==============================================
  ! simulating first individuals from genepool
  ! ==============================================
  i = 3 ! genstart = 3
  j = 1 ! istore = 1
  k = 1 ! withpos = 1
  ! read the files to initialise the genotype of the first generation
  call initialiseGenotypes(verbose, nchr, nanim, i, nloci, nblock, j, &
       genome1, genome2, maxloci, maxblock, ifail, chrL, &
       trim(filename1), trim(filename2))
  ! ==============================================
  ! calculating the cumulative distribution for number of recombination events
  ! ==============================================
  i = 200 ! upto 200 recombinations
  call GetMutRecArray(verbose, i, chrL, mutationRate, nLoci, &
       chiasmaCumP, chiasmacumP0, totalChiasma, maxchiasma, &
       mutationCumP, mutationCumP0, totalMutation, maxmutations)
  ! ================================================================
  ! Getting QTL and SNP list
  ! ================================================================
  !randomQTL = .true.
  call getQTLandSNP(verbose, nChr, nQTL, nSNP, nComp, .true., genome1,&
       QTLlist, SNPlist, vars%corr)
  !  open(1, file = "QTLlist.txt")
  !  do ichr = 1, nchr
  !     do i = 1, nQTL
  !        write(1, *) QTLlist%values(ichr, i, :)
  !     end do
  !  end do
  !  close(1)
  ! ================================================================
  ! Simulating TBV
  ! ================================================================
  ParentGenome => genome1
  Offgenome    => genome2
  ! initialising first generation individuals
  gen: do igen = 1, ngen
     write(STDOUT, '(a,i3)') "generation ", igen

     call SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, ParentGenome,&
          QTLlist, SNPlist, TBV, verbose)

     if (igen == 1) then! this shall run only for the first generation
        do i = 1, nComp
           CTE(i, 1) = sqrt(vars%A(i) / variance(TBV(1:nanim,i), nanim))
           QTLlist%values(1:nChr, 1:nQTL, i) =QTLlist%values(1:nChr, 1:nQTL, i)* &
                CTE(i, 1)
           TBV(1:nAnim, i) = TBV(1:nAnim, i) * CTE(i, 1)
           CTE(i, 2) = sum(TBV(1:nAnim, i)) / nAnim
           TBV(1:nAnim, i) = TBV(1:nAnim, i) - CTE(i, 2)
        end do

        open(1,file = 'tbvinitial')
        do i = 1, nanim
           write(1, *) (tbv(i,j), j = 1, ncomp)
        end do
        close(1)

        open(1, file = 'qtllist')
        do i = 1, nchr
           do j = 1, nqtl
              write(1, *) (qtlList%values(i,j,k), k = 1, ncomp)
           end do
        end do

     else
        do i = 1, nComp
           QTLlist%values(1:nChr,1:nQTL,i) = QTLlist%values(1:nChr, 1:nQTL, i)* &
                CTE(i, 1)
           TBV(1:nAnim, i) = TBV(1:nAnim, i) - CTE(i, 2)
        end do
     end if

     ! ================================================================
     ! Simulating Phenotypes
     ! ================================================================
     call SimulatePhenotype(verbose, nAnim, nComp, indiv, TBV, vars, means, &
          nobs, locations, ids, phenotypes, "cov", X)
     
     ! ================================================================
     ! Genomic Evaluation
     ! TODO: This block (genSel) should eventually be one single
     !       subroutine with name "doGenomicSelection"
     ! ================================================================
     genSel: select case (selectionType)
     case (1) ! random
        if (igen == 1) allocate(raneff(1))
        allocate(raneff(1)%level(nanim))
        call random_number(raneff(1)%level(1:nanim))
        ! selectBySlope needs raneff of size 1
        call selectByIntercept(nanim, indiv, sex, n_m, n_fpm, male, female, raneff)
     case (2, 3) ! slope ebv or interceept ebv
        ! making Gmatrix
        iscaled = 1 !(0:no, 1:yes)
        ivar  = 1!(0:sample, 1:2pq, 2:2p'q')
        j = 0 ! imiss (0:mean, 1:ignore)
        i = 1 ! addDom (1:additive, 2:dominance)
        call getGmatrix(nanim, nChr, nSNP, indiv, ParentGenome, SNPlist, &
             iscaled, ivar, j, i, Amat, verbose)
        ! writing AMAT to disk
        !  i = maxval(indiv)
        !  k = 1
        !  do while (i >= 10)
        !     i = i / 10
        !     k = k + 1
        !  end do
        !  open( 1, FILE = "AMAT.txt", STATUS = 'unknown' )
        !  i = nanim * ( nanim + 1 ) / 2
        !  write(formato,'(a2,i1,a5,i1,a22)')'(i',k,',1x,i',k,&
        !     ',1x,g24.15,i9,g24.15)'
        !  write(6,*) " formato= ", trim(formato)
        !  k = 0
        !  do i = 1, nanim
        !     do j = 1, i
        !        k = k + 1
        !        write(1, formato) indiv(i), indiv(j), amat(k)
        !     end do
        !  end do
        !  close(1)
        ! reading amat if required
        !  open(1, file = 'AMAT.txt')
        !  i = nanim * ( nanim + 1 ) / 2
        !  k = 0
        !  do i = 1, nanim
        !     do j = 1, i
        !        k = k +1
        !        read(1, *) nfix, nobs, amat(k)
        !     end do
        !  end do
        !  close(1)
        !=========================================
        
        if (doreml) then
           call reml(ids, X, phenotypes, nfix, nobs, maxid, Amat, nvar, &
                theta, fixEff, ranEff, verbose, 3, 10)
        else
           call blup(ids, X, phenotypes, nfix, nobs, maxid, Amat, nvar, &
                theta, fixEff, ranEff, verbose)
        end if
        if (verbose) write(STDOUT, *) "Genomic evaluation done"
        ! -----------------------------------
        ! getting ebv accuracy
        ! ----------------------------------
        do i = 1, 2
           write(6, *) 'accuracy for i', i, &
                correlation(raneff(i)%level, TBV(:,i), nanim)
        end do
        !====================================
        ! selection
        !====================================
        ! assuming ids are sorted and uniq(ids) = indiv
        if (selectionType .eq. 2) then
           call selectByIntercept(nanim, indiv, sex, n_m, n_fpm, male, &
                female, ranEff)
        else
           call selectByIntercept(nanim, indiv, sex, n_m, n_fpm, male, &
                female, ranEff)
        end if
     case (4, 5)
        write(6, *) 'not implemented yet'
        stop 1
     end select genSel
     
     ! making pedigree for next generation based on selected parents
     pedigree = makePedigree(n_m, n_fpm, n_opf, male, female, nanim)
     open(1, file = 'ped.txt')
     write(1, '(4x,a2,2x,a4,3x,a3)') "id","sire","dam"
     do i = 1, nanim
        write(1, '(3i6)') (pedigree(i,j), j = 1, 3)
     end do
     close(1)

     ! ========================
     ! mating
     ! ========================
     i = 1 ! istore = 1
     do ichr = 1, nChr
        Parentgen => ParentGenome(iChr)%genotypes
        Offgen => OffGenome(iChr)%genotypes
        nloci = ParentGenome(iChr)%nloci
        nblock = ParentGenome(iChr)%nblock
        do id = 1, nanim
           do igam = 1, 2
              j = pedigree(id, igam+1) !parent
              call sampleGamete(j, id, igam, parentGen, &
                   ParentGenome(iChr)%nloci, ParentGenome(iChr)%nblock, &
                   maxchiasma, chiasmaCumP0, chiasmaCumP, i, &
                   offGenInp = offGen, positions = ParentGenome(iChr)%positions,&
                   samepos = 1)
              call sampleMutation(id, igam, offgen, ParentGenome(iChr)%nloci,&
                   maxMutations, mutationCumP0, mutationCumP, i)
              !===========================================
              !totalmutation(imu) = totalmutation(imu) + 1
              !totalchiasma(ich) = totalchiasma(ich) + 1
              !mutotal = mutotal + imu
              !chtotal = chtotal + ich
              !===========================================
           end do
        end do
     end do
     
!     write(filename1, '(a,i2.2)') "gen", igen
!     open(1, file = filename1)
!     do ichr = 1, nchr
!        write(formato, '(a1, i10, a6)' ) '(', parentGenome(ichr)%nblock, 'i12)'
!        do id =1 , nanim
!           do igam = 1, 2
!              write(1, formato) offGenome(ichr)%genotypes(id, igam, :)
!           end do
!        end do
!     end do
!     close(1)
!     if (igen .eq. 5) stop 3
     
     thisGenome  =>  ParentGenome
     ParentGenome  =>     OffGenome
     OffGenome  =>    thisGenome
     

     write(6, 68) sum(tbv(1:nanim,1))/nanim,sum(tbv(1:nanim,2))/nanim
     
  end do gen
68 format("slope: ", g25.14, "; intercept: ", g25.14)

  open(1, file = 'tbvfinal')
  do i = 1, nanim
     write(1, *) (tbv(i,j), j = 1, ncomp)
  end do
  close(1)

!  call ifinal(seed, startfile)
end program selection
