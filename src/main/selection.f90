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
  character(len=100) :: baseNameFreq, outputfile
  integer :: iunoutput
  character(len=60) :: formato
  type(chromosome), dimension(:), allocatable, target :: genome1, genome2
  type(chromosome), dimension(:), pointer :: Parentgenome, Offgenome, thisGenome
  integer, dimension(:), allocatable :: seed
  real(KINDR) :: chrL, mutationRate
  integer, dimension(:), pointer :: indiv, male, female
  logical, dimension(:), pointer :: sex !true = male, false = female
  integer, dimension(:,:,:), pointer :: parentGen, offGen!, this

  logical :: doreml, randomQTL, reactionNorm

  real(KINDR), dimension(2) :: interval
  real(KINDR) :: farmRange, maf
  real(KINDR), dimension(:,:), allocatable :: locations, farmBounds, X !incidence matrix
  integer, dimension(:), allocatable :: farmInd
  integer :: nlox, nfarm, allocation

  !ncomp is the number of traits by a QTL (slope, intercept --> 2, otherwise 1)
  integer :: nSNP, nQTL, nComp
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array) :: QTLlist

  real(KINDR), allocatable, dimension(:,:) :: TBV, CTE
  real(KINDR), allocatable, dimension(:) :: Amat, means, phenotypes, theta, fixeff
  real(KINDR), allocatable, dimension(:) :: chiasmaCumP, mutationCumP
  real(KINDR), allocatable, dimension(:) :: farmIndReal, farmEffects
  real(KINDR) :: mutationCumP0, chiasmaCumP0
  type(doublePre_Array), dimension(:), allocatable :: raneff
  integer, allocatable, dimension(:) :: ids
  integer, allocatable, dimension(:) :: totalChiasma, totalMutation
  integer, allocatable, dimension(:,:) :: pedigree
  integer :: ivar, iscaled, nobs, maxid, nfix, nvar, nran
  type(variances) :: vars
  integer :: n_m, n_fpm, n_opf, ngen
  integer :: selectionType, analysisType
  integer :: maxchiasma, maxmutations
  ! counters
  integer :: i, j, k, iChr, iGam, id, igen
  real(KINDR) :: val1, val2
  character(len=30) :: status, estatus

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
     write(STDERR, '(a,a)') "error in openning input file ", trim(inputfile)
     stop 2
  end if

  call readInput(inputfile, verbose, nanim, nchr, filename1, filename2,&
       chrL, mutationRate, ncomp, vars, nQTL, nSNP, MAF, baseNameFreq,&
       randomQTL, interval, locations, X, nlox, nFarm, farmBounds,&
       farmInd, farmRange, allocation, selectionType, nobs, means,&
       analysisType, theta, n_m, n_fpm, n_opf, ngen, doreml, reactionNorm,&
       nfix, nvar, nran, outputfile)

  call defineFarms(interval, nfarm, farmRange, farmBounds)

  open(1, file = 'farms.txt')
  do i = 1, nfarm
     write(1, *) farmBounds(i, 1), farmBounds(i, 2)
  end do
  close(1)

  open(newUnit = iunoutput, file = outputfile)
  allocate(indiv(nanim), sex(nanim), CTE(nComp, 2))
  indiv= (/( i, i = 1, nanim )/)
  ! the array sex corresponds to the new generation (before selection)
  sex(1:nanim) = .false. ! all female
  sex(1:nanim:n_opf) = .true. ! for each female one offspring is male
  allocate(pedigree(nanim, 3))
  i = n_m * n_fpm
  allocate(male(n_m), female(i))
  allocate(ids(nobs), phenotypes(nobs))
  allocate(genome1(nChr))
  allocate(genome2(nChr))
  allocate(TBV(nanim, ncomp))
  maxid = nanim ! maxval(indiv)
  ! ==============================================
  ! simulating first individuals from genepool
  ! ==============================================
  i = 3 ! genstart = 3
  j = 1 ! istore = 1
  k = 1 ! withpos = 1
  call initialiseGenotypes(verbose, nchr, nanim, i, nloci, nblock, j, &
       genome1, genome2, maxloci, maxblock, ifail, chrL, &
       trim(filename1), trim(filename2))
  if (verbose) write(STDOUT, '(a)') "Genotypes Initialised"
  ! ==============================================
  ! calculating the cumulative distribution for number of recombination events
  ! ==============================================
  i = 200 ! upto 200 recombinations
  if (verbose) write(STDOUT, '(a)') &
       "Setting total number of mutations and recombinations"
  call GetMutRecArray(verbose, i, chrL, mutationRate, nLoci, &
       chiasmaCumP, chiasmacumP0, totalChiasma, maxchiasma, &
       mutationCumP, mutationCumP0, totalMutation, maxmutations)
  if (verbose) write(STDOUT, '(a)') "Mutation and recombinations set"
  ! ================================================================
  ! Getting QTL and SNP list
  ! ================================================================
  if (verbose) write(STDOUT, '(a)') "Getting SNPs and their values"
  call getQTLandSNP(verbose, nChr, nQTL, nSNP, nComp, randomQTL, genome1,&
       QTLlist, SNPlist, vars%cov, baseNameFreq, maf)
  if (verbose) write(STDOUT, '(a)') "QTL and SNP list simulated"
!  open(1, file = "QTLlist.txt")
!  do ichr = 1, nchr
!     do i = 1, nQTL
!        write(1, *) QTLlist%values(ichr, i, 1:2)
!     end do
!  end do
!  close(1)
  ! ================================================================
  ! Simulating TBV
  ! ================================================================
  if (verbose) write(STDOUT, '(a,i2)') "simulating BV for generation", 0
  call SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, genome1,&
       QTLlist, SNPlist, TBV, verbose)
  if (verbose) write(STDOUT, '(a)') "breeding values simulated"
!  open(1,file = 'tbvinitial1')
!  do i = 1, nanim
!     write(1, *) (tbv(i,j), j = 1, ncomp)
!  end do
!  close(1)
 
  ! getting scaling factor for QTL 
  ! col 1 is scaling, col 2 is shifting
  do i = 1, nComp
     CTE(i, 1) = sqrt(vars%A(i) / variance(TBV(1:nanim,i), nanim))
     QTLlist%values(1:nChr, 1:nQTL, i) = QTLlist%values(1:nChr, 1:nQTL, i)*&
          CTE(i, 1)
     TBV(1:nAnim, i) = TBV(1:nAnim, i) * CTE(i, 1)
     CTE(i, 2) = sum(TBV(1:nAnim, i)) / nAnim
     TBV(1:nAnim, i) = TBV(1:nAnim, i) - CTE(i, 2)
  end do
  ! saving first TBV and the QTLlist
!  open(1,file = 'tbvinitial2')
!  do i = 1, nanim
!     write(1, *) (tbv(i,j), j = 1, ncomp)
!  end do
!  close(1)
!  open(1, file = 'qtllist')
!  do i = 1, nchr
!    do j = 1, nqtl
!        write(1, *) (qtlList%values(i,j,k), k = 1, ncomp)
!     end do
!  end do
!  close(1)

!  call SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, genome1,&
!       QTLlist, SNPlist, TBV, verbose)
!  open(1,file = 'tbvinitial3')
!  do i = 1, nanim
!     write(1, *) (tbv(i,j) - cte(j,2), j = 1, ncomp)
!  end do
!  close(1)

  ! ================================================================
  ! Simulating Phenotypes
  ! ================================================================
  if (verbose) write(STDOUT, '(a)') "allocating individuals"
  i = allocation
  allocation = 1 ! temporary allocation is set to 1 because in the first
  ! generation we do not know the pedigree therefore we can do only random
  call allocateInd(nAnim, nlox, nfarm, allocation, farmBounds, farmRange,&
       farmInd, locations)
  allocation = i
  if (verbose) write(STDOUT, '(a)') "individuals allocated"

  if (verbose) write(STDOUT, '(a)') "simulating phenotypes"
  call SimulatePhenotype(verbose, nAnim, nComp, nFix, nLox, nran, indiv,&
       TBV, vars, means, nobs, locations, ids, phenotypes, "cov", X)
  if (verbose) write(STDOUT, '(a)') "Phenotypes simulated"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !write(formato, '(a,i2.2)') "phen", 0
  !open(1, file = trim(formato))
  !do i = 1, nobs
  !   write(1, *) ids(i), X(i, 1), farmind(i), phenotypes(i)
  !end do
  !close(1)

  ! only for the first generation
  ! a very simple guestimate for theta
  if (nran == 1) then
     val1 = (interval(1) + interval(2)) / 2 ! x=x_middle
     ! nominator
     val2 = val1 * val1 * vars%A(1) + vars%A(2) + 2 * val1 * vars%cov(1,2)
     ! denominator
     val1 = val2 + val1 * val1 * vars%E(1) + vars%E(2)
     val1 = val2 / val1
     write(6, *) 'estimated heritability', val1
     theta(1) = variance(phenotypes, nobs) * val1
     theta(2) = variance(phenotypes, nobs) * (ONE - val1)
  end if
  
  ParentGenome => genome1
  Offgenome    => genome2

  if (verbose) write(STDOUT, '(a,i3)') "generation ", 0 ! must be zero
  val1 = sum(tbv(1:nanim, 1)) / nanim
  val2 = sum(tbv(1:nanim, 2)) / nanim
  write(STDOUT, 68) val1, val2
  write(iunoutput, 200) "gen", "mean(TBVs)", "mean(TBVi)", "var(TBVs)", &
       "var(TBVi)", "Est(mu_s)", "Est(mu_i)", "acc(EBVs)", "acc(EBVi)"
200 format(a3,8(1x,a15))
201 format(i3,4(1x,f15.7))
  write(iunoutput, 201, advance = 'no') 0, val1, val2, &
       variance(tbv(1:nanim, 1), nanim), variance(tbv(1:nanim, 2), nanim)
  ! initialising first generation individuals
  gen: do igen = 1, ngen
     write(STDOUT, '(a, i2)') "generation ", igen
     ! ================================================================
     ! Genomic Evaluation
     ! TODO: This block (genSel) should eventually be one single
     !       subroutine with name "doGenomicSelection"
     ! ================================================================
     genSel: select case (selectionType)
     case (1) ! random
        write(6, *) "selectionType :" , selectionType
        if (igen == 1) then
           allocate(raneff(1))
           allocate(raneff(1)%level(nanim))
        end if
        i = 1
        call random_number(raneff(i)%level(1:nanim))
        ! selectByIntercept needs raneff of size 1 (or size 2 but then effects must be on the second)
        write(iunoutput, '(4(1x,a15))') "NaN", "NaN", "NaN", "NaN"
     case (2, 3, 6) ! slope ebv or interceept ebv
        ! reaction norm means x must be re-written
        if (reactionNorm) then
           if ((selectionType.eq.2).or.(selectionType.eq.3)) then
              if (.not.allocated(farmIndReal)) &
                   allocate(FarmIndReal(nobs), farmEffects(nfarm))
              ! converting to double as leastSquare takes double
              farmIndReal(1:nobs) = dble(farmInd(1:nobs))
              call leastSquare(verbose, nobs, nfarm, ids, farmIndReal,&
                   phenotypes, farmEffects)
              ! scaling farmeffects to [xmin, xmax]
              val1 = minval(farmEffects)
              val2 = maxval(farmEffects)
              farmEffects(1:nfarm) = (val2 - farmEffects(1:nfarm)) /&
                   (val2 - val1)*(interval(2) - interval(1)) + interval(1)
              ! replacing scaled farmeffects with challenge levels
              X(1:nobs, 1) = farmEffects(farmInd(1:nobs))
           elseif (selectionType.eq.6) then
              do i = 1, nobs
                 if (farmInd(i) .eq. 1) cycle
                 X(i, farmInd(i) - 1) = ONE
              end do
              !!!!!!!!!!
              !write(formato, '(a,i2.2)') "incidence", igen - 1
              !open(1, file = trim(formato))
              !333 format(14(f3.1,1x), f3.1)
              !do i = 1, nobs
              !   write(1, 333) X(i,1:nfix) 
              !end do
              !close(1)
           end if
        end if

        ! making Gmatrix
        iscaled = 1 !(0:no, 1:yes)
        ivar  = 1 !(0:sample, 1:2pq, 2:2p'q')
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
           call reml(ids, X, phenotypes, nfix, nobs, maxid, Amat, nvar,&
                nran, theta, fixEff, ranEff, verbose, 0, 20)
        else
           call blup(ids, X, phenotypes, nfix, nobs, maxid, Amat, nvar,&
                nran, theta, fixEff, ranEff, verbose)
        end if
        if (verbose) write(STDOUT, '(a)') "Genomic evaluation done"
        ! -----------------------------------
        ! getting ebv accuracy
        ! ----------------------------------
        if ((selectionType.eq.2).or.(selectionType.eq.3)) then
           write(iunoutput, '(2(1x,f15.7))', advance = 'no') fixEff(1:nfix)
           do i = 1, nFix
              val1 = correlation(raneff(i)%level, TBV(:,i), nanim)
              write(6, *) 'accuracy', i, val1
              write(iunoutput, '(1x,f15.7)', advance = 'no') val1
           end do
           write(iunoutput, *) 
        elseif (selectionType.eq.6) then
           write(iunoutput, '(1x,f15.7,3(1x,a15))') fixEff(nfix), "NaN", "NaN", "NaN"
        end if
        !====================================
        ! selection
        !====================================
        i = merge(selectionType - 1, 1, selectionType .ne. 6)
     case (4, 5)
        write(6, *) 'selectionType :' , selectionType
        if (igen == 1) then
           allocate(raneff(3))
           allocate(raneff(1)%level(nanim), raneff(2)%level(nanim), raneff(3)%level(nobs))
        end if
        raneff(1)%level(1:nanim) = TBV(1:nanim, 1)! random effects are true values
        raneff(2)%level(1:nanim) = TBV(1:nanim, 2)! random effects are true values
        raneff(3)%level(1:nobs) = ZERO
        i = selectionType - 3
        write(iunoutput, '(4(1x,a15))') "NaN", "NaN", "NaN", "NaN"
     end select genSel
     call selectMates(nanim, indiv, sex, n_m, n_fpm, male, female,&
          ranEff(i)%level, verbose)
     if (verbose) write(STDOUT, '(a)') "selection done"
     
     ! making pedigree for next generation based on selected parents
     pedigree = makePedigree(n_m, n_fpm, n_opf, male, female, nanim)
     !write(filename1, '(a, i2.2)') "ped.G", igen
     !open(1, file = trim(filename1))
     !write(1, '(4x,a2,2x,a4,3x,a3)') "id","sire","dam"
     !do i = 1, nanim
     !   write(1, '(3i6)') (pedigree(i,j), j = 1, 3)
     !end do
     !close(1)

     ! ========================
     ! mating
     ! ========================
     i = 1 ! istore = 1
     k = 1 ! samepos = 1
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
                   samepos = k)
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
     
     ! writing all chromosomes
!     do ichr = 1, nchr
!        write(filename1, '(a,i3.3,a1,i2.2)') "chr", iChr, '.', igen
!        open(1, file = trim(filename1))
!        write(formato, '(a1, i10, a6)' ) '(', parentGenome(ichr)%nblock, 'i12)'
!        do id =1 , nanim
!           do igam = 1, 2
!              write(1, formato) offGenome(ichr)%genotypes(id, igam, 1:offGenome(ichr)%nblock)
!           end do
!        end do
!        close(1)
!     end do
 
     ! ===========================================
     ! swapping pointers for offspring and parents
     ! ===========================================    
     thisGenome  =>  ParentGenome
     ParentGenome  => OffGenome
     OffGenome  => thisGenome
     
     ! ===========================================
     ! simulating data for the new generation
     ! ===========================================
     ! getting breeding values
     if (verbose) write(STDOUT, '(a,i2)') "simulating BV for generation", igen
     call SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, ParentGenome,&
          QTLlist, SNPlist, TBV, verbose)
     if (verbose) write(STDOUT, '(a)') "breeding values simulated"

     ! shifting true breeding values
     do i = 1, nComp
        TBV(1:nAnim, i) = TBV(1:nAnim, i) - CTE(i, 2)
     end do
!     write(filename1, '(a,i2.2)') "TBV.G", igen
!     open(1, file = trim(filename1))
!     do i = 1, nanim
!        write(1, *) (tbv(i,j), j = 1, ncomp)
!     end do
!     close(1)

     if (verbose) write(STDOUT, '(a)') "allocating individuals"
     call allocateInd(nAnim, nlox, nfarm, allocation, farmBounds, &
          farmRange, farmInd, locations, pedigree = pedigree, nm = n_m,&
          male = male)
     if (verbose) write(STDOUT, '(a)') "individuals allocated"

     ! getting phenotypes
     if (verbose) write(STDOUT, '(a)') "simulating phenotypes"
     call SimulatePhenotype(verbose, nAnim, nComp, nFix, nLox, nran, indiv,&
          TBV, vars, means, nobs, locations, ids, phenotypes, "cov", X)
     if (verbose) write(STDOUT, '(a)') "Phenotypes simulated"

     !write(filename1, '(a,i2.2)') "phen", igen
     !open(1, file = trim(filename1))
     !do i = 1, nobs
     !   write(1, *) ids(i), X(i, 1), farmind(i), phenotypes(i)
     !end do
     !close(1)

     write(6, 68) sum(tbv(1:nanim,1))/nanim,sum(tbv(1:nanim,2))/nanim
68   format("slope: ", g25.14, "; intercept: ", g25.14)
     
     write(iunoutput, '(i3,4(1x,f15.7))', advance = 'no') igen, &
          sum(tbv(1:nanim,1))/nanim, sum(tbv(1:nanim,2))/nanim, &
          variance(tbv(1:nanim,1), nanim), variance(tbv(1:nanim,2), nanim)

  end do gen
  write(iunoutput, '(4(1x,a15))') "NaN",  "NaN",  "NaN",  "NaN"
  close(iunoutput)
  call ifinal(seed, startfile)
end program selection
