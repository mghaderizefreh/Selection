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
  integer :: nanim, nloci, nblock, maxloci, maxblock, ifail, nChr, nelement
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

  logical :: randomQTL, reactionNorm
  integer :: varEst
  real(KINDR), allocatable, dimension(:) :: P, V
  real(KINDR), allocatable, dimension(:,:) :: Vhat, temp
  real(KINDR), dimension(2) :: interval
  real(KINDR) :: farmRange, maf
  real(KINDR), dimension(:,:), allocatable :: locations, farmBounds, X !incidence matrix
  integer, dimension(:), allocatable :: farmInd
  integer :: nlox, nfarm, allocation

  !ncomp is the number of traits by a QTL (slope, intercept --> 2, otherwise 1)
  integer :: nSNP, nQTL, nComp
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array) :: QTLlist

  integer, allocatable, dimension(:) :: chr_nlocibefore, ipiv
  real(KINDR), allocatable, dimension(:) :: Py

  real(KINDR), allocatable, dimension(:,:) :: TBV, CTE
  real(KINDR), allocatable, dimension(:) :: Amat, means, phenotypes, theta, fixeff
  real(KINDR), allocatable, dimension(:) :: chiasmaCumP, mutationCumP
  real(KINDR), allocatable, dimension(:) :: farmIndReal, farmEffects
  real(KINDR) :: mutationCumP0, chiasmaCumP0
  real(KINDR), allocatable, dimension(:) :: weight
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
  real(KINDR) :: val1, val2, val3, val4
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

  call readInput(inputfile, verbose, nchr, filename1, filename2,&
       chrL, mutationRate, nQTL, nSNP, randomQTL, MAF, baseNameFreq, ncomp, &
       vars, nanim, n_m, n_fpm, n_opf, interval, nlox, nFarm, farmRange, &
       allocation, means, nobs, selectionType, weight, ngen, VarEst, &
       reactionNorm, analysisType, nfix, nvar, nran, outputfile)

  ! the rest of allocations
  allocate(farmBounds(nfarm, 2))
  farmBounds(1:nfarm, 2) = ZERO
  allocate(locations(nanim, nlox))
  locations(1:nanim, 1:nlox) = ZERO
  allocate(X(nobs, nfix))
  X(1:nobs, 1:nfix) = ZERO
  allocate(farmInd(nobs))
  farmInd(1:nobs) = 0
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
  allocate(SNPlist(nChr, nSNP))
  allocate(QTLlist%indices(nChr, nQTL), QTLlist%values(nChr, nQTL, nComp))
  allocate(TBV(nanim, ncomp))
  i = nAnim * (nAnim + 1) / 2 
  allocate(AMat(i))
  allocate(fixeff(nfix))
  allocate(raneff(nran))
  maxid = nanim ! maxval(indiv)
  allocate(raneff(1)%level(maxid)) ! slope effect (genetic)
  if (nran == 3) then
     allocate(raneff(2)%level(maxid)) ! intercept effect (genetic)
     allocate(raneff(3)%level(nobs))   ! environment slope effect (diagonal)
  elseif (nran == 1) then
  else
     write(STDERR, *) " ERROR"
     write(STDERR, *) " not implemented for nran != 1 or 3"
     stop 2
  end if
  allocate(chr_nlocibefore(nchr))
  allocate(ipiv(nobs),Py(nobs))
  nelement = nObs * (nObs + 1) / 2 ! size V matrix, P, etc.
  allocate(P(nelement), V(nelement))
  allocate(Vhat(nfix, nobs))
  allocate(temp(nobs, nfix))
  ! for covariate analysis, 2 step RN is required
  allocate(FarmIndReal(nobs), farmEffects(nfarm))
  allocate(theta(nvar+1))
     
  call defineFarms(interval, nfarm, farmRange, farmBounds)
  open(1, file = 'farms.txt')
  do i = 1, nfarm
     write(1, *) farmBounds(i, 1), farmBounds(i, 2)
  end do
  close(1)

  open(newUnit = iunoutput, file = outputfile)
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
  ! saving QTL list
  !open(1, file = "QTLlist.txt")
  !do ichr = 1, nchr
  !   do i = 1, nQTL
  !      write(1, *) QTLlist%values(ichr, i, 1:2)
  !   end do
  !end do
  !close(1)
  !! saving SNP list
  !open(1, file = "SNPlist.txt")
  !do ichr = 1, nchr
  !   do i = 1, nSNP
  !      write(1, '(2i8)') ichr, SNPlist(ichr, i)
  !   end do
  !end do
  !close(1)
  ! ================================================================
  ! Simulating TBV
  ! ================================================================
  if (verbose) write(STDOUT, '(a,i2)') "simulating BV for generation", 0
  call SimulateTBV(nAnim, nChr, nComp, indiv, genome1, chr_nlocibefore,&
       QTLlist, TBV, verbose)
  if (verbose) write(STDOUT, '(a)') "breeding values simulated"
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
  !! saving first TBV and the QTLlist
  !open(1,file = 'tbvinitial')
  !do i = 1, nanim
  !   write(1, *) (tbv(i,j), j = 1, ncomp)
  !end do
  !close(1)
  !open(1, file = 'qtllist')
  !do i = 1, nchr
  !  do j = 1, nqtl
  !      write(1, *) i,qtlList%indices(i,j),(qtlList%values(i,j,k), &
  !           k = 1, ncomp)
  !   end do
  !end do
  !close(1)

  ! ================================================================
  ! Simulating Phenotypes
  ! ================================================================
  if (verbose) write(STDOUT, '(a)') "allocating individuals"
  i = allocation
  allocation = 1 ! temporary allocation is set to 1 because in the first
  ! generation we do not know the pedigree therefore we can do only random
  call allocateInd(nAnim, nlox, nobs, nfarm, allocation, farmBounds, &
       farmRange, farmInd, locations)
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
  call getGen0Variance(nvar, nran, nanim, nobs, interval, &
       vars, phenotypes, theta)
  
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
        call random_number(raneff(1)%level(1:nanim))
        ! selectByIntercept needs raneff of size 1 (or size 2 but then effects must be on the second)
        write(iunoutput, '(2(1x,a15),2(1x,f15.7))') "NaN", "NaN", &
             correlation(raneff(1)%level, TBV(:,1), nanim),&
             correlation(raneff(1)%level, TBV(:,2), nanim)
     case (2,3) ! requires analysis
        RN: if (reactionNorm) then
           if (selectionType.eq.3) then
              ! converting to double as leastSquare takes double
              farmIndReal(1:nobs) = dble(farmInd(1:nobs))
              ! step 1: farm effects
              call leastSquare(verbose, nobs, nfarm, ids, farmIndReal,&
                   phenotypes, farmEffects, ipiv, Py)
              ! scaling farmeffects to [xmin, xmax]
              val1 = minval(farmEffects)
              val2 = maxval(farmEffects)
              farmEffects(1:nfarm) = (val2 - farmEffects(1:nfarm)) /&
                   (val2 - val1)*(interval(2) - interval(1)) + interval(1)
              ! replacing scaled farmeffects with challenge levels
              X(1:nobs, 1) = farmEffects(farmInd(1:nobs))
           elseif (selectionType.eq.2) then
              do i = 1, nobs ! removing the first farm
                 if (farmInd(i) .eq. 1) cycle
                 X(i, farmInd(i) - 1) = ONE
              end do
              !!!!!!!!!!
              !! writing the incidence matrix
              !write(formato, '(a,i2.2)') "incidence", igen - 1
              !open(1, file = trim(formato))
              !333 format(14(f3.1,1x), f3.1)
              !do i = 1, nobs
              !   write(1, 333) X(i,1:nfix) 
              !end do
              !close(1)
           end if
        end if RN
        !! writing all chromosomes
        !do ichr = 1, nchr
        !   write(filename1, '(a,i2.2,a1,i3.3)') "genchr", igen, '.', iChr
        !   open(1, file = trim(filename1))
        !   write(formato, '(a1, i10, a6)' ) '(', parentGenome(ichr)%nblock, 'i12)'
        !   write(1, *) nanim, parentGenome(ichr)%nloci, parentGenome(ichr)%nblock, 0
        !   do id =1 , nanim
        !      do igam = 1, 2
        !         write(1, formato) (parentGenome(ichr)%genotypes(id, igam, i), &
        !              i = 1, parentGenome(ichr)%nblock)
        !      end do
        !   end do
        !   close(1)
        !end do

        ! making Gmatrix
        iscaled = 1 !(0:no, 1:yes)
        ivar  = 1 !(0:sample, 1:2pq, 2:2p'q')
        j = 0 ! imiss (0:mean, 1:ignore)
        i = 1 ! addDom (1:additive, 2:dominance)
        call getGmatrix(nanim, nChr, nSNP, indiv, ParentGenome, SNPlist,&
             chr_nlocibefore, iscaled, ivar, j, i, Amat, verbose)
        !writing AMAT to disk
        !i = maxval(indiv)
        !k = 1
        !do while (i >= 10)
        !   i = i / 10
        !   k = k + 1
        !end do
        !write(formato, '(a4,i2.2,a4)') "AMAT", igen - 1, ".txt"
        !open( 1, FILE = trim(formato), STATUS = 'unknown' )
        !i = nanim * ( nanim + 1 ) / 2
        !write(formato,'(a2,i1,a5,i1,a22)')'(i',k,',1x,i',k,&
        !     ',1x,g24.15,i9,g24.15)'
        !write(6,*) " formato= ", trim(formato)
        !k = 0
        !do i = 1, nanim
        !   do j = 1, i
        !      k = k + 1
        !      write(1, formato) indiv(i), indiv(j), amat(k)
        !   end do
        !end do
        !close(1)
        !=========================================
        varianceEstimation: select case(varEst)
        case(1) ! do reml
           i = 0 ! number of em iterations
           j = 20 ! number of ai iterations
           call reml(ids, X, phenotypes, nfix, nobs, maxid, nelement, Amat,&
                nvar, nran, theta, verbose, ipiv, Py, P, V, Vhat, temp, i, j)
           ! reml may fail however
           i = nvar +1
           if ( any(theta(1:2)<0) .or. ( (nran==3) .and. &
                ( (theta(3)<0) .or. (theta(i)<0) ) ) ) then
              write(STDERR, '(a)') "Error:"
              write(STDERR, *) "REML failed to obtain variance components"
              stop 2
           end if
        case(2) ! use from generation 0
           call getGen0Variance(nvar, nran, nanim, nobs, interval, vars, &
                phenotypes, theta)
        case(3) ! use true values
           call getTrueVariance(nvar, nran, nanim, ncomp, nobs, interval, &
                vars, tbv, phenotypes, theta)
        end select varianceEstimation
        
        ! now calculating effects
        call blup(ids, X, phenotypes, nfix, nobs, maxid, nelement, Amat,&
             nvar, nran, theta, fixeff, raneff, verbose, ipiv, Py, P, V, &
             temp, Vhat)

        !! write ebvs
        !write(formato, '(a3,i2.2)') 'ebv', igen
        !open(1, file = trim(formato))
        !k = 1
        !if (nran == 3) k = k + 1
        !do i = 1, nanim
        !   write(1, *) (raneff(j)%level(i), j = 1, k)
        !end do
        !close(1)
        if (verbose) write(STDOUT, '(a)') "Genomic evaluation done"
        
        ! -----------------------------------
        ! getting ebv accuracy
        ! ----------------------------------
        if (selectionType.eq.3) then
           write(iunoutput, '(2(1x,f15.7))', advance = 'no') fixEff(1:nfix)
        elseif (selectionType.eq.2) then
           write(iunoutput, '(a15,1x,f15.7)', advance= 'no') "NaN", fixEff(nfix)
        end if
        do i = 1, 2
           j = i
           if ((i == 2) .and. (selectionType == 2)) j = 1
           val1 = correlation(raneff(j)%level, TBV(:,i), nanim)
           write(6, *) 'accuracy', i, val1
           write(iunoutput, '(1x,f15.7)', advance = 'no') val1
        end do
        write(iunoutput, *)
        !====================================
        ! selection
        !====================================
        ! i is used in selectMates, (new implementation ==> i = 1 always)
        if (selectionType.eq.3) then
           ranEff(1)%level(1:nanim) = weight(1) * ranEff(1)%level(1:nanim) +&
                weight(2) * ranEff(2)%level(1:nanim)
        end if
     end select genSel
     i = 1 ! new implementation i = 1 always, (todo: use Py instead)
     call selectParents(nanim, indiv, sex, n_m, n_fpm, male, female,&
          ranEff(i)%level)
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
     call SimulateTBV(nAnim, nChr, nComp, indiv, ParentGenome, &
          chr_nlocibefore, QTLlist, TBV, verbose)
     if (verbose) write(STDOUT, '(a)') "breeding values simulated"

     ! shifting true breeding values
     do i = 1, nComp
        TBV(1:nAnim, i) = TBV(1:nAnim, i) - CTE(i, 2)
     end do
     !! writing tbv
     !write(filename1, '(a,i2.2)') "TBV.G", igen
     !open(1, file = trim(filename1))
     !do i = 1, nanim
     !   write(1, *) (tbv(i,j), j = 1, ncomp)
     !end do
     !close(1)

     if (verbose) write(STDOUT, '(a)') "allocating individuals"
     call allocateInd(nAnim, nlox, nobs,nfarm, allocation, farmBounds, &
          farmRange, farmInd, locations, pedigree = pedigree, nm = n_m,&
          male = male)
     if (verbose) write(STDOUT, '(a)') "individuals allocated"

     ! getting phenotypes
     if (verbose) write(STDOUT, '(a)') "simulating phenotypes"
     call SimulatePhenotype(verbose, nAnim, nComp, nFix, nLox, nran, indiv,&
          TBV, vars, means, nobs, locations, ids, phenotypes, "cov", X)
     if (verbose) write(STDOUT, '(a)') "Phenotypes simulated"

     !! phenotyp file
     !write(filename1, '(a,i2.2)') "phen", igen
     !open(1, file = trim(filename1))
     !do i = 1, nobs
     !   write(1, *) ids(i), X(i, 1), farmind(i), phenotypes(i)
     !end do
     !close(1)

     val1 = sum(tbv(1:nanim,1))/nanim
     val2 = sum(tbv(1:nanim,2))/nanim
     val3 = variance(tbv(1:nanim,1), nanim)
     val4 = variance(tbv(1:nanim,2), nanim)
     !write(6, 68) val1, val2

68   format("slope: ", g25.14, "; intercept: ", g25.14)

     write(iunoutput, '(i3,4(1x,f15.7))', advance = 'no') igen, &
          val1, val2, val3, val4

  end do gen
  write(iunoutput, '(4(1x,a15))') "NaN",  "NaN",  "NaN",  "NaN"
  close(iunoutput)
  call ifinal(seed, startfile)
end program selection
