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
  character(len=100) :: baseNameFreq, outputfile, logfile
  integer :: iunoutput, iunlog, iunFarm, iuniter
  character(len=60) :: formato
  type(chromosome), dimension(:), allocatable, target :: genome1, genome2
  type(chromosome), dimension(:), pointer :: Parentgenome, Offgenome, thisGenome
  integer, dimension(:), allocatable :: seed
  real(KINDR) :: chrL, mutationRate
  integer, dimension(:), allocatable :: indiv, male, female
  logical, dimension(:), allocatable :: sex !true = male, false = female
  integer, dimension(:,:,:), pointer :: parentGen, offGen!, this

  logical :: randomQTL, reactionNorm
  integer :: varEst
  real(KINDR), allocatable, dimension(:) :: P, V, values
  real(KINDR), allocatable, dimension(:,:) :: Vhat, temp
  real(KINDR), dimension(2) :: interval, chalvl
  real(KINDR) :: farmRange, maf
  real(KINDR), dimension(:,:), allocatable :: locations, farmBounds, X !incid. mat.
  integer, dimension(:), allocatable :: farmInd
  integer :: nlox, nfarm, allocation

  !ncomp is the number of traits by a QTL (slope, intercept --> 2, otherwise 1)
  integer :: nSNP, nQTL, nComp
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array) :: QTLlist

  integer, allocatable, dimension(:) :: chr_nlocibefore, ipiv
  real(KINDR), allocatable, dimension(:) :: Py

  real(KINDR), allocatable, dimension(:,:) :: tbv, CTE
  real(KINDR), allocatable, dimension(:) :: Amat, means, phenotypes,&
       theta, fixeff, old_theta, farmIndReal, farmEffects, weight
  real(KINDR), allocatable, dimension(:) :: chiasmaCumP, mutationCumP
  real(KINDR) :: mutationCumP0, chiasmaCumP0
  type(Jarr), dimension(:), allocatable :: raneff
  integer, allocatable, dimension(:) :: totalChiasma, totalMutation, ids
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
     write(STDERR, 100) "reading/setting seed faild"
     stop 2
  end if
  call random_seed(put = seed)

  status = "old"
  call askFileName(inputfile, "input file: ", status, estatus)
  if (status(1:1) .eq. "x") then
     write(STDERR, 101) "error in openning input file ", trim(inputfile)
     stop 2
  end if

  call readInput(inputfile, verbose, nchr, filename1, filename2,&
       chrL, mutationRate, nQTL, nSNP, randomQTL, MAF, baseNameFreq, ncomp, &
       vars, nanim, n_m, n_fpm, n_opf, interval, nlox, nFarm, farmRange, &
       allocation, means, nobs, selectionType, weight, ngen, VarEst, &
       reactionNorm, analysisType, chalvl, nfix, nvar, nran, outputfile,&
       logfile)

  ! the rest of allocations
  call alloc2D(farmBounds, nfarm, 2, "farmBounds", "main")
  farmBounds(1:nfarm, 2) = ZERO
  call alloc2D(locations, nanim, nlox, "locations", "main")
  locations(1:nanim, 1:nlox) = ZERO
  call alloc2D(X, nobs, nfix, "X", "main")
  X(1:nobs, 1:nfix) = ZERO
  call alloc1I(farmInd, nobs, "farmInd", "main")
  farmInd(1:nobs) = 0
  call alloc1I(indiv, nanim, "indiv", "main")
  call alloc1L(sex, nanim, "sex", "main")
  call alloc2D(CTE, nComp, 2, "CTE", "main")
  indiv= (/( i, i = 1, nanim )/)
  ! the array sex corresponds to the new generation (before selection)
  sex(1:nanim) = .false. ! all female
  sex(1:nanim:n_opf) = .true. ! for each female one offspring is male
  call alloc2I(pedigree, nanim, 3, "pedigree", "main")
  pedigree(1:nanim, 1:3) = 0
  i = n_m * n_fpm
  call alloc1I(male, n_m, "male", "main")
  call alloc1I(female, i, "female", "main")
  call alloc1I(ids, nobs, "ids", "main")
  call alloc1D(phenotypes, nobs, "phenotypes", "main")
  allocate(genome1(nChr))
  allocate(genome2(nChr))
  call alloc2I(SNPlist, nChr, nSNP, "SNPlist", "main")
  call alloc2I(QTLlist%indices, nChr, nQTL, "QTLlist%indices", "main")
  call alloc3D(QTLlist%values, nChr, nQTL, nComp, "QTLlist%values", "main")
  call alloc2D(tbv, nanim, ncomp, "tbv", "main")
  i = nAnim * (nAnim + 1) / 2 
  call alloc1D(AMat, i, "AMat", "main")
  call alloc1D(fixeff, nfix, "fixEff", "main")
  allocate(ranEff(nRan))
  call alloc1D(values, nanim, "values", "main")
  maxid = nanim ! maxval(indiv)
  call alloc1D(ranEff(1)%array, maxid, "ranEff(1)%array", "main") ! slope effect (genetic)
  if (nran == 3) then
     call alloc1D(raneff(2)%array,maxid,"ranEff(2)%array","main")! intercept effect (genetic)
     call alloc1D(raneff(3)%array,nobs,"ranEff(3)%array","main")! env slope effect (diagonal)
  elseif (nran == 1) then
  else
     write(STDERR, *) " ERROR"
     write(STDERR, *) " not implemented for nran != 1 or 3"
     stop 2
  end if
  call alloc1I(chr_nlocibefore, nchr, "chr_nlociBefore", "main")
  call alloc1I(ipiv, nobs, "ipiv", "main")
  call alloc1D(Py, nobs, "Py", "main")
  nelement = nObs * (nObs + 1) / 2 ! size V matrix, P, etc.
  call alloc1D(p, nelement, "p", "main")
  call alloc1D(V, nelement, "v", "main")
  call alloc2D(Vhat, nfix, nobs, "Vhat", "main")
  call alloc2D(temp, nobs, nfix, "temp", "main")
  ! for covariate analysis, 2 step RN is required
  call alloc1D(FarmIndReal, nobs, "farmIndReal", "main")
  call alloc1D(farmEffects, nfarm, "farmEffects", "main")
  i = nvar + 1
  call alloc1D(theta, i, "theta", "main")
  call alloc1D(old_theta, i, "old_theta", "main")
     
  call defineFarms(interval, nfarm, farmRange, farmBounds)
  open(newUnit = iunFarm, file = 'farms.txt')
  do i = 1, nfarm
     write(iunFarm, *) farmBounds(i, 1), farmBounds(i, 2)
  end do
  close(iunFarm)

  open(newUnit = iunoutput, file = outputfile)
  open(newUnit = iunlog, file = logfile)
  open(newUnit = iunIter, file = "iteration.var")
  write(iunlog, 204, advance = 'no') "g", "anim", "sire", "dam",&
       "nlox", "tbv_s", "tbv_i", "ebv_s", "ebv_i", "value"
  do i = 1, nlox 
     write(iunlog, 205, advance = 'no') i, i, i, i
  end do
  write(iunlog, *) 
  write(iunIter, '(a)') "iteration estimation"
  ! ==============================================
  ! simulating first individuals from genepool
  ! ==============================================
  i = 3 ! genstart = 3
  j = 1 ! istore = 1
  k = 1 ! withpos = 1
  call initialiseGenotypes(verbose, nchr, nanim, i, nloci, nblock, j, &
       genome1, genome2, maxloci, maxblock, ifail, chrL, &
       trim(filename1), trim(filename2))
  if (verbose) write(STDOUT, 100) "Genotypes Initialised"
  ! ==============================================
  ! calculating the cumulative distribution for number of recombination events
  ! ==============================================
  i = 200 ! upto 200 recombinations
  if (verbose) write(STDOUT, 100) &
       "Setting total number of mutations and recombinations"
  call GetMutRecArray(verbose, i, chrL, mutationRate, nLoci, &
       chiasmaCumP, chiasmacumP0, totalChiasma, maxchiasma, &
       mutationCumP, mutationCumP0, totalMutation, maxmutations)
  if (verbose) write(STDOUT, 100) "Mutation and recombinations set"
  ! ================================================================
  ! Getting QTL and SNP list
  ! ================================================================
  if (verbose) write(STDOUT, 100) "Getting SNPs and their values"
  call getQTLandSNP(verbose, nChr, nQTL, nSNP, nComp, randomQTL, genome1,&
       QTLlist, SNPlist, vars%cov, baseNameFreq, maf)
  if (verbose) write(STDOUT, 100) "QTL and SNP list simulated"
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
  if (verbose) write(STDOUT, 206) "simulating BV for generation", 0
  call SimulateTBV(nAnim, nChr, nComp, indiv, genome1, chr_nlocibefore,&
       QTLlist, TBV, verbose)
  if (verbose) write(STDOUT, 100) "breeding values simulated"
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
  if (verbose) write(STDOUT, 100) "allocating individuals"
  i = allocation
  allocation = 1 ! temporary allocation is set to 1 because in the first
  ! generation we do not know the pedigree therefore we can do only random
  call allocateInd(nAnim, nlox, nobs, nfarm, allocation, farmBounds, &
       farmRange, farmInd, locations)
  allocation = i
  if (verbose) write(STDOUT, 100) "individuals allocated"

  if (verbose) write(STDOUT, 100) "simulating phenotypes"
  call SimulatePhenotype(verbose, nAnim, nComp, nFix, nLox, nran, indiv,&
       TBV, vars, means, nobs, locations, ids, phenotypes, "cov", X)
  if (verbose) write(STDOUT, 100) "Phenotypes simulated"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !write(formato, '(a,i2.2)') "phen", 0
  !open(1, file = trim(formato))
  !do i = 1, nobs
  !   write(1, *) ids(i), X(i, 1), farmInd(i), phenotypes(i)
  !end do
  !close(1)

  ! only for the first generation
  ! a very simple guestimate for theta
  call getGen0Variance(nvar, nran, nanim, nobs, interval, chalvl,&
       vars, phenotypes, theta)
  old_theta(:) = theta(:)
  
  ParentGenome => genome1
  Offgenome    => genome2

  if (verbose) write(STDOUT, 206) "generation ", 0 ! must be zero
  val1 = sum(tbv(1:nanim, 1)) / nanim
  val2 = sum(tbv(1:nanim, 2)) / nanim
  write(iunoutput, 200) "gen", "mean_TBVs", "mean_TBVi", "var_TBVs", &
       "var_TBVi", "mu_s", "mu_i", "acc_EBVs", "acc_EBVi",&
       "acc_scoreS", "acc_scoreI", "A_s", "A_i","rho", "E_s", "E_i"
100 format(a)
101 format(2a)
200 format(a3,15(1x,a15))
201 format(i0,4(1x,f15.7))
202 format(i0,1x,3i5,i3,4(f15.7,1x),f15.7)
203 format(1x,i0,3(1x,f15.7))
204 format(11(a,1x))
205 format(" farmInd",i2.2," farmEffectSc",i2.2," env",i2.2, " phen", i2.2)
206 format(a, i2)
207 format(2(1x,a15),4(1x,f15.7))
208 format(2(1x,f15.7))
209 format(1x,a15,1x,f15.7)
210 format(1x,f15.7)
211 format(6(1x,a15))
212 format(5(1x,f15.7))
213 format(2(1x,a15))
214 format(5(1x,a15))
  write(iunoutput, 201, advance = 'no') 0, val1, val2, &
       variance(tbv(1:nanim, 1), nanim), variance(tbv(1:nanim, 2), nanim)
  ! initialising first generation individuals
  gen: do igen = 1, ngen
     write(STDOUT, 206) "generation ", igen
     ! ================================================================
     ! Genomic Evaluation
     ! TODO: This block (genSel) should eventually be one single
     !       subroutine with name "doGenomicSelection"
     ! ================================================================
     genSel: select case (selectionType)
     case (1) ! random
        call random_number(values)
        ! selectByIntercept needs raneff of size 1 (or size 2 but then effects must be on the second)
        val1 = correlation(values, TBV(:,1), nanim)
        val2 = correlation(values, TBV(:,2), nanim)
        val3 = variance(values, nanim)
        val4 = variance(values, nanim)
        write(iunoutput, 207) "NaN", "NaN", val1, val2, val1, val2
        write(iunoutput, 214) "NaN", "NaN", "NaN", "NaN", "NaN" 
        write(iunIter, '(i4,a)') igen, "NA"
     case (2,3) ! requires analysis
        RN: if (reactionNorm) then
           if (selectionType.eq.3) then
              ! converting to double as leastSquare takes double
              farmIndReal(1:nobs) = dble(farmInd(1:nobs))
              ! step 1: farm effects
              call leastSquare(verbose, nobs, nfarm, ids, farmIndReal,&
                   phenotypes, farmEffects, ipiv, Py, ifail)
              if (ifail /= 0) then
                 write(STDERR, *) "stopping"
                 stop 2
              end if
              ! scaling farmeffects to [xmin, xmax]
              val1 = minval(farmEffects)
              val2 = maxval(farmEffects)
              farmEffects(1:nfarm) = (val2 - farmEffects(1:nfarm)) /&
                   (val2 - val1)*(chalvl(2) - chalvl(1)) + chalvl(1)
              ! replacing scaled farmeffects with challenge levels
              X(1:nobs, 1) = farmEffects(farmInd(1:nobs))
           elseif (selectionType.eq.2) then
              do i = 1, nobs ! removing the first farm
                 if (farmInd(i) .eq. 1) cycle
                 X(i, farmInd(i) - 1) = ONE
              end do
              ! and then including a column of ones in x
              X(1:nobs, nFix) = ONE
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

        !!writing AMAT to disk
        !i = maxval(indiv)
        !k = 1
        !do while (i >= 10)
        !i = i / 10
        !k = k + 1
        !end do
        !write(formato, '(a4,i2.2,a4)') "AMAT", igen - 1, ".txt"
        !open( 1, FILE = trim(formato), STATUS = 'unknown' )
        !i = nanim * ( nanim + 1 ) / 2
        !write(formato,'(a2,i1,a5,i1,a22)')'(i',k,',1x,i',k,&
        !',1x,g24.15,i9,g24.15)'
        !write(6,*) " formato= ", trim(formato)
        !k = 0
        !do i = 1, nanim
        !do j = 1, i
        !k = k + 1
        !write(1, formato) indiv(i), indiv(j), amat(k)
        !end do
        !end do
        !close(1)

        !=========================================
        varianceEstimation: select case(varEst)
        case(1) ! do reml
           i = 0 ! number of em iterations
           j = 20 ! number of ai iterations
           call reml(ids, X, phenotypes, nfix, nobs, maxid, nelement, Amat,&
                nvar, nran, theta, verbose, ipiv, Py, P, V, Vhat, temp, &
                ifail, i, j)
           ! if ifail /= 0 then reml failed for some reasons
           ! then try again with em iterations first
           if (ifail /= 0) then
              write(STDOUT, '(A)') " REML failed to find variance components."
              write(STDOUT, '(A)') "   trying again with 8 em iterations first"
              theta(:) = old_theta(:)
              i = 8 ! number of em iterations
              j = 25 ! number of ai iterations
              call reml(ids, X, phenotypes, nfix, nobs, maxid, nelement, Amat,&
                   nvar, nran, theta, verbose, ipiv, Py, P, V, Vhat, temp, &
                   ifail, i, j)
              if (ifail /= 0) then
                 write(STDOUT, '(A)') " REML failed again"
                 write(STDOUT, '(A)') "   doing a blup with previous (or 1st) var comps"
                 theta(:) = old_theta(:)
                 i = 900! code for reml failure
              else
                 i = 998! code for reml success at second attempt
              end if
           else
              i = 999 ! code for reml success at first attempt
           end if
        case(2) ! use from generation 0
           call getGen0Variance(nvar, nran, nanim, nobs, interval, chalvl, vars, &
                phenotypes, theta)
           i = 997 ! code for gen0
        case(3) ! use true values
           call getTrueVariance(nvar, nran, nanim, ncomp, nobs, interval, &
                vars, tbv, phenotypes, theta)
           i = 996 ! code for true vars
        end select varianceEstimation
        
        ! now calculating effects
        call blup(ids, X, phenotypes, nfix, nobs, maxid, nelement, Amat,&
             nvar, nran, theta, fixeff, raneff, verbose, ipiv, Py, P, V, &
             Vhat, temp, ifail)

        ! failure may be due to invalid variance components from reml.
        ! so apart from failure, check whether reml ran succesfully
        blupfail: if (ifail /= 0) then
           if ((varEst == 1) .and. (i >= 998)) then
              write(STDOUT, '(a,/,a)') " warning: blup failed.",&
                 "Trying with previous variances"
              theta(:) = old_theta(:)
              call blup(ids, X, phenotypes, nfix, nobs, maxid, nelement,&
                 Amat, nvar, nran, theta, fixeff, raneff, verbose, ipiv,&
                 Py, p, v, vhat, temp, ifail)
              i = 900
              if (ifail /= 0) then
                 write(STDERR, *) "blup failed several times. stopping"
                 stop 2
              end if
           else
              write(STDERR, *) "blup should not have failed. stopping"
              stop 2
           end if
        end if blupfail

        recordingVarEstScen: select case(i)
        case(999)
           write(iunIter,'(i4,1x,a)') igen, "AI20EM00"
        case(998) 
           write(iunIter,'(i4,1x,a)') igen, "AI25EM08"
        case(997)
           write(iunIter,'(i4,1x,a)') igen, "Gen_Zero"
        case(996)
           write(iunIter,'(i4,1x,a)') igen, "True_Val"
        case(900)
           write(iunIter,'(i4,1x,a)') igen, "Previous"
        end select recordingVarEstScen

        ! finally record what was used by keeping it in old_theta
        old_theta(:) = theta(:)

        !! write ebvs
        !write(formato, '(a3,i2.2)') 'ebv', igen
        !open(1, file = trim(formato))
        !k = 1
        !if (nran == 3) k = k + 1
        !do i = 1, nanim
        !   write(1, *) (raneff(j)%array(i), j = 1, k)
        !end do
        !close(1)
        if (verbose) write(STDOUT, 100) "Genomic evaluation done"
        
        !====================================
        ! selection
        !====================================
        if (selectionType.eq.2) then
           values(1:nanim) = raneff(1)%array(1:nanim)
        elseif (selectionType.eq.3) then
           values(1:nanim) = weight(1) * ranEff(1)%array(1:nanim) +&
                weight(2) * ranEff(2)%array(1:nanim)
        end if

        ! -----------------------------------
        ! getting ebv accuracy
        ! ----------------------------------
        if (selectionType.eq.3) then
           write(iunoutput, 208, advance = 'no') fixEff(1:nfix)
        elseif (selectionType.eq.2) then
           write(iunoutput, 209, advance= 'no') "NaN", fixEff(nfix)
        end if
        do i = 1, 2
           j = i
           if ((i == 2) .and. (selectionType == 2)) j = 1
           val1 = correlation(raneff(j)%array, TBV(:,i), nanim)
           write(6, *) 'accuracy', i, val1
           write(iunoutput, 210, advance = 'no') val1
        end do
        do i = 1, 2
           val1 = correlation(values, tbv(:,i), nanim)
           write(iunoutput, 210, advance = 'no') val1
        end do
        ! -----------------------------------
        ! getting variance of ebv
        ! -----------------------------------
        if (selectionType.eq.2) then
           write(iunoutput, 210, advance = 'no') theta(1)
           write(iunoutput, 213, advance = 'no') "NaN", "NaN"
           write(iunoutput, 210, advance = 'no') theta(2)
           write(iunoutput, '(1x,a15)') "NaN"
        elseif (selectionType .eq. 3) then
           write(iunoutput, 212) theta(1), theta(2), theta(4), &
              theta(3), theta(5)
        end if
        !====================================
        ! recording logs
        !====================================
        k = merge(2, 1, selectionType.eq.3)

        do id = 1, nanim
           write(iunlog, 202, advance = 'no') &
                igen, id, (pedigree(j,2), j=2,3), nlox, tbv(id,1), tbv(id,2), &
                raneff(1)%array(id), raneff(k)%array(id), values(id)
           i = (id - 1) * nlox
           do j = ((id -1)*nlox+1), (id*nlox)
              write(iunlog, 203, advance = 'no') &
                   farmInd(j), X(j,1), locations(id,j-i), phenotypes(j)
           end do
           write(iunlog, *) 
        end do
     end select genSel

     call selectParents(nanim, indiv, sex, n_m, n_fpm, male, female, values)
     if (verbose) write(STDOUT, 100) "selection done"
     
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
     if (verbose) write(STDOUT, 206) "simulating BV for generation", igen
     call SimulateTBV(nAnim, nChr, nComp, indiv, ParentGenome, &
          chr_nlocibefore, QTLlist, TBV, verbose)
     if (verbose) write(STDOUT, 100) "breeding values simulated"

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

     if (verbose) write(STDOUT, 100) "allocating individuals"
     call allocateInd(nAnim, nlox, nobs,nfarm, allocation, farmBounds, &
          farmRange, farmInd, locations, pedigree = pedigree, nm = n_m,&
          male = male)
     if (verbose) write(STDOUT, 100) "individuals allocated"

     ! getting phenotypes
     if (verbose) write(STDOUT, 100) "simulating phenotypes"
     call SimulatePhenotype(verbose, nAnim, nComp, nFix, nLox, nran, indiv,&
          TBV, vars, means, nobs, locations, ids, phenotypes, "cov", X)
     if (verbose) write(STDOUT, 100) "Phenotypes simulated"

     !! phenotyp file
     !write(filename1, '(a,i2.2)') "phen", igen
     !open(1, file = trim(filename1))
     !do i = 1, nobs
     !   write(1, *) ids(i), X(i, 1), farmInd(i), phenotypes(i)
     !end do
     !close(1)

     val1 = sum(tbv(1:nanim,1))/nanim
     val2 = sum(tbv(1:nanim,2))/nanim
     val3 = variance(tbv(1:nanim,1), nanim)
     val4 = variance(tbv(1:nanim,2), nanim)

     write(iunoutput, 201, advance = 'no') igen, val1, val2, val3, val4

  end do gen
  write(iunoutput, 211, advance = 'no') "NaN",  "NaN",  "NaN",  "NaN", "NaN", "NaN"
  write(iunoutput, 211) "NaN",  "NaN",  "NaN",  "NaN", "NaN"
  close(iunoutput)
  close(iunlog)
  close(iunIter)
  call ifinal(seed, startfile)
end program selection
