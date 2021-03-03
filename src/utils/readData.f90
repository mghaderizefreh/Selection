subroutine readInput(inputfile, verbose, nanim, nchr, genepoolfile, &
     geneposfile, chrL, mu, ncomp, vars, nQTL, nSNP, MAF, baseNameFreq,&
     randomQTL, interval, locations, X, nlox, nFarm, farmBounds, &
     farmInd, farmRange, allocation, selectionType, nobs, means, &
     analysisType, theta, n_m, n_fpm, n_opf, ngen, doreml, reactionNorm,&
     nfix, nvar, nran, output)
  use constants
  implicit none
  character(len=*), intent(in) :: inputfile ! list of all inputs

  logical, intent(out) :: verbose
  integer, intent(out) :: nanim
  integer, intent(out) :: nchr
  character(len=100), intent(out) :: genepoolfile
  character(len=100), intent(out) :: geneposfile
  real(KINDR), intent(out) :: chrL
  real(KINDR), intent(out) :: mu
  integer, intent(out) :: nComp
  type(variances), intent(out) :: vars
  integer, intent(out) :: nQTL
  integer, intent(out) :: nSNP
  real(KINDR), intent(out) :: MAF
  character(len=100), intent(out) :: baseNameFreq
  logical, intent(out) :: randomQTL
  real(KINDR), dimension(2), intent(out) :: interval
  real(KINDR), allocatable, dimension(:,:), intent(out) :: locations
  real(KINDR), allocatable, dimension(:,:), intent(out) :: X
  integer, intent(out) :: nlox
  integer, intent(out) :: nFarm
  real(KINDR), allocatable, dimension(:,:), intent(out) :: farmBounds
  integer, allocatable, dimension(:), intent(out) :: farmInd
  real(KINDR), intent(out) :: farmRange
  integer, intent(out) :: allocation ! allocation scenario
  ! 1 = random, 
  integer, intent(out) :: selectionType 
  integer, intent(out) :: nobs
  real(KINDR), allocatable, dimension(:), intent(out) :: means
  integer, intent(out) :: analysisType ! 1 = rr, 2 = single location
  real(KINDR), allocatable, dimension(:), intent(out) :: theta ! for analysis
  ! 1 = random, 2 = slopeEBV, 3 = interceptEBV, 4 = slopeTBV, 5 = interceptTBV
  integer, intent(out) :: n_m ! number of males
  integer, intent(out) :: n_fpm ! number of females per male
  integer, intent(out) :: n_opf ! number of offsprings per female
  integer, intent(out) :: ngen ! number of generations
  logical, intent(out) :: doreml ! whether to do a reml or only a blup
  logical, intent(out) :: reactionNorm ! whether to estimate farm effects
  integer, intent(out) :: nfix ! number of fixed effects
  integer, intent(out) :: nvar ! number of variance components
  integer, intent(out) :: nran ! number of random effects
  character(len=100), intent(out) :: output

  character(len = 200) :: line, formato
  integer :: iinput, stat, iun, lno, j, i
  real(KINDR) :: rinput
  logical :: bool
  lno = 0

  open(newUnit = iun, file = trim(inputfile))

33 format(a34,": ", l1)
34 format(a34,": ", i0)
35 format(a34,": ", g0.15)
36 format(a34,": ", a)

  ! verbosity
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for verbose", lno)
  verbose = iinput .eq. 1
  write(STDOUT, 33) "verbose?", verbose

  ! nanim
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for nanim", lno)
  nanim = iinput
  write(STDOUT, 34) "number of animals/gen", nanim  

  ! number of chrosomomes
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for nchr", lno)
  call assert(iinput.gt.0, "nchr must be greater than 0", lno)
  nchr = iinput
  write(STDOUT, 34) "number of chromosomes", nchr

  ! base-name for genepool filename
  call nextInput(iun, line, lno)
  genepoolfile = trim(line)
  do i = 1, nchr
     write(line, '(a,i3.3)') trim(genepoolfile),i
     inquire(file=line, exist = bool)
     call assert(bool, &
          "error in finding a file for genepool with as many as chromosomes", lno)
  end do
  write(STDOUT, 36) "file for genepool", genepoolfile

  ! base-name for position filename
  call nextInput(iun, line, lno)
  geneposfile = trim(line)
  do i = 1, nchr
     write(line, '(a,i3.3)') trim(geneposfile), i
     inquire(file=line, exist = bool)
     call assert(bool, &
          "error in finding a file for genepos with as many as chromosomes", lno)
  end do
  write(STDOUT, 36) "file for genepositions", geneposfile

  ! chrL
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read input for chrL", lno)
  call assert(rinput.gt.ZERO, "Chromosome length must be > 0.0", lno)
  chrL = rinput
  write(STDOUT, 35) "chromosome length", chrL

  ! mu (mutation rate)
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read input for mu", lno)
  call assert(rinput.gt.ZERO, "mutation rate must be > 0.0", lno)
  mu = rinput
  write(STDOUT, 35) "mutation rate", mu

  ! nComp
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of components", lno)
  call assert(iinput.le.2, "number of components cannot be > 2 (atm)", lno)
  ncomp = iinput
  write(STDOUT, 34) "number of components", nComp

  ! variance (A)
  allocate(vars%A(ncomp))
  do i = 1, nComp
     call nextInput(iun, line, lno)
     read(line, *, iostat = stat) rinput
     call assert(stat.eq.0, "failed to read vA", lno)
     call assert(rinput.gt.ZERO, "variance must be > 0.0", lno)
     vars%A(i) = rinput
     write(formato, '(a7,i1,a1)') "vars%A(", i, ")"
     write(STDOUT, 35) trim(formato), vars%A(i)
  end do

  ! variance (E)
  allocate(vars%E(ncomp))
  do i = 1, nComp
     call nextInput(iun, line, lno)
     read(line, *, iostat = stat) rinput
     call assert(stat.eq.0, "failed to read vE", lno)
     call assert(rinput.gt.ZERO, "variance must be > 0.0", lno)
     vars%E(i) = rinput
     write(formato, '(a7,i1,a1)') "vars%E(", i, ")"
     write(STDOUT, 35) trim(formato), vars%E(i)
  end do

  allocate(vars%PE(ncomp))
  vars%PE = ZERO

  ! correlations
  allocate(vars%corr(ncomp,ncomp), vars%cov(ncomp,ncomp))
  do i = 1, ncomp
     do j = 1, i
        call nextInput(iun, line, lno)
        read(line, *, iostat = stat) rinput
        call assert(stat.eq.0, "failed to read a correlation", lno)
        call assert((rinput.ge.-1._KINDR).and.(rinput.le.ONE),&
             "correlation must be in [-1,+1]", lno)
        if (i.eq.j) call assert(rinput.eq.ONE, "self correlation must be 1.0",&
             lno)
        vars%corr(i,j) = rinput
        vars%corr(j,i) = rinput
        vars%cov(i,j) = vars%corr(i,j) * sqrt(vars%A(i)) * sqrt(vars%A(j))
        vars%cov(j,i) = vars%cov(i,j)
        write(formato, '(a10,2(i1,a1))') "vars%corr(", i, ",", j, ")"
        write(STDOUT, 35) trim(formato), vars%corr(i,j)
        write(formato, '(a9,2(i1,a1))') "vars%cov(", i, ",", j, ")"
        write(STDOUT, 35) trim(formato), vars%cov(i,j)
     end do
  end do

  ! nQTL
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of QTLs", lno)
  call assert(iinput.gt.0, "nQTL must be > 0", lno)
  nQTL = iinput
  write(STDOUT, 34) "number of QTL/chromosome", nQTL

  ! nSNP
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of SNPs", lno)
  call assert(iinput.gt.0, "number of SNP cannot be < 0", lno)
  nSNP = iinput
  write(STDOUT, 34) "number of SNP/chromosome", nSNP

  ! randomQTL
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for randomQTL", lno)
  randomQTL = iinput .eq. 1
  write(STDOUT, 33) "Are QTLs random?", randomQTL

  ! baseNameFreq
  call nextInput(iun, line, lno)
  baseNameFreq = trim(line)
  if (.not.randomQTL) then
     do i = 1, nchr
        write(line, '(a,i3.3)') trim(baseNameFreq),i
        inquire(file=line, exist = bool)
        call assert(bool, "error in finding a file for frequenecies with&
             &as many as chromosomes", lno)
     end do
     write(STDOUT, 36) "file for frequnecy", baseNameFreq
  else
     write(STDOUT, 36) "file for frequency (is ignored)", baseNameFreq
  end if

  ! MAF
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read MAF for finding QTLs", lno)
  if (.not.randomQTL) then
     call assert((rinput.ge.ZERO).and.(rinput.lt.HALF), "maf must be &
          &between 0.0 and 0.5 (incl,excl., respecitvely)",lno)
     maf = rinput
     write(STDOUT, 35) "MAF cutoff for QTL", MAF
  else
     maf = rinput
     write(STDOUT, 35) "MAF cutoff for QTL (is ignored)", MAF
  end if

  ! interval (boundaries of X; xmin, xmax)
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read xmin", lno)
  interval(1) = rinput
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read xmamx", lno)
  interval(2) = rinput
  call assert(interval(2)>interval(1), "xmax must be > xmin", lno)
  write(STDOUT, 35) "xmin", interval(1)
  write(STDOUT, 35) "xmax", interval(2)

  ! nlox
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of locations per indiv", lno)
  call assert((iinput.gt.0).and.(iinput.lt.nanim)&
       , "number of locations must be > 0 and < nanim", lno)
  nlox = iinput
  allocate(locations(nanim, nlox))
  write(STDOUT, 34) "number of locations per individual", nlox

  ! number of farms  
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of farms", lno)
  nFarm = iinput
  allocate(farmBounds(nfarm, 2))
  write(STDOUT, 34) "number of farms", nFarm

  ! farm range
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read range of each farm", lno)
  call assert((rinput.gt.ZERO).and.(rinput.lt.ONE), &
       "farm range must be betwen 0.0 and 1.0 (excl.)", lno)
  farmRange = rinput
  write(STDOUT, 35) "range of each farm", farmRange

  ! allocation scenario
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read allocation scenario", lno)
  call assert(iinput.eq.1, "so far only random allocation is allowed", lno)
  allocation = iinput
  write(STDOUT, 34, advance = 'no') "allocation scenario", allocation
  select case(allocation)
  case(1)
     write(STDOUT, *) "(random)"
  end select

  ! selection type
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read selection type", lno)
  call assert((iinput.gt.0).and.(iinput.lt.7), &
       "selection type must be an integer from 1 to 6 incl.", lno)
  selectionType = iinput
  write(STDOUT, 34, advance = 'no') "selection type", selectionType
  select case (selectionType)
  case(1) 
     write(STDOUT, *) "(random)"
  case(2) 
     write(STDOUT, *) "(EBV slope)"
  case(3) 
     write(STDOUT, *) "(EBV intercept)"
  case(4) 
     write(STDOUT, *) "(TBV slope)"
  case(5) 
     write(STDOUT, *) "(TBV intercept)"
  case(6)
     write(STDOUT, *) "(overall performance)"
  end select

  ! nobs
  nobs = size(locations, 2) * nanim
  write(STDOUT, 34) "number of records", nobs

  ! number of m, fpm, and opf
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read n_m", lno)
  n_m = iinput
  write(STDOUT, 34) "number of males", n_m
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read n_fpm", lno)
  n_fpm = iinput
  write(STDOUT, 34) "number of females per male", n_fpm
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read n_opf", lno)
  n_opf = iinput
  write(STDOUT, 34) "number of offspring per female", n_opf
  call assert(n_m*n_fpm*n_opf.eq.nanim, &
       "n_m x n_fpm x n_opf != nanim", lno)

  ! number of generations
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of generations", lno)
  ngen = iinput
  write(STDOUT, 34) "number of generations", ngen

  ! analysis type
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read analysis type", lno)
  call assert(iinput.eq.1.or.iinput.eq.2, &
       "analysis type may be 1 or 2 only", lno)
  analysisType = iinput
  write(STDOUT, 34, advance = 'no') "analysis type", analysisType
  if (analysisType .eq. 1) write(STDOUT, *) "(random regression)"
  if (analysisType .eq. 2) write(STDOUT, *) "(single trait)"
  if (analysisType .eq. 2) then
     call assert((selectionType.eq.1).or.(selectionType.eq.6), &
       "For single trait analysis, only selection type 1(random) and 6 (o&
       &verall performace) is allowed", lno)
  else     
  end if

  ! reml
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to understand wether to do reml", lno)
  doreml = iinput == 1
  write(STDOUT, 33) "reml required?" , doreml

  ! reactionNorm
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to understand wether to estimate farm &
       &effects", lno)
  reactionNorm = iinput == 1
  if ( (selectionType .eq. 1).or.(selectionType .eq. 4).or.&
       (selectionType .eq. 5)) then
     write(STDOUT, 33) "estimating farm eff (is ignored)?", reactionNorm
  elseif ((selectionType.eq.2).or.(selectionType.eq.3)) then
     write(STDOUT, 33) "estimating farm effects (RN)?", reactionNorm
  elseif ((selectionType .eq. 6)) then
     write(STDOUT, 33) "estimating farm effects (ST)?", reactionNorm
  end if

  ! means
  allocate(means(ncomp))
  do i = 1, nComp
     call nextInput(iun, line, lno)
     read(line, *, iostat = stat) rinput
     call assert(stat.eq.0, "failed to read mean for a component", lno)
     means(i) = rinput
     write(formato, '(a6, i1, a1)') "means(",i,")=" 
     write(STDOUT, 35) trim(formato), means(i)
  end do

  ! outputfile
  call nextInput(iun, line, lno)
  output = trim(line)
  write(STDOUT, 36) "output file", trim(output)

  ! making theta
  if (analysisType .eq. 2) then
     allocate(theta(2))
     theta(1:2) = ZERO
     nvar = 1
     nran = 1
     if (reactionNorm) then
        nfix = nfarm ! +1 (for mu) - 1 (for stability)
     else
        nfix = 1 ! mu
     end if
  elseif (analysisType .eq. 1) then
     nran = 3 ! Ai, As, Es
     nfix = 2 ! mu_i, mu_s
     if (vars%corr(1,2) .eq. ZERO) then
        allocate(theta(4))
        nvar = 3 ! var_Ai, var_As, var_Es, (excl. var_Ei)
     else
        allocate(theta(5))
        nvar = 4 ! var_Ai, var_As, var_Es, cov, (excl. var_Ei)
        theta(4) = vars%cov(1,2)
     end if
     theta(1:2) = vars%A(1:2)
     theta(3) = vars%E(1)
     theta(nvar+1) = vars%E(2)
  else
  end if
  write(STDOUT, 34) "number of fix effects", nfix
  write(STDOUT, 34) "number of variance components", nvar
  do i = 1, nvar + 1
     write(formato, '(a6, i1,a1)') "theta(", i,")"
     write(STDOUT, 35) trim(formato), theta(i)
  end do

  ! dealing with X
  allocate(X(nobs, nfix))
  X(1:nobs, 1:nfix) = ZERO

  ! farm ind
  allocate(farmInd(nobs))

  write(STDOUT, '(a,i4,a)') "all inputs read in", lno, " lines"
  close(iun)
end subroutine readInput
subroutine nextInput(iun, output, lno)
  use constants
  implicit none
  character(len = 200), intent(out) :: output
  integer, intent(in) :: iun
  integer, intent(inout) :: lno
  character(len = 200) :: line
  do
     read(iun, '(A)', end = 100) line
     lno = lno + 1
     if(line(1:1) .eq. '!') cycle
     if(trim(line) .eq. '') cycle
     output = trim(line)
     return
100  write(STDERR, *) "ERROR: unexpected end of input file" 
     write(STDERR, *) lno, "lines were read"
     stop 2
  end do
end subroutine nextInput
subroutine assert(iostat, msg, no)
  use constants
  implicit none
  character(len = *), intent(in) :: msg
  logical, intent(in) :: iostat
  integer, intent(in) :: no
  if (iostat) return
  write(STDERR, '(a21,i3)') "ERROR in line ", no
  write(STDERR, *) msg
  stop 2
end subroutine assert
