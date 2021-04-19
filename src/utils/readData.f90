subroutine readInput(inputfile, verbose, nchr, genepoolfile, geneposfile,&
     chrL, mu, nQTL, nSNP, randomQTL, MAF, baseNameFreq, ncomp, vars, nanim,&
     n_m, n_fpm, n_opf, interval, nlox, nFarm, farmRange, allocation, means,&
     nobs, selectionType, weight, ngen, VarEst, reactionNorm, analysisType,&
     nfix, nvar, nran, output)
  use constants
  implicit none
  character(len=*), intent(in) :: inputfile ! list of all inputs
  logical, intent(out) :: verbose
  !!!!!!!!!!!!!!!!!!!! genomic !!!!!!!!!!!!!!!!!!!!
  integer, intent(out) :: nchr
  character(len=100), intent(out) :: genepoolfile
  character(len=100), intent(out) :: geneposfile
  real(KINDR), intent(out) :: chrL
  real(KINDR), intent(out) :: mu
  integer, intent(out) :: nQTL
  integer, intent(out) :: nSNP
  logical, intent(out) :: randomQTL
  real(KINDR), intent(out) :: MAF
  character(len=100), intent(out) :: baseNameFreq
  !!!!!!!!!!!!!!!!!!!! genetic !!!!!!!!!!!!!!!!!!!!
  integer, intent(out) :: nComp
  type(variances), intent(out) :: vars
  !!!!!!!!!!!!!!!!!!!! population !!!!!!!!!!!!!!!!!!!!
  integer, intent(out) :: nanim
  integer, intent(out) :: n_m ! number of males
  integer, intent(out) :: n_fpm ! number of females per male
  integer, intent(out) :: n_opf ! number of offsprings per female
  !!!!!!!!!!!!!!!!!!!! phenotypic !!!!!!!!!!!!!!!!!!!!
  real(KINDR), dimension(2), intent(out) :: interval
  integer, intent(out) :: nlox
  integer, intent(out) :: nFarm
  real(KINDR), intent(out) :: farmRange
  integer, intent(out) :: allocation ! allocation scenario 1,random;2,clust
  real(KINDR), allocatable, dimension(:), intent(out) :: means
  integer, intent(out) :: nobs
  !!!!!!!!!!!!!!!!!!!! selection !!!!!!!!!!!!!!!!!!!!
  integer, intent(out) :: selectionType !1 = random, 2 = local, 3 = index
  real(KINDR), allocatable, dimension(:), intent(out) :: weight
  integer, intent(out) :: ngen ! number of generations
  !!!!!!!!!!!!!!!!!!!! analysis !!!!!!!!!!!!!!!!!!!!
  integer, intent(out) :: VarEst !how vars are est.(1=reml,2=gen0,3=true)
  logical, intent(out) :: reactionNorm ! whether to estimate farm effects
  integer, intent(out) :: analysisType !
  integer, intent(out) :: nfix ! number of fixed effects
  integer, intent(out) :: nvar ! number of variance components
  integer, intent(out) :: nran ! number of random effects
  !!!!!!!!!!!!!!!!!!!! output !!!!!!!!!!!!!!!!!!!!
  character(len=100), intent(out) :: output
  !!!!!!!!!!!!!!!!!!!! dummy variables !!!!!!!!!!!!!!!!!!!!
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

  !!!!!!!!!!!!!!!!!!!! genomic !!!!!!!!!!!!!!!!!!!!
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

  !!!!!!!!!!!!!!!!!!!! genetic !!!!!!!!!!!!!!!!!!!!
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

  !!!!!!!!!!!!!!!!!!!! population !!!!!!!!!!!!!!!!!!!!
  ! nanim
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read input for nanim", lno)
  nanim = iinput
  write(STDOUT, 34) "number of animals/gen", nanim  

  ! number of males (n_m)
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read n_m", lno)
  n_m = iinput
  write(STDOUT, 34) "number of males", n_m

  ! number of females per males (n_fpm)
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read n_fpm", lno)
  n_fpm = iinput
  write(STDOUT, 34) "number of females per male", n_fpm
  call nextInput(iun, line, lno)

  ! number of offspring fer female (n_opf)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read n_opf", lno)
  n_opf = iinput
  write(STDOUT, 34) "number of offspring per female", n_opf

  ! sanity check
  call assert(n_m*n_fpm*n_opf.eq.nanim, &
       "n_m x n_fpm x n_opf != nanim", lno)

  !!!!!!!!!!!!!!!!!!!! phenotypic !!!!!!!!!!!!!!!!!!!!
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

  ! nlox (number of locations or records per individual)
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of locations per indiv", lno)
  call assert((iinput.gt.0).and.(iinput.lt.nanim)&
       , "number of locations must be > 0 and < nanim", lno)
  nlox = iinput
  write(STDOUT, 34) "number of locations per individual", nlox

  ! number of farms
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of farms", lno)
  nFarm = iinput
  write(STDOUT, 34) "number of farms", nFarm

  ! farm range (input as proportion)
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) rinput
  call assert(stat.eq.0, "failed to read range of each farm", lno)
  call assert((rinput.gt.ZERO).and.(rinput.lt.ONE), &
       "farm range must be betwen 0.0 and 1.0 (excl.)", lno)
  farmRange = rinput * (interval(2) - interval(1))
  write(STDOUT, 35) "range of each farm", farmRange

  ! allocation scenario
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read allocation scenario", lno)
  call assert((iinput.eq.1).or.(iinput.eq.2), "so far only random and &
       &clustered allocation are allowed", lno)
  allocation = iinput
  write(STDOUT, 34, advance = 'no') "allocation scenario", allocation
  select case(allocation)
  case(1)
     write(STDOUT, *) "(random)"
  case(2)
     write(STDOUT, *) "(clustered)"
  case default
     write(STDOUT, *) "UNSPECIFIED!!"
  end select

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
  
  ! nobs
  nobs = nlox * nanim
  write(STDOUT, 34) "number of records", nobs

  !!!!!!!!!!!!!!!!!!!! selection !!!!!!!!!!!!!!!!!!!!
  ! selection type
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read selection type", lno)
  call assert((iinput.gt.0).and.(iinput.le.3), &
       "selection type must be an integer from 1 to 3 incl.", lno)
  selectionType = iinput
  write(STDOUT, 34, advance = 'no') "selection type", selectionType
  select case (selectionType)
  case(1) 
     write(STDOUT, *) "(random)"
  case(2) 
     write(STDOUT, *) "(overall performance)"
  case(3)
     write(STDOUT, *) "(manual index)"
  case(4)
     write(STDOUT, *) "(UNSPECIFIED!!)"
  end select
     ! weight
  allocate(weight(ncomp))
  do i = 1, ncomp
     call nextInput(iun, line, lno) ! for slope
     read(line, *, iostat = stat) rinput
     call assert(stat.eq.0, "failed to read weight", lno)
     weight(i) = rinput
     write(formato, '(a, i1)') "weight for comp ", i
     if (selectionType.ne.3) &
          write(formato, '(a,a)') trim(formato), "(is ignored)"
     write(STDOUT, 35) trim(formato), weight(i)
  end do

  ! number of generations
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of generations", lno)
  ngen = iinput
  write(STDOUT, 34) "number of generations", ngen

  !!!!!!!!!!!!!!!!!!!! analysis !!!!!!!!!!!!!!!!!!!!
  ! varEst
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to understand how to compute variances",&
       lno)
  call assert(((iinput.ge.1).and.(iinput.le.3)), &
       "varEst can be 1, 2, 3", lno)
  varEst = iinput
  write(STDOUT, 34, advance = 'no') "variance estimate", varEst
  if (varEst.eq.1) write(STDOUT, *) "(reml)"
  if (varEst.eq.2) write(STDOUT, *) "(generation 0)"
  if (varEst.eq.3) write(STDOUT, *) "(true value)"

  ! reactionNorm
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to understand wether to estimate farm &
       &effects", lno)
  reactionNorm = iinput == 1
  if (selectionType .eq. 1) then
     write(STDOUT, 33) "estimating farm eff (is ignored)?", reactionNorm
  elseif (selectionType.eq.2) then
     write(STDOUT, 33) "estimating farm effects (ST)?", reactionNorm
  elseif (selectionType .eq.3) then
     write(STDOUT, 33) "estimating farm effects (RN)?", reactionNorm
  end if

  ! analysis type
  if ((selectionType.eq.1) .or. (selectionType.eq.2)) then
     analysisType = 2 ! single trait
  elseif (selectionType .eq. 3) then
     analysisType = 1 ! covariate
  end if
  write(STDOUT, 34, advance = 'no') "analysis type", analysisType
  if (analysisType .eq. 1) write(STDOUT, *) "(random regression)"
  if (analysisType .eq. 2) write(STDOUT, *) "(single trait)"

  ! setting number of fixed and random effects and variance compoentns
  if (analysisType .eq. 2) then
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
        nvar = 3 ! var_Ai, var_As, var_Es, (excl. var_Ei)
     else
        nvar = 4 ! var_Ai, var_As, var_Es, cov, (excl. var_Ei)
     end if
  else
  end if
  write(STDOUT, 34) "number of fix effects", nfix
  write(STDOUT, 34) "number of variance components", nvar
  write(STDOUT, 34) "number of random effects", nran

  !!!!!!!!!!!!!!!!!!!! output !!!!!!!!!!!!!!!!!!!!
  ! outputfile
  call nextInput(iun, line, lno)
  output = trim(line)
  write(STDOUT, 36) "output file", trim(output)

  !!!!!!!!!!!!!!!!!!!! FINISHED !!!!!!!!!!!!!!!!!!!!
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
