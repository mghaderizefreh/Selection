subroutine readInput(inputfile, verbose, nanim, nchr, genepoolfile, &
     geneposfile, chrL, mu, ncomp, vars, nQTL, nSNP, locations, selectionType,&
     nobs, means, analysisType, theta, n_m, n_fpm, n_opf, ngen, doreml, nfix, nvar)
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
  real(KINDR), allocatable, dimension(:,:), intent(out) :: locations
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
  integer, intent(out) :: nfix
  integer, intent(out) :: nvar

  character(len = 200) :: line, formato
  integer :: iinput, stat, iun, lno, j, i, nlox
  real(KINDR) :: rinput
  logical :: bool, locRandom, bool2
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
  allocate(vars%corr(ncomp,ncomp))
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
        write(formato, '(a10,2(i1,a1))') "vars%corr(", i, ",", j, ")"
        write(STDOUT, 35) trim(formato), vars%corr(i,j)
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

  ! nlox
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read number of locations per indiv", lno)
  call assert((iinput.gt.0).and.(iinput.lt.nanim)&
       , "number of locations must be > 0 and < nanim", lno)
  nlox = iinput
  write(STDOUT, 34) "number of locations per individual", nlox

  ! random locations  
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to read randomness of locations", lno)
  locRandom = iinput == 1
  write(STDOUT, 33) "random location for individuals?", locRandom
  ! if locations are random
  if (locRandom) then
     allocate(locations(nanim, nlox))
     call random_number(locations)
  else ! if not, further information is required
     call nextInput(iun, line, lno)
     read(line, *, iostat = stat) iinput
     call assert(stat.eq.0, "failed to read wether locations are common", lno)
     bool = iinput == 1
     write(STDOUT, 33) "common locations for all individuals?", bool
     if (bool) then! if locations are common, j locations are required
        allocate(locations(1,nlox))
        do i = 1, nlox
           call nextInput(iun, line, lno)
           read(line, *, iostat = stat) rinput
           call assert(stat.eq.0, "failed to read one of the locations", lno)
           locations(1,i) = rinput
           write(formato, '(a15,1x,i2,1x,a2)') "common location", i, "is"
           write(STDOUT, 35) trim(formato), locations(1,i)
        end do
     else ! otherwise a filename is required
        allocate(locations(nanim, nlox))
        call nextInput(iun, line, lno)
        inquire(file = trim(line), exist = bool2)
        call assert(bool2, "file for locations does not exist", lno)
        open(1, file = trim(line))
        do i = 1, nanim
           read(1, *, err = 101) locations(i, 1:j)
        end do
        write(STDOUT, 36) trim(line), " was successfully read"
     end if
  end if

  ! selection type
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, &
       "selection type must be an integer from 1 to 5 incl.", lno)
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
     if (Locrandom) then
        call assert(.false., &
             "cannot do single trait analysis when locations are random", lno)
     elseif (.not.bool) then
        call assert(.false.,&
             "cannot do single trait analysis when locations are not common", lno)
     end if
  elseif (analysisType .eq. 1) then
     if ((.not.LocRandom).and.bool.and.(size(locations,2).eq.1)) then
        call assert(.false., &
             "RR analysis not possible when animals  phenotyped at 1 location",&
             lno)
     end if
  end if

  ! reml
  call nextInput(iun, line, lno)
  read(line, *, iostat = stat) iinput
  call assert(stat.eq.0, "failed to understand wether to do reml", lno)
  doreml = iinput == 1
  write(STDOUT, 33) "reml required?" , doreml

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

  ! making theta
  if ((.not.locRandom).and.bool.and.(size(locations,2).eq.1)) then
     allocate(theta(2))
     nvar = 1
     nfix = 1
     theta(1) = vars%A(1)* locations(1,1)**2 + vars%A(2) + 2* vars%corr(1,2)*&
          sqrt(vars%A(1) * vars%A(2))
     theta(2) = theta(1) + vars%E(1)* locations(1,1)**2 + vars%E(2) 
  else
     nfix = 2
     if (vars%corr(1,2) .eq. ZERO) then
        allocate(theta(4))
        nvar = 3
     else
        allocate(theta(5))
        nvar = 4
        theta(4) = vars%corr(1,2) * sqrt(vars%A(1) * vars%A(2))
     end if
     theta(1:2) = vars%A(1:2)
     theta(3) = vars%E(1)
     theta(nvar+1) = vars%E(2)
  end if
  write(STDOUT, 34) "number of fix effects", nfix
  write(STDOUT, 34) "number of variance components", nvar
  do i = 1, nvar + 1
     write(formato, '(a6, i1,a1)') "theta(", i,")"
     write(STDOUT, 35) trim(formato), theta(i)
  end do

  write(STDOUT, '(a,i4,a)') "all inputs read in", lno, " lines"
  close(iun)
  return
101 call assert(.false.,&
       " location file must be in format nanim x nlox real values", lno)
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
  write(STDERR, '(a21,i2)') "ERROR in line number ", no
  write(STDERR, *) msg
  stop 2
end subroutine assert
