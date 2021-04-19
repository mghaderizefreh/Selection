! order of matrices, and therefore variances shall be 
!   1     ZsGZs               genetic slope
!   2     ZiGZi               genetic intercept
!   3     ZsZs                environmental slope (stored as diagonal)
!   4     ZiGZs+ZsGZi         genetic slope-intercept covariance
!   5     ZsZs                permanent environment effect slope (SAME AS NUMBER 3: but counted)
!   6     ZiZi                permanent environment effect intercept
!   7     ZsZi+ZiZs           permanent environment effect slope-intercept covariance
! (last)  Identity            environmental intercept (NOT COUNTED IN theZGZ and it is always LAST one)
program analysis
  use constants
  use global_module
  use blup_module
  use reml_module
  implicit none
  !! ================ variable definitions  ================ !!
  character(LEN=256):: phenFile, AmatFile, fixEffFile, ranEffFile, varFile, msg
  character(len=30) :: status, eStatus, formato
  logical :: verbose
  integer :: ifail, doreml
  integer :: i, j, k, maxid, nvar, nobs, nfix, nran, nelement
  integer :: phenFileID, AmatFileID, iunFix, iunRan, iunVar
  integer :: lines, empties, emiteration, maxIter
  integer, dimension(:), allocatable :: id ! real(KINDR) id of animals
  integer, dimension(:), allocatable :: levels
  integer, allocatable, dimension(:) :: ipiv
  real(KINDR), allocatable, dimension(:) :: Py

  real(KINDR) :: val1, val2
  real(KINDR), dimension(:), allocatable :: y ! phenotypes
  integer, dimension(:,:), allocatable :: xtemp
  real(KINDR), dimension(:,:), allocatable :: x !incid. matrix for fixed effects
  real(KINDR), dimension(:), allocatable :: temAmat
  real(KINDR), dimension(:), allocatable :: theta, oldtheta
  type (doublePre_Array), dimension(:), allocatable :: raneff
  real(KINDR), allocatable, dimension(:) :: fixEff
  real(KINDR), allocatable, dimension(:) :: P, V
  real(KINDR), allocatable, dimension(:,:) :: Vhat, temp
  !! ================ No defintion after this line ================ !!

  call askYesNoInteger(i, " should the program be verbose? (0:No, 1:Yes) ", 0)
  verbose = (i == 1)
  
  ! getting phenotype file name and reading it
  eStatus = "old"
  call askFileName(phenfile, " filename for phenotypes", status, eStatus)
  if (status(1:1) .eq. "x") then
     write(STDERR, *) "error in openning file ", phenFile
     stop 1
  end if

  call askInteger(nran, "Number of random effects (3 if covariate analysis)")
  if (nran .eq. 3) then
     call askInteger(nfix, "Number of fix effects (atm only 2 is accepted)")
  elseif (nran .eq. 1) then
     call askInteger(nfix, "Number of fix effects (excluding mean)")
     if (nfix > 0) then
        allocate(levels(nfix))
        j = 0
        do i = 1, nfix
           write(msg, '(a,1x,i2)') "Number of levels for effect", i
           call askInteger(k, trim(msg))
           levels(i) = k
           j = j + k
        end do
     end if
     nfix = nfix + 1
  else
     write(STDERR, '(a)') "Error"
     write(STDERR, *) "not implemented for nran = ", nran
     stop 3
  end if
  if ((nran .eq. 3) .and. (nfix .ne. 2)) then
     write(STDERR, '(a)') "Error"
     write(STDERR, *) " unexpected config for nfix and nran"
     stop 3
  end if

  ! counting number of lines
  j = 0 ! number of skipped lines
  empties = 1
  call countNumberLines(phenFile, j, lines, empties, ifail)
  if (ifail .ne. 0) stop 1
  nobs = lines - empties

  ! allocating y (phenotypes) and id (real(KINDR) id of animals) and incidience matrix
  allocate(y(nobs), id(nobs), Xtemp(nobs, (nfix-1)))

  ! reading the data
  open(newUnit = phenFileID, file = phenFile, status = 'old')
  maxid = 0
  if ((nran == 1).and.(nfix == 1)) then
     allocate(X(nobs,nfix))
     do i = 1, nobs
        read(phenFileID, *) id(i), y(i)
        if (maxid < id(i)) maxid = id(i)
     end do
     X(1:nobs, 1) = ONE
  elseif ((nran == 1).and.(nfix > 1)) then
     j = sum(levels)
     allocate(X(nobs, j))
     do i = 1, nobs
        read(phenFileID, *) id(i), Xtemp(i, 1:(nfix-1)), y(i)
        k = 0
        do j = 1, (nfix-1)
           if (xtemp(i,j) .eq. 1) cycle
           X(i, Xtemp(i, j) + k - 1) = ONE
           k = k + levels(j)
        end do
        if (maxid < id(i)) maxid = id(i)
     end do
     nfix = sum(levels)
     X(1:nobs, nfix) = ONE
  elseif (nran == 3) then
     allocate(X(nobs,nfix))
     do i = 1, nobs
        read(phenFileID,*) id(i), X(i,1), y(i)
        X(i, nfix) = ONE
        if (maxid < id(i)) maxid = id(i)
     end do
  else
     write(STDERR, '(a)') 'Error:'
     write(STDERR, *) "Not implemented"
     stop 2
  end if
  close(phenFileID)

  allocate(fixeff(nfix))
  allocate(raneff(nran))
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

  write(msg, '(a28)') "file for relationship matrix"
  eStatus = "old"
  call askFileName(AmatFile, trim(msg), status, eStatus)
  if (status(1:1) .eq. "x") then
     write(STDERR, *) "error in openning file ", AmatFile
     stop 2
  end if

  ! counting number of lines
  open(newunit=AmatFileID, file=AmatFile, err=73, status='old')
  do 
     read(AmatFileID,*,end=74) i, j, val1
     if (i > maxid) maxid = i
     if (j > maxid) maxid = j
  end do
73 write(STDERR, *) "error in reading file ", AmatFile
  stop 1
74 continue
  close(AmatFileID)

  i = (maxid + 1) * maxid / 2
  allocate(temAmat(i))

  nelement = (nobs + 1) * nobs / 2
  allocate(P(nelement), V(nelement))

  j = 0 ! file is not binary
  k = 0 ! number of lines to skip
  call trsmReadMat(AmatFile, temAmat, maxid, k, ifail, j)
  if (verbose) write(STDOUT, *) " end reading files"

  THETACOND: if (nran .eq. 3) then
     allocate(oldtheta(5))
     write(STDOUT, '(a27)') "initial guess for variances"
     write(STDOUT, '(3x, a23)', advance = 'no') "genetic part of slope: "
     read(STDIN, *) oldtheta(1)
     write(STDOUT, '(3x, a27)', advance = 'no') "genetic part of intercept: "
     read(STDIN, *) oldtheta(2)
     write(STDOUT, '(3x, a63)', advance = 'no') &
          "correlation between slope and intercept (0.0 if there is not): "
     read(STDIN, *) val1
     if ((val1 > 1.d0 .or. val1 < -1.d0)) then
        write(STDERR, *) "invalid value for correlation. The program will stop"
        stop 1
     end if
     ! theta contains variances and covaraince only; 
     ! hence the correlation must be converted to covariance
     oldtheta(4) = val1 * sqrt(oldtheta(1) * oldtheta(2))
     write(STDOUT, '(3x, a30)', advance = 'no') "environmental variance slope: "
     read(STDIN, *) oldtheta(3)
     write(STDOUT, '(3x, a34)', advance = 'no') &
          "environmental variance intercept: "
     read(STDIN, *) oldtheta(5)

     ! depending on the given correlation, we may need 3 or 4 ZGZ matrices. 
     ! So better to check that because one matrix makes a lot of difference in 
     ! using the amount of memory
     if (oldtheta(4) == 0.d0) then 
        nvar = 3
        allocate(theta(nvar + 1))
        theta(1 : nvar) = oldtheta(1 : nvar)
        theta(nvar + 1) = oldtheta(nvar + 2)
        deallocate(oldtheta)
        allocate(oldtheta(nvar + 1))
        write(STDOUT, '(2x,a22)') "no correlation assumed"
     else
        nvar = 4
        allocate(theta(nvar + 1))
        theta(1 : (nvar + 1)) = oldtheta(1 : (nvar + 1))
        write(STDOUT, '(2x,a30)') "correlation taken into account"
     end if
  elseif (nran == 1) then
     allocate(oldtheta(2), theta(2))
     val1 = sum(y) / size(y)
     val2 = sum((y - val1) ** 2) / (size(y) - 1)
     write(STDOUT, '(a27)') "initial guess for variances"
     write(STDOUT, '(3x, a20)', advance = 'no') "heritability"
     read(STDIN, *) val1
     oldtheta(1) = val2 * val1
     oldtheta(2) = (1-val1) * val2
     nvar = 1
     theta(1 : (nvar + 1)) = oldtheta(1 : (nvar + 1))
  end if THETACOND


  estatus = 'unknown'
  call askFileName(fixEfffile, " filename for fixed effects", status, eStatus)
  call askFileName(ranEfffile, " filename for random effects", status, eStatus)
  call askFileName(varFile, " filename for variances", status, eStatus)

  call askYesNoInteger(doreml, " is a reml required (1:Yes, 0:No)?", 1)

  allocate(Vhat(nfix, nobs))
  allocate(temp(nobs, nfix))
  allocate(ipiv(nobs))
  allocate(py(nobs))

  if (doreml == 1) then
     call askInteger(emiteration, " number of EM iterations (in case reml): ")
     call askInteger(maxIter, &
          " max number of iterations (in case reml, > 5 + em): ")
  end if

  if (doreml == 1) then
     call Reml(id, X, y, nfix, nobs, maxid, nelement, temAmat, nvar, nran,&
          theta, verbose, ipiv, Py, P, V, Vhat, temp, &
          emIterations = emIteration, maxIters = maxIter)
  end if
  call Blup(id, X, y, nfix, nobs, maxid, nelement, temAmat, nvar, nran,&
       theta, fixEff, ranEff, verbose, ipiv, Py, P, V, Vhat, temp)


  if (verbose) write(STDOUT, *) 'fixed effects: ' , fixeff(1 : nfix)

  write(formato, 268) "((", (nfix-1), "(g24.15, 1x), g24.15))"
268 format(a2, i1, a22)
270 format(a2, i1, a21)
271 format(a1, i1, a4)
272 format(i12, 1x, g24.15)
  open(newUnit = iunfix, file = fixEffFile)
  open(newUnit = iunvar, file = varFile)
  write(iunfix, '(2a24)') "fixed effects"
  write(iunfix, trim(formato)) fixeff(1 : nfix)
  close(iunfix)

  write(formato, 271) "(", (nvar+1), "a24)"
  if (nfix == 2) then
     if (nvar .eq. 4) then
        write(iunvar, formato) "var_A_slope","var_A_intercept", &
             "corr(A_int,A_slope)", "var_E_slope", "var_E_intercept"
        write(formato, 270) "(", nvar, "(f24.15, 1x), f24.15)"
        write(iunvar, formato)  theta(1:2), theta(4) / sqrt(theta(1) &
             * theta(2)) , theta(3), theta(5)
     elseif (nvar .eq. 3) then
        write(iunvar, formato) "var_A_slope","var_A_intercept", &
             "var_E_slope", "var_E_intercept"
        write(formato, 270) "(", nvar, "(f24.15, 1x), f24.15)"
        write(iunvar, formato)  theta(1:4)
     end if
  elseif (nfix == 1) then
     write(iunvar, formato) "genetic","residual"
     write(formato, 270) "(", nvar, "(g24.15, 1x), g24.15)"
     write(iunvar, formato)  theta(1:2)
  end if
  close(iunvar)

  open(newUnit = iunran, file = raneffFile)
  do i = 1, maxid
     write(iunran, 272) i, raneff(1)%level(i)
  end do
  if (nfix == 2) then
     do i = 1, maxid
        write(iunran, 272) i, raneff(2)%level(i)
     end do
     do i = 1, nobs
         write(iunran, 272) i, raneff(3)%level(i)
     end do
  end if
  close(iunran)
end program analysis
