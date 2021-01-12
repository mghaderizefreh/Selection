! order of matrices, and therefore variances shall be 
!   1     ZsGZs               genetic slope
!   2     ZiGZi               genetic intercept
!   3     ZsZs                environmental slope (stored as diagonal)
!   4     ZiGZs+ZsGZi         genetic slope-intercept covariance
!   5     ZsZs                permanent environment effect slope (SAME AS NUMBER 3: but counted)
!   6     ZiZi                permanent environment effect intercept
!   7     ZsZi+ZiZs           permanent environment effect slope-intercept covariance
! (last)  Identity            environmental intercept (NOT COUNTED IN theZGZ and it is always LAST one)
program doRRReml
  use constants
  use global_module
  use blup_module
  use reml_module
  use RR_reml
  implicit none
  !! ================ variable definitions  ================ !!
  character(LEN=256)                                  :: phenFile, AmatFile, fixEffFile, ranEffFile, varFile, msg
  character(len=30)                                   :: status, eStatus, formato
  logical                                             :: verbose = .false.
  integer                                             :: ifail
  integer                                             :: i, j, k, maxid, nvar, nobs, nfix
  integer                                             :: phenFileID, AmatFileID, iunFix, iunRan, iunVar
  integer                                             :: lines, empties, emiteration
  integer, dimension(:), allocatable                  :: id ! real id of animals

  double precision                                    :: val1
  double precision, dimension(:), allocatable         :: y ! phenotypes
  double precision, dimension(:,:), allocatable       :: x ! incidence matrix for fixed effects
  double precision, dimension(:), allocatable         :: temAmat
  double precision, dimension(:), allocatable         :: theta, oldtheta
  type (doublePre_Array), dimension(:), allocatable   :: raneff
  double precision, allocatable, dimension(:)         :: fixEff
  !! ================ No defintion after this line ================ !!

  ! getting phenotype file name and reading it
  eStatus = "old"
  call askFileName(phenfile, " filename for phenotypes", status, eStatus)
  if (status(1:1) .eq. "x") then
     write(stderr, *) "error in openning file ", phenFile
     stop 1
  end if

  nfix = 2

  ! counting number of lines
  j = 0 ! number of skipped lines
  empties=1
  call countNumberLines(phenFile, j, lines, empties, ifail)
  if (ifail .ne. 0) stop 1
  nobs=lines-empties

  ! allocating y (phenotypes) and id (real id of animals) and incidience matrix
  allocate(y(nobs), id(nobs), X(nobs,nfix))

  ! reading the data
  open(newUnit = phenFileID, file = phenFile, status = 'old')
  maxid = 0
  do i = 1, nobs
     read(phenFileID,*) id(i), X(i,1), y(i)
     X(i,2) = 1.d0
     if (maxid < id(i)) maxid = id(i)
  end do
  close(phenFileID)

  write(msg, '(a28)') "file for relationship matrix"
  eStatus = "old"
  call askFileName(AmatFile, trim(msg), status, eStatus)
  if (status(1:1) .eq. "x") then
     write(stderr, *) "error in openning file ", AmatFile
     stop 1
  end if

  ! counting number of lines
  open(newunit=AmatFileID, file=AmatFile, err=73, status='old')
  do 
     read(AmatFileID,*,end=74) i, j, val1
     if (i > maxid) maxid = i
     if (j > maxid) maxid = j
  end do
73 write(stderr, *) "error in reading file ", AmatFile
  stop 1
74 continue
  close(AmatFileID)

  i = (maxid + 1) * maxid / 2 
  allocate(temAmat(i))

  j = 0 ! file is not binary
  k = 0 ! number of lines to skip
  call trsmReadMat(AmatFile, temAmat, maxid, k, ifail, j)

  if (verbose) write(stdout, *) " end reading files"

  allocate(oldtheta(5))
  write(stdout, '(a27)') "initial guess for variances"
  write(stdout, '(3x, a23)', advance = 'no') "genetic part of slope: "
  read(stdin, *) oldtheta(1)
  write(stdout, '(3x, a27)', advance = 'no') "genetic part of intercept: "
  read(stdin, *) oldtheta(2)
  write(stdout, '(3x, a63)', advance = 'no') &
       "correlation between slope and intercept (0.0 if there is not): "
  read(stdin, *) val1
  if ((val1 > 1.d0 .or. val1 < -1.d0)) then
     write(stderr, *) "invalid value for correlation. The program will stop"
     stop 1
  end if
  ! theta contains variances and covaraince only; 
  ! hence the correlation must be converted to covariance
  oldtheta(4) = val1 * sqrt(oldtheta(1) * oldtheta(2))
  write(stdout, '(3x, a30)', advance = 'no') "environmental variance slope: "
  read(stdin, *) oldtheta(3)
  write(stdout, '(3x, a34)', advance = 'no') &
       "environmental variance intercept: "
  read(stdin, *) oldtheta(5)

  ! depending on the given correlation, we may need 3 or 4 ZGZ matrices. So better to check that
  ! because one matrix makes a lot of difference in using the amount of memory
  if (oldtheta(4) == 0.d0) then 
     nvar = 3
     allocate(theta(nvar + 1))
     theta(1 : nvar) = oldtheta(1 : nvar)
     theta(nvar + 1) = oldtheta(nvar + 2)
     deallocate(oldtheta)
     allocate(oldtheta(nvar + 1))
     write(stdout, '(2x,a22)') "no correlation assumed"
  else
     nvar = 4
     allocate(theta(nvar + 1))
     theta(1 : (nvar + 1)) = oldtheta(1 : (nvar + 1))
     write(stdout, '(2x,a30)') "correlation taken into account"
  end if

  estatus = 'unknown'
  call askFileName(fixEfffile, " filename for fixed effects", status, eStatus)
  call askFileName(ranEfffile, " filename for random effects", status, eStatus)
  call askFileName(varFile, " filename for variances", status, eStatus)
  !  fixEffFile = "fixedEffects"
  !  ranEffFile = "randomEffects"
  !  varFile = "variances"
  call askInteger(emiteration, "number of EM iterations: ")

  allocate(fixEff(2), raneff(3))
  allocate(raneff(1)%level(maxid)) ! slope effect (genetic)
  allocate(raneff(2)%level(maxid)) ! intercept effect (genetic)
  allocate(raneff(3)%level(nobs))   ! environment slope effect (diagonal)

  call RRReml(id, X, y, nfix, nobs, maxid, temAmat, nvar, theta, &
       fixEff, ranEff, verbose, emIterations = emIteration)

  if (verbose) write(stdout, *) 'fixed effects: ' , fixeff(1 : nfix)
268 format(a2, i1, a22)
  write(formato, 268) "((", (nfix-1), "(g24.15, 1x), g24.15))"
  open(newUnit = iunfix, file = fixEffFile)
  write(iunfix, '(2a24)') "slope", "intercept"
  write(iunfix, trim(formato)) fixeff(1 : nfix)
  close(iunfix)

  open(newUnit = iunvar, file = varFile)
270 format(a2, i1, a21)
271 format(a1, i1, a4)
  write(formato, 271) "(", (nvar+1), "a24)"
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
  close(iunvar)

272 format(i12, 1x, g24.15)

  open(newUnit = iunran, file = raneffFile)
  do i = 1, maxid
     write(iunran, 272) i, raneff(1)%level(i)
  end do
  do i = 1, maxid
     write(iunran, 272) i, raneff(2)%level(i)
  end do
  do i = 1, nobs
     write(iunran, 272) i, raneff(3)%level(i)
  end do
  close(iunran)

end program DORRREML
