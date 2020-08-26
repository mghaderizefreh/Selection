program blup
  use constants
  use global_module
  use blup_module
  use reml_module
  implicit none


  character(LEN=256)                                                  :: phenFile, AmatFile, fixEffFile, ranEffFile, varFile, msg
  character(len=20)                                                   :: status, eStatus
  logical                                                             :: verbose = .false.
  integer                                                             :: i, j, k, maxid, nvar, nobs, nfix, ifail
  integer                                                             :: phenFileID, AmatFileID
  integer                                                             :: lines, empties
!  integer, dimension(8)                                               :: clock_beginning, clock_elements1, clock_elements2, diff_elements   ! array must be length 8
  integer, dimension(:), allocatable                                  :: id ! real id of animals
  double precision                                                    :: val1, val2
  double precision, dimension(:), allocatable                         :: y ! phenotypes
  double precision, dimension(:,:), allocatable                       :: Vhat, x ! incidence matrix for fixed effects
  double precision, dimension(:), allocatable                         :: temAmat, Py
  double precision, dimension(:), allocatable                         :: theta, oldtheta

  type (doublePre_Array), dimension(:), allocatable                   :: theZGZ
  integer, dimension(:), allocatable, save                            :: ipiv
  double precision, dimension(:), allocatable, save                   :: V, P, work
  double precision, external                                          :: dnrm2, ddot, dasum

  nfix = 2 ! this is because I look for the mean intercept and the mean slope

  ! getting phenotype file name and reading it
  eStatus = "old"
  call askFileName(phenfile, " phenotypic filename", status, eStatus)
  if (status(1:1) .eq. "x") then
     write(stderr, *) "error in openning file ", phenFile
     stop 1
  end if

  ! counting number of lines
  j = 0 ! number of skipped lines
  empties=1
  call countNumberLines(phenFile, j, lines, empties, ifail)
  if (ifail .ne. 0) stop 1
  nobs=lines-empties
  write(6, *) nobs
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

  if (verbose) write(6, *) " end reading files"
  allocate(Py(nobs), Vhat(nfix, nobs))
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

  ! making G* matrices
  allocate(theZGZ(nvar))
  i = nobs * (nobs + 1) / 2
  do j = 1, nvar 
     if (j .eq. 3) then
        allocate(theZGZ(j)%level(nobs))
     else        
        allocate(theZGZ(j)%level(i))
     end if
  end do
  if (nvar .eq. 3) then
     call getMatricesUncorrelated(verbose, nobs, X, temAmat, id, &
          theZGZ(1)%level, theZGZ(2)%level, theZGZ(3)%level)
  else
     call getMatricesCorrelated(verbose, nobs, X, temAmat, id, &
          theZGZ(1)%level, theZGZ(2)%level, theZGZ(3)%level, &
          theZGZ(4)%level)
  end if
  deallocate(temAmat)

69 format(a12, i3)
70 format(1x, a9)
71 format(1x, a10, 1x, f24.15, a10, 1x, f24.15)
  if (verbose) then
     write(stdout, *) 
     write(stdout, 69) "iteration: " ,0
     write(stdout, 70, advance='no') " theta_0:"
     write(stdout, *) theta(1 : (nvar + 1))
  end if

  eStatus = "u"
  call askFileName(fixEfffile, " filename for fixed effects", status, eStatus)
  call askFileName(ranEfffile, " filename for random effects", status, eStatus)
  call askFileName(varFile, " filename for variances", status, eStatus)
  !  fixEffFile = "fixedEffects"
  !  ranEffFile = "randomEffects"
  !  varFile = "variances"

  I = nobs * (nobs + 1) / 2
  allocate(P(I),V(I))
  I = nobs * nobs
  allocate(work(I),ipiv(nobs))

  call calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose)
  if (verbose) write(6, *) " V is calculated"

  call detInv(nobs, V, val1, ipiv, work, verbose)
  if (verbose) write(6, *) " V is replaced by its inverse"

  call calculateP(nobs, nfix, V, X, P, val2, Vhat, verbose)
  if (verbose) write(6, *) " P is calcuated"

  call dspmv('u', nobs, 1.d0, P, y, 1, 0.d0, Py, 1)
  if (verbose) write(6, *) "  DSPMV finished calculating Py (=P * y)"

  do i = 1, nvar
     deallocate(theZGZ(i)%level)
  end do
  deallocate(theZGZ)

  call getEffects(nobs, maxid, nfix, nvar, fixeffFile, raneffFile, &
       varFile, theta, AmatFile, Vhat, Py, y, X, id, verbose)

end program blup
