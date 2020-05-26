! order of matrices, and therefore variances shall be 
!   1     ZsGZs               genetic slope
!   2     ZiGZi               genetic intercept
!   3     ZsZs                environmental slope (stored as diagonal)
!   4     ZiGZs+ZsGZi         genetic slope-intercept covariance
!   5     ZsZs                permanent environment effect slope (SAME AS NUMBER 3: but counted)
!   6     ZiZi                permanent environment effect intercept
!   7     ZsZi+ZiZs           permanent environment effect slope-intercept covariance
! (last)  Identity            environmental intercept (NOT COUNTED IN theZGZ and it is always LAST one)
program RRREML
  use constants
  use global_module
  use blup_module
  use reml_module
  use iteration
  implicit none
  !! ================ variable definitions  ================ !!
  character(LEN=256)                                  :: phenFile, AmatFile, fixEffFile, ranEffFile, varFile, msg
  character(len=20)                                   :: status, eStatus
  logical                                             :: verbose = .false.
  integer                                             :: ifail, maxiter = 20
  integer                                             :: i, j, k, maxid, nvar, nobs, nfix
  integer                                             :: phenFileID, AmatFileID
  integer                                             :: lines, empties
!  integer, dimension(8)                               :: clock_beginning, clock_elements1, clock_elements2,&
!                                                                         diff_elements   ! array must be length 8
  integer, dimension(:), allocatable                  :: id ! real id of animals

  double precision                                    :: logl, epsilon = 1.d-6
  double precision                                    :: val1, val2
  double precision, dimension(:), allocatable         :: y ! phenotypes
  double precision, dimension(:,:), allocatable       :: Vhat ,x ! incidence matrix for fixed effects
  double precision, dimension(:), allocatable         :: temAmat, Py
  double precision, dimension(:), allocatable         :: theta, oldtheta

  type (doublePre_Array), dimension(:), allocatable   :: theZGZ

  double precision, external                          :: dnrm2, ddot, dasum
  !! ================ No defintion after this line ================ !!

  nfix = 2 ! this is because I look for the mean intercept and the mean slope

  ! getting phenotype file name and reading it
  eStatus = "old"
  call askFileName(phenfile, " filename for phenotypes", status, eStatus)
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

  oldtheta(1 : (nvar + 1)) = theta(1 : (nvar + 1))
  i = 0
  do 
!     if (verbose) call date_and_time(values = clock_elements1)

     i = i + 1

     call iterate(nobs, nvar, nfix, theZGZ, y, x, logl , theta, Py, Vhat, verbose)

     val1 = dnrm2(nvar + 1, oldtheta, 1)
     oldtheta(1 : (nvar + 1)) = oldtheta(1 : (nvar + 1)) - theta(1 : (nvar + 1))
     val2 = dnrm2(nvar + 1, oldtheta, 1) / val1
     val1 = dasum(nvar + 1, oldtheta, 1) / (nvar + 1)

     write(stdout, *) 
     write(stdout, 69) "iteration: ",i
     write(stdout, 71, advance='no') " l1 error:", val1 ,"; l2 error:", val2

!     if (verbose) then
!        call date_and_time(values = clock_elements2)
!        write(stdout, 72, advance='no') " iteration time: "
!        call getTimeDiff(clock_elements1,clock_elements2,diff_elements)
!        write(stdout, 100, advance = 'no') diff_elements(5:8)
!     end if
     write(stdout, *) 

     if ((val1 < sqrt(epsilon)) .or. (val2 < epsilon)) then
        write(stdout, '(a10)') "converged!"
        exit
     elseif (i > maxiter) then
        write(stdout, '(a16)') "did not converge"
        stop 1
     end if
     oldtheta(1 : (nvar + 1)) = theta(1 : (nvar + 1))
  end do
!100 format (i2,":",i2.2,":",i2.2,".",i3.3)

  do i = 1, nvar
     deallocate(theZGZ(i)%level)
  end do
  deallocate(theZGZ)

  call getEffects(nobs, maxid, nfix, nvar, fixeffFile, raneffFile, &
       varFile, theta, AmatFile, Vhat, Py, y, X, id, verbose)

  !  if (verbose) then
!  call printingDateTime(6,1,clock_elements1)
!  call printingElapseTime( 6, clock_beginning, clock_elements1)
  ! end if
end program RRREML
