program AIREML
  use constants
  use global_module
  use quickSort
  use blup_module
  use reml_module
  implicit none
  !! ================ variable definitions  ================ !!
  character(LEN=256)                                        :: phenFile, AmatFile, msg
  character(len=20)                                         :: status, eStatus

  logical                                                   :: verbose

  integer                                                   :: ifail
  integer                                                   :: i, j, nobs, nfix, maxid!,niter, iter, ivar 
  integer                                                   :: phenFileID!, inciFileID, AmatFileID
  integer                                                   :: skip, lines, empties
  !  integer, dimension(8)                                     :: clock_start, clock_beginning, clock_elements1, clock_elements2
  integer, dimension(:), allocatable                        :: id, tempInd
  double precision                                          :: phenVar, val1, val2
  !  double precision                                          :: ypy
  !  double precision, dimension(:), allocatable               :: gvariances
  double precision, dimension(:), allocatable               :: y, Amatrix, temAmat
  double precision, dimension(:,:), allocatable             :: x 
  !  double precision, dimension(:), allocatable               :: Py, temamat1, raneff, xtvinv
  !  double precision, dimension(:), pointer                   :: a_trmat1,a_trmat2
  !  double precision, dimension(:), allocatable               :: yext, theta, newtheta

  !  type (doublePre_Array), dimension(:), allocatable, target :: theGy

  !  double precision                                          :: determinantV, logl
  !  integer                                                   :: ineg, info
  integer, dimension(:), allocatable                        :: env, temp!, ipiv
  !  double precision, dimension(:,:), allocatable             :: temp1, temp2
  !  integer, dimension(:), allocatable                        :: temp1shape, temp2shape
  !  integer                                                   :: temp1dim, temp2dim
  !  double precision, external                                :: dnrm2, ddot

  !! ================ No defintion after this line ================ !!
  ! getting phenotype file name and reading it
  eStatus = "old"
  call askFileName(phenfile, " filename for phenotypes", status, eStatus)
  if (status(1:1) .eq. "x") then
     write(STDERR, *) "error in openning file ", phenFile
     stop 1
  end if

  call askInteger(skip, " Number of lines to skip: ")
  ! counting number of lines
  j = 0 ! number of skipped lines
  empties=1
  call countNumberLines(phenFile, skip, lines, empties, ifail)
  if (ifail .ne. 0) stop 1
  nobs=lines-empties

  ! allocating y (phenotypes) and id (real id of animals)
  allocate(y(nobs), id(nobs), env(nobs), tempInd(nobs), temp(nobs))

  open(newUnit = phenFileID, file = phenFile)

  ! reading the data and calculating variance
  val1 = 0.d0
  val2 = 0.d0
  maxid = 0
  do i = 1, nobs
     read(phenFileID,*) id(i), y(i), env(i)
     val1 = val1 + y(i)
     val2 = val2 + y(i) * y(i)
     if (maxid < id(i)) maxid = id(i)
  end do
  phenVar = (val2 - (val1 * val1) / dble(nobs)) / dble(nobs-1)
  close(phenFileID)

  ! counting distinct values in the environment (nfix)
  call sortix(nobs, env, tempInd)

  i = 1
  j = 1
  temp(j) = env(tempInd(i))
  do i = 2, nobs
     if (env(tempInd(i)) .ne. env(tempInd(i-1))) then
        j = j + 1
        temp(j) = env(tempInd(i))
     end if
  end do
  nfix = j

  write(STDOUT, *) "number of fixed effects are ", nfix

  allocate(X(nobs, nfix))
  X(1 : nobs, 1 : nfix) = 0

  do i = 1, nobs
     do j = 1, nfix
        if (temp(j) == env(i)) then
           X(i, j) = 1
           exit
        end if
     end do
  end do

  ! reading g matrix(matrices)
  i = maxid * (maxid + 1) / 2
  allocate(temAmat(i))
  i = nobs * (nobs + 1) / 2
  write(msg, *) "file for relationship matrix:"
  eStatus = "old"
  call askFileName(AmatFile, trim(msg), status, eStatus)
  j = 0 ! the file is not binary
  skip = 0
  call trsmReadMat(AmatFile, temAmat, maxid, skip, ifail, j)
  allocate(Amatrix(i))
  call getMatrices(verbose, nobs, X, temAmat, id, Amatrix)
  write(STDOUT, *) " end reading files"
  
  ! allocating some variables based on the nvar

  !
  !  ! time printing for some reason
  !  call printingDateTime(6,1,clock_elements2)
  !  call printingElapseTime(6,clock_elements1,clock_elements2)
  !
  !  ! initial value for ratio of variances
  !  
  !  !! fixing variances based on the options
  !  val1 = 0.1d0
  !  gvariances(:) = phenVar * 0.1d0
  !  
  !  write(6,*) " starting values "
  !  do i = 1, nvar
  !     write(6,*) " gvariances ",i,gvariances(i)
  !  end do
  !  write(6,*) " V_P (phenVar): ", 0, phenVar
  !  write(6,*) " total", 0, phenVar + sum(gvariances)
  !  write(6,*)
  !
  !  write(6,*) " calculatng logL for staring value"
  !  call printingDateTime(6,0,clock_elements1)
  !  clock_start=clock_elements1
  !
  !  i = nobs * nobs
  !
  !  allocate(theta(nvar + 1), newtheta(nvar + 1))
  !  allocate(raneff(nobs), Py(nobs), xtvinv(nobs * nfix))
  !  theta(1) = merge((phenvar - sum(gvariances)) , phenvar  , phenvar > sum(gvariances))
  !  theta(2:) = gvariances(:)
  !  newtheta(:) = theta(:)
  !  verbose = .true.
  !  i = 0
  !  do 
  !     i = i + 1
  !     theta(:) = newtheta(:)
  !
  !     call iterate(nobs, nvar, nfix, theZGZ, y, x, logl , theta, Py, xtvinv, verbose)
  !     newtheta = newtheta - theta
  !     val1 = dnrm2(nvar + 1, newtheta, 1)
  !     write(6, *) ' iteration:',i, 'error is', val1
  !     if (val1 < 1.d-6) then
  !        write(6, *) 'converged!'
  !        exit
  !     elseif (i > 100) then
  !        write(6, *) ' not going to coverge :('
  !        exit
  !     end if
  !     newtheta(:) = theta(:)
  !  end do
  !  
  !  ! Quick and Dirty implementation for getting random and fixed effects
  !  call dspmv('u', nobs, theta(1 + 1),  theZGZ(1)%level, Py, 1, 0.d0, raneff, 1)
  !
  !  
  !  open(1, file = 'raneff.txt')
  !  do i = 1, nobs
  !     write(1, '(i4,2x, g25.16)') id(i), raneff(i)
  !  end do
  !  close(1)
  !  open(1, file = 'fixeff.txt')
  !  write(6,  *) " fixed effect value is ", ddot(nobs, xtvinv, 1, y, 1) / ddot(nobs, xtvinv, 1, x, 1)
  !  write(1,  *) " fixed effect value is ", ddot(nobs, xtvinv, 1, y, 1) / ddot(nobs, xtvinv, 1, x, 1)
  !  close(1)
  !  
end program AIREML
