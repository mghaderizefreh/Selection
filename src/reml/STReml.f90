module AIREMLmodule
  implicit none

  contains
    include "createZGtZ.f90"
    include "AIRemlCalcVmat.f90"
    include "iteration.f90"

end module AIREMLmodule

program AIREML
  use constants
  use global_module
  use trsm_module
  use AIREMLmodule
  implicit none
  !! ================ variable definitions  ================ !!
  character(LEN=256)                                                  :: phenFile, inciFile, AmatFile
  character(len=20)                                                   :: status, eStatus, msg

  logical                                                             :: isorted, verbose

  integer                                                             :: ifail
  integer                                                             :: i, j, k, niter, iter, maxid, nvar, ivar, nobs, nfix
  integer                                                             :: phenFileID, inciFileID, AmatFileID
  integer                                                             :: skip, lines, empties
  integer, dimension(8)                                               :: clock_start, clock_beginning, clock_elements1, clock_elements2   ! array must be length 8
  integer, dimension(:), allocatable                                  :: id, tempInd

  double precision                                                    :: phenVar
  double precision                                                    :: val1, val2, ypy
  double precision, dimension(:), allocatable                         :: gvariances
  double precision, dimension(:), allocatable                         :: y , env, tempEnv
  double precision, dimension(:), allocatable                         :: x ! incidence matrix
  double precision, dimension(:), allocatable                         :: temAmatFirst, Py, temamat1, raneff, xtvinv
  double precision, dimension(:), pointer                             :: a_trmat1,a_trmat2
  double precision, dimension(:), allocatable                         :: error, yext, theta, newtheta

  type (doublePre_Array), dimension(:), allocatable, target           :: theZGZ, theGy

  double precision                                                    :: determinantV, logl
  integer                                                             :: ineg, info
  integer, dimension(:), allocatable                                  :: ipiv
  double precision, dimension(:,:), allocatable                       :: temp1, temp2
  integer, dimension(:), allocatable                                  :: temp1shape, temp2shape
  integer                                                             :: temp1dim, temp2dim
  double precision, external                                          :: dnrm2, ddot

  !! ================ No defintion after this line ================ !!
  ! getting phenotype file name and reading it
  eStatus = "old"
  call askFileName(phenfile, " filename for phenotypes", status, eStatus)
  open(newUnit = phenFileID, file = phenfile, status = eStatus)
  call askInteger(skip, " Number of lines to skip: ")
  if (status(1:1) .eq. "x") then
     write(stderr, *) "error in openning file ", phenFile
     stop 1
  end if

  ! counting number of lines
  j = 0 ! number of skipped lines
  empties=1
  call countNumberLines(phenFile, skip, lines, empties, ifail)
  if (ifail .ne. 0) stop 1
  nobs=lines-empties

  ! allocating y (phenotypes) and id (real id of animals)
  allocate(y(nobs), id(nobs), env(nobs), tempInd(nobs), tempEnv(nobs))

  ! reading the data and calculating variance
  val1 = 0.d0
  val2 = 0.d0
  isorted = .true.
  maxid = 0
  do i = 1, nobs
     read(phenFileID,*) id(i), env(i), y(i)
     val1 = val1 + y(i)
     val2 = val2 + y(i) * y(i)
     if (id(i) .ne. i) then
        isorted = .false.  !i.e., the phenotypes are not sorted
     end if
     if (maxid < id(i)) maxid = id(i)
  end do
  phenVar = (val2 - (val1 * val1) / dble(nobs)) / dble(nobs-1)
  close(phenFileID)
  
  ! counting distinct values in the environment (nfix)
  call sortdx(nobs, env, tempInd)
  j = 1
  tempEnv(j) = env(tempInd(j))
  do i = 2, nobs
     if (env(tempInd(i)) .ne. env(tempInd(i-1))) then
        j = j + 1
        tempEnv(j) = env(tempInd(j))
     end if
  end do
  nfix = j
  write(6, *) "number of fixed effects are ", nfix

  allocate(X(nobs, nfix))

!!here here here -> to do: popultate x array; uncomment rest and debug; re-use reml subroutines to solve the problem
  

  ! reading the incidence matrix (same way as Ricardo's code)
  eStatus = "old"
  call askFileName(inciFile," incidence matrix filename",status, eStatus)
  open(newUnit = inciFileID, file = inciFile, status = eStatus)

  read(inciFileID,*) i, nfix
  j = i * nfix
  write(6,*) " nobs i nfix j", nobs, i, nfix, j
  
  ! allocating the incidence matrix
  ! TODO: shouldn't the incidince matrix be two dimensional? what is going on here
  allocate(X(j))
  read(inciFileID,*) X(:)  !correct as element are stored row first
  close(inciFileID)

  ! reading g matrix(matrices)
  call askInteger(nvar, " nvar: ")
  allocate(theZGZ(nvar))
  if (isorted) then
     write(6,*) " phenotype  sorted.... G matrix is ZGZ"
     i = trsmCalcSize(nobs) ! i = (nobs + 1) * nobs / 2
     do ivar = 1, nvar
        write(msg, '(a13,i1,a2)') "file for ZGZ(",ivar,"):"
        eStatus = "old"
        call askFileName(AmatFile,trim(msg),status, eStatus)
        call askYesNoInteger(j, " binary file?", 0)
        allocate(theZGZ(ivar)%level(i))
        skip = 0
        call trsmReadMat(AmatFile, temAmatFirst, nobs, skip, ifail, j)
        theZGZ(ivar)%level = temAmatFirst
     end do
  else
     write(6,*) " phenotype not sorted...  G matrix is not ZGZ"
     i = trsmCalcSize(maxid)
     allocate(temAmatFirst(i))
     i = trsmCalcSize(nobs)
     do ivar = 1, nvar
        write(msg, '(a11,i1,a2)') "file for G(", ivar, "):"
        eStatus = "old"
        call askFileName(AmatFile, trim(msg), status, eStatus)
        call askYesNoInteger(j, " binary file?", 0)
        skip = 0
        call trsmReadMat(AmatFile, temAmatFirst, maxid, skip, ifail, j)
        allocate(theZGZ(ivar)%level(i))
        call createZGtZ(id, temAmatFirst, nobs, theZGZ(ivar)%level)
     end do
  end if
  write(6,*) " end reading files"

  ! allocating some variables based on the nvar
  allocate(gvariances(nvar))

  ! time printing for some reason
  call printingDateTime(6,1,clock_elements2)
  call printingElapseTime(6,clock_elements1,clock_elements2)

  ! initial value for ratio of variances
  
  !! fixing variances based on the options
  val1 = 0.1d0
  gvariances(:) = phenVar * 0.1d0
  
  write(6,*) " starting values "
  do i = 1, nvar
     write(6,*) " gvariances ",i,gvariances(i)
  end do
  write(6,*) " V_P (phenVar): ", 0, phenVar
  write(6,*) " total", 0, phenVar + sum(gvariances)
  write(6,*)

  write(6,*) " calculatng logL for staring value"
  call printingDateTime(6,0,clock_elements1)
  clock_start=clock_elements1

  i = nobs * nobs

  allocate(theta(nvar + 1), newtheta(nvar + 1))
  allocate(raneff(nobs), Py(nobs), xtvinv(nobs * nfix))
  theta(1) = merge((phenvar - sum(gvariances)) , phenvar  , phenvar > sum(gvariances))
  theta(2:) = gvariances(:)
  newtheta(:) = theta(:)
  verbose = .true.
  i = 0
  do 
     i = i + 1
     theta(:) = newtheta(:)

     call iterate(nobs, nvar, nfix, theZGZ, y, x, logl , theta, Py, xtvinv, verbose)
     newtheta = newtheta - theta
     val1 = dnrm2(nvar + 1, newtheta, 1)
     write(6, *) ' iteration:',i, 'error is', val1
     if (val1 < 1.d-6) then
        write(6, *) 'converged!'
        exit
     elseif (i > 100) then
        write(6, *) ' not going to coverge :('
        exit
     end if
     newtheta(:) = theta(:)
  end do
  
  ! Quick and Dirty implementation for getting random and fixed effects
  call dspmv('u', nobs, theta(1 + 1),  theZGZ(1)%level, Py, 1, 0.d0, raneff, 1)

  
  open(1, file = 'raneff.txt')
  do i = 1, nobs
     write(1, '(i4,2x, g25.16)') id(i), raneff(i)
  end do
  close(1)
  open(1, file = 'fixeff.txt')
  write(6,  *) " fixed effect value is ", ddot(nobs, xtvinv, 1, y, 1) / ddot(nobs, xtvinv, 1, x, 1)
  write(1,  *) " fixed effect value is ", ddot(nobs, xtvinv, 1, y, 1) / ddot(nobs, xtvinv, 1, x, 1)
  close(1)
  
end program AIREML
