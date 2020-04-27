subroutine getMatricesCorrelated(verbose, nobs, X, G, id, Z1GZ1, Z2GZ2, Z1Z1, Z1GZ2)!, Z1Z1pe, Z2Z2, Z1Z2)
  ! This subroutine calculates ZGZ for different terms
  ! Z1GZ1  slope (=Zs*G*Zs') so Z1 is the incidence matrix with challenge levels
  ! Z2GZ2  intercept (=Zi*G*Zi') so Z2 is the incidnece matrix with ones 
  ! Z1Z1   heterogeneous residual (environmenal slope). It is diagonl (=diag(env)^2) (? and also permanent env effect slope)
  ! Z1GZ2  genetic correlation betweeen slope and intercept (=Zi*G*Zs'+Zs*G*Zi')
  ! Z2Z2   permanent environmental effect (??? it is block diagonal with ones [not sure]) (?=Zi*Zi')
  ! Z1Z2   correlation between permanent environmental effect slope and intercept (=Z1*Z2'+Z2*Z1')
  ! Written by M. Ghaderi Zefreh
  implicit none
  logical                                                               :: verbose
  integer, intent(in)                                                   :: nobs
  integer, dimension(:), intent(in)                                     :: id
  double precision, dimension(:), intent(in)                            :: G
  double precision, dimension(:,:), intent(in)                          :: X
  double precision, dimension(:), intent(out)                           :: Z1GZ1, Z2GZ2, Z1Z1, Z1GZ2
  !  double precision, dimension(:), intent(out), optional                 :: Z1Z1pe, Z2Z2, Z1Z2
  integer                                                               :: i, j, k, ipos, ipos1, ipos2
  intrinsic                                                             :: max
  integer, external                                                     :: lowerpos
  ! logical                                                               :: IsPeIntPresent = .false., IsPeSloPresent = .false. , IsPeCovPresent = .false.

  j = nobs
  i = (j + 1) * j / 2
  !allocate(eye(i))
  !eye(:) = 0.d0
  Z1Z1(1:j) = 0.d0
  Z1GZ1(1:i) = 0.d0
  Z2GZ2(1:i) = 0.d0
  Z1GZ2(1:i) = 0.d0

  !  if (present(Z1Z1pe)) then
  !     Z1Z1pe(:) = 0.d0
  !     IsPeSloPresent = .true.
  !  end if
  !
  !  if (present(Z2Z2)) then
  !     Z2Z2(:) = 0.d0
  !     IsPeIntPresent = .true.
  !  end if
  !
  !  if (present(Z1Z2)) then
  !     Z1Z2(:) = 0.d0
  !     IsPeCovPresent = .true.
  !  end if

  ipos = 0
  k = 0
  do i = 1, nobs
     k = k + i
     !   eye(k) = 1.d0
     Z1Z1(i) = X(i,1) * X(i,1)
     do j = 1, i
        ipos = ipos + 1
        ipos1 = lowerpos(id(i), id(j))
        ipos2 = lowerpos(id(j), id(i))

        Z1GZ1(ipos) = X(i,1) * G(ipos1) * X(j,1)
        Z2GZ2(ipos) = G(ipos1)

        Z1GZ2(ipos) = X(j,1) * G(ipos1) + X(i,1) * G(ipos2)    ! A + transpose(A) 
        !      if (IsPeSloPresent) Z1Z1pe(ipos) = do something
        !      if (IsPeIntPresent) Z2Z2(ipos) = eye(ipos1) ! this is wrong
        !      if (IsPeCovPresent) Z1Z2(ipos) = X(j,1) * eye(ipos1) + X(i,1) * eye(ipos2) ! A + transpose(A) 
     end do
  end do

  if (verbose) write(6, '(a)') " End of getMatrices subroutine"
end subroutine getMatricesCorrelated


subroutine getMatricesUncorrelated(verbose, nobs, X, G, id, Z1GZ1, Z2GZ2, Z1Z1)
  ! This subroutine calculates ZGZ for different terms
  ! Z1GZ1  slope (=Zs*G*Zs') so Z1 is the incidence matrix with challenge levels
  ! Z2GZ2  intercept (=Zi*G*Zi') so Z2 is the incidnece matrix with ones 
  ! Z1Z1   heterogeneous residual (environmenal slope). It is diagonl (=diag(env)^2) (? and also permanent env effect slope)
  ! Z1GZ2  genetic correlation betweeen slope and intercept (=Zi*G*Zs'+Zs*G*Zi')
  ! Z2Z2   permanent environmental effect (??? it is block diagonal with ones [not sure]) (?=Zi*Zi')
  ! Z1Z2   correlation between permanent environmental effect slope and intercept (=Z1*Z2'+Z2*Z1')
  ! Written by M. Ghaderi Zefreh
  implicit none
  logical                                                               :: verbose
  integer, intent(in)                                                   :: nobs
  integer, dimension(:), intent(in)                                     :: id
  double precision, dimension(:), intent(in)                            :: G
  double precision, dimension(:,:), intent(in)                          :: X
  double precision, dimension(:), intent(out)                           :: Z1GZ1, Z2GZ2, Z1Z1
  integer                                                               :: i, j, k, ipos, ipos1, ipos2
  intrinsic                                                             :: max
  integer, external                                                     :: lowerpos

  j = nobs
  i = (j + 1) * j / 2
  Z1Z1(1:j) = 0.d0
  Z1GZ1(1:i) = 0.d0
  Z2GZ2(1:i) = 0.d0

  ipos = 0
  k = 0
  do i = 1, nobs
     k = k + i
     !   eye(k) = 1.d0
     Z1Z1(i) = X(i,1) * X(i,1)
     do j = 1, i
        ipos = ipos + 1
        ipos1 = lowerpos(id(i), id(j))
        ipos2 = lowerpos(id(j), id(i))

        Z1GZ1(ipos) = X(i,1) * G(ipos1) * X(j,1)
        Z2GZ2(ipos) = G(ipos1)
     end do
  end do

  if (verbose) write(6, '(a)') " End of getMatrices subroutine"
end subroutine getMatricesUncorrelated

!================================================================================

subroutine trsmReadMat(matfile,amat,nrank,skip,ifail,ibin)
  implicit none
  ! written by R. Pong-Wong
  ! edited by M. Ghaderi Zefreh (minor edits for stdout/stderr output)
  character(len=*)   , intent(IN) :: matfile
  double precision   , dimension(:), intent(inout) :: amat
  integer            , intent(inout) :: nrank
  integer            , intent(in) :: skip
  integer            , intent(out) :: ifail

  integer :: irow, icol, ipos,i,ibin,k
  integer :: irank,iun
  double precision :: val1
  integer, external :: lowerpos

  if (ibin .eq. 1) then
     open(newunit=iun, file=matfile, status='old', form='unformatted')
     read(iun) irank
     write(6,'(1x,a24,i8)') "rank in unformated file", irank
  else
     open (newunit=iun, file=matfile, status='old')
     do i=1,skip
        read(iun, *)
     end do

     if (nrank > 0) then
        irank = nrank
     else
        read(iun, *) irank
     endif
  endif
  ipos = (irank + 1) * irank / 2
  icol = size(amat, dim = 1)
  if (icol < ipos) then
     write(0, *) " array is not enough to hold matrix of rank"
     write(0, *) "ipos, icol, irank", ipos, icol, irank
     ifail = 1
     close(iun)
     return
  end if
  nrank=irank
  !---------------------------------------------------
  k=0
  i=0
  amat(1:ipos)= 0.d0
  do
     if(ibin==1)then
        read(iun,end=100)irow,icol, val1
     else
        read(iun,*,end=100)irow,icol, val1
     end if
     if(irow > nrank .or. icol > nrank) then
        write(0,*) "error array contains elements in row/column greater than rank"
        i=i+1
        write(0,'(a,i11,5x,i11,i11,g25.15)') 'line element',i,irow, icol, val1
        ifail=1
        close(iun)
        return
     endif
     amat(lowerPos(irow,icol))=val1
     i=i+1
  end do
100 continue
  close(iun)
  !  WRITE(*,*)' element read ', i,DBLE(i)/DBLE(ipos)
  return
end subroutine trsmReadMat

!================================================================================

function LowerPos(irow,icol) result(ipos)
  ! Written by R. Pong-Wong
  implicit none
  integer, intent(in) :: icol, irow
  integer :: ipos
  if(icol <= irow) then             ! position enquired in lower diagonal
     ipos=(irow-1)*(irow)/2 + icol
  else
     ipos=(icol-1)*(icol)/2 + irow   ! position enquired in upper diagonal swap row and columns
  endif
  return
end function LowerPos

