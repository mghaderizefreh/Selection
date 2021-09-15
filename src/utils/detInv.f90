!*
!*CALCULATE THE LOG DETERMINANT OF packed stored symmetric matrix
!* after it has been factorized the Bunch-Kaufman diagonal pivoting method:
!* implemented in the lapcak subroutine  DSPTRF computes 
!*     A = U*D*U**T  or  A = L*D*L**T
!---------------------------------------------------
!  let A = UDU'
!     where U is ls upper triangular unit diagonal 
!           D is a (block) diagonal
!    then det(A) = Det(U)*det(D) * Det(U')
!
!  Det(U)=1 because it is triangular unit diagonal
!  
! since the factorization is done using the Bunch-Kaufman diagonal pivoting method
!  d is block diagonals  which each block if only 1x1 or 2x2 
! 
! so det(D) = det(Blocks(1))*... * det(Blocks(i)) ...*... * ... det(block(i+1) * det(block(n))
!
!    in log scale = sum(log(det(block(1..n)))   .  
!   to avoid problem with negative value being log tranformed 
!   do the sumation in absolute value and count the number of negative pivots
! 
! when block = 2x2 det(block) = b[1,1]*b[2,2]- b[1,2]*b[2,1]
!---------------------------------------------------
!
!
! to be sure the factorization is correct, the input of the subroutine 
! should be calculated using the lapcak subroutine DSPTRF
!!
!====================================================

subroutine dspdrf_Ldet ( BLASuplo, N, AP, IPIV, logDet, signDet, INFO )
  use constants
  implicit none
  character (len = *), intent(in) :: BLASuplo
  integer, intent(in) :: N
  real(KINDR), dimension(1:(N*(N+1)/2)), intent(inout) :: AP
  integer, dimension(1:N), intent(in) :: IPIV
  real(KINDR), intent(out) :: logDet
  integer, intent(out) :: signDet
  integer, intent(out) :: info

  real(KINDR) :: val1

  integer :: k, i, ipos, ipos1

  logDet = ZERO
  signDet = -1
  info = -1
101 format("  warning: error in dspdrf_Ldet, info: ", i5)
  if (BLASuplo(1:1) == 'U' .or. BLASuplo(1:1) == 'u') then

     logDet = ZERO
     k = 0
     ipos1 = 0
     do i = 1, N
        ipos = ipos1 + i
        if (ipiv(i) > 0) then
           if (AP(ipos) == ZERO) then
              info = i
              logDet = ZERO
              signDet = - 1
              write(STDOUT, 101) info
              return
           else
              if (AP(ipos) < 0) k = k + 1
              logDet = logDet + log(abs(AP(ipos)))
           end if
        else if (i > 1) then
           if (ipiv(i) < 0 .AND. ipiv(i - 1) .EQ. ipiv(i)) then
              val1 = AP(ipos1) * AP(ipos) - AP(ipos - 1) * AP(ipos - 1)
              if (val1 == ZERO) then
                 info = i
                 logDet = ZERO
                 signDet = - 1
                 write(STDOUT, 101) info
                 return
              else
                 if (val1 < 0) k = k + 1
                 logDet = logDet + log(abs(val1))
              end if
           end if
        end if
        ipos1 = ipos
     end do
     signDet = k
     !==========================================
     info = 0
  else
     write(STDERR, '(a)') "Error:"
     write(STDERR,*) " so far subroutine has been implemeneted only using&
          &the (blas) Col-Upper stored packed matrix"
     stop 3
  endif

end subroutine dspdrf_Ldet


subroutine detInv(nobs, V, detV, ipiv, work, verbose, info)
  use constants
  implicit none
  logical, intent(in) :: verbose
  integer, intent(in) :: nobs
  integer, intent(inout) :: info
  integer, dimension(1:nobs), intent(inout) :: ipiv
  real(KINDR), dimension(1:nobs), intent(inout) :: work
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(inout) :: V
  real(KINDR), intent(out) :: detV

  integer :: ineg
  external :: dsptri, dsptrf, dspdrf_Ldet

  if (verbose) write(STDOUT, *) "  In the subroutine detinv"
  call dsptrf('u', nobs, V, ipiv, info)
  if (verbose) write (STDOUT, *) "  info after DSPTRF", info
  if ( info /= 0 ) then
     write(STDOUT, 100) "dsptrf", info
     return
  end if

  call dspdrf_Ldet ( 'u', nobs, V, ipiv, detV, ineg, info)
  if (verbose) write(STDOUT, *) "  log detV ineg", detV, ineg
  if ( info /= 0) then
     write(STDOUT, 99) 
     return
  end if

  call dsptri('u', nobs, V, ipiv, work, info)
  if (verbose) write(STDOUT, *) "  info after DSPTRI", info
  if ( info /= 0 ) then
     write(STDOUT, 100) "dsptri", info
     return
  end if

99 format("  warning: error in detinv ")
100 format("  warning: error in ",a ," info:",i5,/,"  warning: error in detinv")
  if (verbose) write(STDOUT, *) "  detinv return successfully"
end subroutine detInv


