subroutine calculateP(nobs, nfix, Vinv, X, P, det_xt_vinv_x, Vhat, Vinvfull,&
     temp, verbose, info)
  use constants, only: KINDR, STDOUT, ONE, ZERO
  use global_module, only: dunpack, detinv
  implicit none
  logical, intent(in) :: verbose
  integer, intent(in) :: nobs, nfix
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(in) :: Vinv
  real(KINDR), dimension(1:nobs,1:nfix), intent(in) :: X
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(out) :: P
  real(KINDR), intent(out) :: det_xt_vinv_x
  real(KINDR), dimension(1:nfix,1:nobs), intent(inout) :: Vhat
  real(KINDR), dimension(1:nobs,1:nobs), intent(inout) :: Vinvfull
  real(KINDR), dimension(1:nobs,1:nfix), intent(inout) :: temp
  integer, intent(inout) :: info
  integer :: I
  integer, dimension(nfix) :: ipiv2

  real(KINDR), dimension(nfix,nfix) :: mat
  real(KINDR), dimension(1:((nfix+1)*nfix/2)) :: vec
  real(KINDR), dimension((nfix * nfix)) :: work2

  real(KINDR) :: val0, val1
  external :: dgemm, dtrttp, daxpy

401 format("  warning: error in calculateP")
402 format("  warning: error in ",a," info: ", i5,/,&
         "  warning: error in caculate P")
  if (verbose) write(STDOUT,*) "  In the subroutine CalculateP"

  if (verbose) write(STDOUT,'(3x,a6,i2,a34)') "using ", nfix, "-dimensional matrix for X"
  call dunpack('u', nobs, Vinv, Vinvfull, nobs, info)
  if (info .ne. 0) then
     write(STDOUT, 401)
     return
  end if

  val0 = ZERO
  val1 = ONE
  call dgemm('n', 'n', nobs, nfix, nobs, val1, VinvFull, nobs, X, nobs, val0, temp, nobs)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating temp ( = Vinv * X )"

  call dgemm('t', 'n', nfix, nfix, nobs, val1, X, nobs, temp, nobs, val0, mat, nfix)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating mat (= X' * VinvX )"

  call dtrttp('u', nfix, mat, nfix, vec, info)
  if (info .ne. 0) then
     write(STDOUT, 402) "dtrttp", info
     return
  end if
  if (verbose) write(STDOUT, *) "  info after DTRTTP", info

  call detinv(nfix, vec, det_xt_vinv_x, ipiv2, work2,  verbose, info)
  if (info /= 0) then
     write(STDOUT, 401) 
     return
  end if

  call dunpack('u', nfix, vec, mat, nfix, info)
  if (info /= 0) then
     write(STDOUT, 401) 
     return
  end if
  if (verbose) write(STDOUT, *) "  DUNPACK finished unpacking vec to mat"
  
  call dgemm('n', 't', nfix, nobs, nfix, val1, mat, nfix, temp, nobs, val0, Vhat, nfix)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating Vhat (= [X'VinvX]inv * (VinvX)' )"

  val1 = -ONE
  call dgemm('n', 'n', nobs, nobs, nfix, val1, temp, nobs, Vhat, nfix, val0, Vinvfull, nobs)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating Vinvfull (= - VinvX * ([X'VinvX]inv*X'Vinv))"

  call dtrttp('u', nobs, Vinvfull, nobs, P, info)
  if (info .ne. 0) then
     write(STDOUT, 402) "dtrttp", info
     return
  end if
  if (verbose) write(STDOUT, *) "  info after DTRTTP", info

  I = nobs * (nobs + 1) / 2
  val1 = ONE
  info = 1 ! info is the increment (not success state)
  call daxpy(I, val1, Vinv, info, P, info)
  if (verbose) write(STDOUT, *) "  DAXPY finished calculating P (=vinv + P_old)"
  info = 0 ! because it was set to 1 to reuse a variable
  if (verbose) write(STDOUT,*) "  calculateP returned successfully"
end subroutine calculateP

