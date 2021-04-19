subroutine calculateP(nobs, nfix, Vinv, X, P, det_xt_vinv_x, Vhat, Vinvfull,&
     temp, verbose)
  use constants
  use global_module
  implicit none
  logical, intent(in) :: verbose
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(in) :: Vinv
  real(KINDR), dimension(1:nobs,1:nfix), intent(in) :: X
  integer, intent(in) :: nobs, nfix
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(out) :: P
  real(KINDR), intent(out) :: det_xt_vinv_x
  real(KINDR), dimension(1:nfix,1:nobs), intent(inout) :: Vhat
  real(KINDR), dimension(1:nobs,1:nobs), intent(inout) :: Vinvfull
  real(KINDR), dimension(1:nobs,1:nfix), intent(inout) :: temp
  integer :: info, I
  integer, dimension(nfix) :: ipiv2

  real(KINDR), dimension(nfix,nfix) :: mat
  real(KINDR), dimension(1:((nfix+1)*nfix/2)) :: vec
  real(KINDR), dimension((nfix * nfix)) :: work2

  real(KINDR) :: val0, val1
  external :: dgemm, dtrttp, daxpy

  if (verbose) write(STDOUT,*) "  In the subroutine CalculateP"

  if (verbose) write(STDOUT,'(3x,a6,i2,a34)') "using ", nfix, "-dimensional matrix for X"
  call dunpack('u', nobs, Vinv, Vinvfull, nobs, info)
  if (verbose) write(STDOUT, *) "  info after DUNPACK", info
  if (info .ne. 0) then
     write(STDERR,*) "  error. DUNPACK failed for some reason. Error inside calculateP"
     stop 2
  else
     if (verbose) write(STDOUT, *) "  DUNPACK finished unpacking vinv to Vinvfull"
  end if

  val0 = ZERO
  val1 = ONE
  call dgemm('n', 'n', nobs, nfix, nobs, val1, VinvFull, nobs, X, nobs, val0, temp, nobs)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating temp ( = Vinv * X )"
!  call dgemm('t', 'n', nfix, nobs, nobs, val1, X, nobs, VinvFull, nobs, val0, temp, nfix)
!  if (verbose) write(STDOUT, *) "  DGEMM finished calculating temp ( = X' * Vinv )"

  call dgemm('t', 'n', nfix, nfix, nobs, val1, X, nobs, temp, nobs, val0, mat, nfix)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating mat (= X' * VinvX )"
!  call dgemm('n', 'n', nfix, nfix, nobs, val1, temp, nfix, X, nobs, val0, mat, nfix)
!  if (verbose) write(STDOUT, *) "  DGEMM finished calculating mat (= X'Vinv * X)"

  call dtrttp('u', nfix, mat, nfix, vec, info)
  if (verbose) write(STDOUT, *) "  info after DTRTTP", info
  if (info .ne. 0) then
     write(STDERR, *) "  error. DTRTTP failed to pack mat for some reason. Error inside calculateP"
     stop 2
  else
     if (verbose) write(STDOUT, *) "  DTRTTP finished packing mat to vec"
  end if

  call detinv(nfix, vec, det_xt_vinv_x, ipiv2, work2,  verbose)
  if (verbose) write(STDOUT, *) "  DETINV finished calculating inverse and determinant of vec (=[X'VinvX]inv)"

  call dunpack('u', nfix, vec, mat, nfix, info)
  if (info .ne. 0) then
     write(STDOUT, *) "  error. DUNPACK failed to unpack vec for some reasons. Error inside calculateP"
     stop 2
  else
     if (verbose) write(STDOUT, *) "  DUNPACK finished unpacking vec to mat"
  end if

  call dgemm('n', 't', nfix, nobs, nfix, val1, mat, nfix, temp, nobs, val0, Vhat, nfix)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating Vhat (= [X'VinvX]inv * (VinvX)' )"

!  call dgemm('n', 'n', nfix, nobs, nfix, val1, mat, nfix, temp, nfix, val0, Vhat, nfix)
!  if (verbose) write(STDOUT, *) "  DGEMM finished calculating Vhat (= [X'VinvX]inv * X'Vinv )"

  val1 = -ONE
  call dgemm('n', 'n', nobs, nobs, nfix, val1, temp, nobs, Vhat, nfix, val0, Vinvfull, nobs)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating Vinvfull (= - VinvX * ([X'VinvX]inv*X'Vinv))"
!  call dgemm('t', 'n', nobs, nobs, nfix, val1, temp, nfix, Vhat, nfix, val0, Vinvfull, nobs)
!  if (verbose) write(STDOUT, *) "  DGEMM finished calculating Vinvfull (= - VinvX * ([X'VinvX]inv*X'Vinv))"

  call dtrttp('u', nobs, Vinvfull, nobs, P, info)
  if (verbose) write(STDOUT, *) "  info after DTRTTP", info
  if (info .ne. 0) then
     write(STDERR, *) "  error. DTRTTP failed to pack Vinvfull to P for some reason. Error inside calculateP"
     stop 2
  else
     if (verbose) write(STDOUT, *) "  Vinvfull packed to P"
  end if
  I = nobs * (nobs + 1) / 2
  val1 = ONE
  info = 1
  call daxpy(I, val1, vinv, info, P, info)
  if (verbose) write(STDOUT, *) "  DAXPY finished calculating P (=vinv + P_old)"

  if (verbose) write(STDOUT,*) "  calculateP returned successfully"
end subroutine calculateP

