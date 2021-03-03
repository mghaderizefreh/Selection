subroutine calculateP(nobs, nfix, Vinv, X, P, det_xt_vinv_x, Vhat, verbose)
  use constants
  use global_module
  implicit none
  logical, intent(in) :: verbose
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(in) :: Vinv
  real(KINDR), dimension(1:nobs,1:nfix), intent(in) :: X
  integer, intent(in) :: nobs, nfix
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(out) :: P
  real(KINDR), intent(out) :: det_xt_vinv_x
  real(KINDR), dimension(1:nfix,1:nobs), intent(out) :: Vhat

  integer :: info, I
  integer, dimension(:), allocatable, save :: ipiv2
  real(KINDR), dimension(:,:), allocatable, save :: Vinvfull, mat, temp
  real(KINDR), dimension(:), allocatable, save :: vec, work2

  external :: dgemm, dtrttp, daxpy 

  if (verbose) write(STDOUT,*) "  In the subroutine CalculateP"
  if (.not.allocated(Vinvfull)) then
     I = nfix * nfix
     allocate(Vinvfull(nobs,nobs), ipiv2(nfix), work2(I), temp(nfix,nobs))
  end if

  if (verbose) write(STDOUT,'(3x,a6,i2,a34)') "using ", nfix, "-dimensional matrix for X"
  if (.not.allocated(mat)) allocate(mat(nfix, nfix), vec((nfix+1)*nfix/2))
  call dunpack('u', nobs, Vinv, Vinvfull, nobs, info)
  if (verbose) write(STDOUT, *) "  info after DUNPACK", info
  if (info .ne. 0) then
     write(STDERR,*) "  error. DUNPACK failed for some reason. Error inside calculateP"
     stop 2
  else
     if (verbose) write(STDOUT, *) "  DUNPACK finished unpacking vinv to Vinvfull"
  end if

  call dgemm('t', 'n', nfix, nobs, nobs, 1.d0, X, nobs, VinvFull, nobs, 0.d0, temp, nfix)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating temp ( = X' * Vinv )"

  call dgemm('n', 'n', nfix, nfix, nobs, 1.d0, temp, nfix, X, nobs, 0.d0, mat, nfix)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating mat (= X'Vinv * X)"

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

  call dgemm('n', 'n', nfix, nobs, nfix, 1.d0, mat, nfix, temp, nfix, 0.d0, Vhat, nfix)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating Vhat (= [X'VinvX]inv * X'Vinv )"

  call dgemm('t', 'n', nobs, nobs, nfix, -1.d0, temp, nfix, Vhat, nfix, 0.d0, Vinvfull, nobs)
  if (verbose) write(STDOUT, *) "  DGEMM finished calculating Vinvfull (= - VinvX[X'VinvX]inv * X'Vinv )"

  call dtrttp('u', nobs, Vinvfull, nobs, P, info)
  if (verbose) write(STDOUT, *) "  info after DTRTTP", info
  if (info .ne. 0) then
     write(STDERR, *) "  error. DTRTTP failed to pack Vinvfull to P for some reason. Error inside calculateP"
     stop 2
  else
     if (verbose) write(STDOUT, *) "  Vinvfull packed to P"
  end if
  I = nobs * (nobs + 1) / 2
  call daxpy(I, 1.d0, vinv, 1, P, 1)
  if (verbose) write(STDOUT, *) "  DAXPY finished calculating P (=vinv + P_old)"

  if (verbose) write(STDOUT,*) "  calculateP returned successfully"
end subroutine calculateP

