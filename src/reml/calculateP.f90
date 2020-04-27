subroutine calculateP(nobs, nfix, Vinv, X, P, det_xt_vinv_x, Vhat, verbose)
  use constants
  implicit none
  logical, intent(in)                                                   :: verbose
  double precision, dimension(:), intent(in)                            :: Vinv
  double precision, dimension(:,:), intent(in)                          :: X
  integer, intent(in)                                                   :: nobs, nfix
  double precision, dimension(:), intent(out)                           :: P
  double precision, intent(out)                                         :: det_xt_vinv_x
  double precision, dimension(:,:), intent(out)                         :: Vhat

  integer                                                               :: info, I
  integer, dimension(:), allocatable, save                              :: ipiv2
  double precision, dimension(:,:), allocatable, save                   :: Vinvfull, mat, temp
  double precision, dimension(:), allocatable, save                     :: vec, work2

  double precision, external                                            :: ddot   !blas function

  if (verbose) write(stdout,*) "  In the subroutine CalculateP"
  if (.not.allocated(Vinvfull)) then
     I = nfix * nfix
     allocate(Vinvfull(nobs,nobs), ipiv2(nfix), work2(I), temp(nfix,nobs))
  end if

  if (verbose) write(stdout,'(3x,a6,i1,a34)') "using ", nfix, "-dimensional matrix for X"
  if (.not.allocated(mat)) allocate(mat(nfix, nfix), vec((nfix+1)*nfix/2))
  call dunpack('u', nobs, Vinv, Vinvfull, nobs, info)
  if (verbose) write(stdout, *) "  info after DUNPACK", info
  if (info .ne. 0) then
     write(stderr,*) "  error. DUNPACK failed for some reason. Error inside calculateP"
     stop 1
  else
     if (verbose) write(stdout, *) "  DUNPACK finished unpacking vinv to Vinvfull"
  end if

  call dgemm('t', 'n', nfix, nobs, nobs, 1.d0, X, nobs, VinvFull, nobs, 0.d0, temp, nfix)
  if (verbose) write(stdout, *) "  DGEMM finished calculating temp ( = X' * Vinv )"

  call dgemm('n', 'n', nfix, nfix, nobs, 1.d0, temp, nfix, X, nobs, 0.d0, mat, nfix)
  if (verbose) write(stdout, *) "  DGEMM finished calculating mat (= X'Vinv * X)"

  call dtrttp('u', nfix, mat, nfix, vec, info)
  if (verbose) write(stdout, *) "  info after DTRTTP", info
  if (info .ne. 0) then
     write(stderr, *) "  error. DTRTTP failed to pack mat for some reason. Error inside calculateP"
     stop 1
  else
     if (verbose) write(stdout, *) "  DTRTTP finished packing mat to vec"
  end if

  call detinv(nfix, vec, det_xt_vinv_x, ipiv2, work2,  verbose)
  if (verbose) write(stdout, *) "  DETINV finished calculating inverse and determinant of vec (=[X'VinvX]inv)"

  call dunpack('u', nfix, vec, mat, nfix, info)
  if (info .ne. 0) then
     write(stdout, *) "  error. DUNPACK failed to unpack vec for some reasons. Error inside calculateP"
     stop 1
  else
     if (verbose) write(stdout, *) "  DUNPACK finished unpacking vec to mat"
  end if

  call dgemm('n', 'n', nfix, nobs, nfix, 1.d0, mat, nfix, temp, nfix, 0.d0, Vhat, nfix)
  if (verbose) write(stdout, *) "  DGEMM finished calculating Vhat (= [X'VinvX]inv * X'Vinv )"

  call dgemm('t', 'n', nobs, nobs, nfix, -1.d0, temp, nfix, Vhat, nfix, 0.d0, Vinvfull, nobs)
  if (verbose) write(stdout, *) "  DGEMM finished calculating Vinvfull (= - VinvX[X'VinvX]inv * X'Vinv )"

  call dtrttp('u', nobs, Vinvfull, nobs, P, info)
  if (verbose) write(stdout, *) "  info after DTRTTP", info
  if (info .ne. 0) then
     write(stderr, *) "  error. DTRTTP failed to pack Vinvfull to P for some reason. Error inside calculateP"
     stop 1
  else
     if (verbose) write(stdout, *) "  Vinvfull packed to P"
  end if
  I = nobs * (nobs + 1) / 2
  call daxpy(I, 1.d0, vinv, 1, P, 1)
  if (verbose) write(stdout, *) "  DAXPY finished calculating P (=vinv + P_old)"

  if (verbose) write(stdout,*) "  calculateP returned successfully"
end subroutine calculateP

