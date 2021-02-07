subroutine calculateLogL(nobs, detV, det_xt_vinv_x, P, y, LogL,Py, yPy, verbose)
  use constants
  implicit none
  logical, intent(in) :: verbose
  integer, intent(in) :: nobs
  real(KINDR), intent(in) :: detV, det_xt_vinv_x
  real(KINDR), dimension(1:(nobs*(nobs+1)/2)), intent(in) :: P
  real(KINDR), dimension(1:nobs), intent(in) :: y
  real(KINDR), dimension(1:nobs), intent(out) :: Py
  real(KINDR), intent(out) :: yPy, LogL
  real(KINDR), external :: ddot

  external :: dspmv

  if (verbose) write(STDOUT, *) "  In the subroutine calculateLogL"

  call dspmv('u', nobs, 1.d0, P, y, 1, 0.d0, Py, 1)
  if (verbose) write(STDOUT, *) "  DSPMV finished calculating Py (=P * y)"

  yPy = ddot(nobs, y, 1, Py, 1)
  if (verbose) write(STDOUT, *) "  DDOT finished calculating yPy"

  LogL = -.5 * (detV + det_xt_vinv_x + yPy)
  if (verbose) write(STDOUT, *) "  calculateLogL returned successfully"

end subroutine calculateLogL

