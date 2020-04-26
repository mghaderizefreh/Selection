subroutine calculateLogL(nobs, detV, det_xt_vinv_x, P, y, LogL,Py, yPy, verbose)
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs
  double precision, intent(in)                                        :: detV, det_xt_vinv_x
  double precision, dimension(:), intent(in)                          :: P, y
  double precision, dimension(:), intent(out)                         :: Py
  double precision, intent(out)                                       :: yPy, LogL
  double precision, external                                          :: ddot

  if (verbose) write(6, *) "  In the subroutine calculateLogL"

  call dspmv('u', nobs, 1.d0, P, y, 1, 0.d0, Py, 1)
  if (verbose) write(6, *) "  DSPMV finished calculating Py (=P * y)"

  yPy = ddot(nobs, y, 1, Py, 1)
  if (verbose) write(6, *) "  DDOT finished calculating yPy"

  LogL = -.5 * (detV + det_xt_vinv_x + yPy)
  if (verbose) write(6, *) "  calculateLogL returned successfully"

end subroutine calculateLogL

