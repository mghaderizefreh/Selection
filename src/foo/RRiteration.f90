subroutine iterate(nobs, nvar, nfix, theZGZ, y, x, logl, theta, Py, Vhat, verbose)
  use constants
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nvar, nobs, nfix
  type (doublePre_array), dimension(:), intent(in)                    :: theZGZ
  double precision, dimension(:,:), intent(in)                        :: x
  double precision, dimension(:), intent(in)                          :: y
  double precision, intent(out)                                       :: logl
  double precision, dimension(:), intent(out)                         :: Py
  double precision, dimension(:,:), intent(out)                       :: Vhat
  double precision, dimension(:), intent(inout)                       :: theta

  integer                                                             :: i, ifail
  double precision                                                    :: detV, det_xt_vinv_x, yPy
  integer, dimension(:), allocatable, save                            :: ipiv
  double precision, dimension(:), allocatable, save                   :: V, P, AI, rhs, work
  double precision, external                                          :: ddot

  if (.not.allocated(P)) then
     I = nobs * (nobs + 1) / 2
     allocate(P(I),V(I))
     I = (nvar + 1) * (nvar + 2) / 2
     allocate(AI(I), rhs(nvar + 1))
     I = nobs * nobs
     allocate(work(I),ipiv(nobs))
  end if

  call calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose)
  if (verbose) write(6, *) " V is calculated"

  call detInv(nobs, V, detV, ipiv, work, verbose)
  if (verbose) write(6, *) " V is replaced by its inverse"

  call calculateP(nobs, nfix, V, X, P, det_xt_vinv_x, Vhat, verbose)
  if (verbose) write(6, *) " P is calcuated"

  call calculateLogL(nobs, detV, det_xt_vinv_x, P, y, LogL,Py, yPy, verbose)
  if (verbose) write(6, *) " LogL is calculated"
  write(6, '(1x,a8,g25.16)') " LogL = ", logl

  call calculateAImatrix(nobs, nvar, theZGZ, P, Py, AI, verbose)
  if (verbose) write(6, *) " AI matrix is calcuated"

  call calculaterhs(nobs, nvar, theZGZ, P, Py, rhs, work, verbose)
  if (verbose) write(6, *) " Right hand side is calculated"

  call updatetheta(nvar, AI, rhs, theta, verbose)
  if (verbose) write(6, *) " theta is updated"

  if (verbose) then
     write(6, *) " new variance vector:"
     write(6, *) "gen. slope", theta(1)
     write(6, *) "gen. inter", theta(2)
     write(6, *) "env. slope", theta(3)
     if (nvar == 4) write(6, *) "covariance", theta(4)
     write(6, *) "env. inter", theta(nvar+1)
  end if
end subroutine iterate

