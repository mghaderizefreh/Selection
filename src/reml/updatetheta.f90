subroutine updatetheta(nvar, AI, rhs, theta, verbose)
  use constants
  implicit none
  integer, intent(in)                                                 :: nvar
  logical, intent(in)                                                 :: verbose
  double precision, dimension(:), intent(in)                          :: AI, rhs
  double precision, dimension(:), intent(inout)                       :: theta
  integer                                                             :: n, info
  integer, dimension(:), allocatable, save                            :: ipiv

  if (verbose) write(stdout, *) "  In the subroutine solve"
  n = nvar + 1
  if (.not.allocated(ipiv)) allocate(ipiv(n))
  call dspsv('u', n, 1, AI, ipiv, rhs, n, info)
  if (info.ne.0) then
     write(stdout, *) "  error in solving AI*x = rhs"
     if (info.gt.0) then
        write(stdout, '(2x,a3,i1,a1,i1,a9)') "AI(", info, ",", info, ") is zero"
     else 
        write(stderr, *) "  AI had illegal value"
     end if
     stop 1
  else
     if (verbose) write(stdout, *) "  DSPSV finished solving AI*x=rhs"
  end if

  theta(1 : (nvar +1)) = theta(1 : (nvar + 1)) + rhs(1 : (nvar + 1))

  if (verbose) write(stdout, *) "  solve returned successfully"
end subroutine updatetheta
