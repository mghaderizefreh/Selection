subroutine updatetheta(nvar, AI, rhs, theta, verbose)
  use constants
  implicit none
  integer, intent(in)                             :: nvar
  logical, intent(in)                             :: verbose
  double precision, dimension(:), intent(in)      :: AI, rhs
  double precision, dimension(:), intent(inout)   :: theta
  integer                                         :: n, info
  integer, dimension(:), allocatable, save        :: ipiv

  external                                        :: dspsv

  if (verbose) write(STDOUT, *) "  In the subroutine solve"
  n = nvar + 1
  if (.not.allocated(ipiv)) allocate(ipiv(n))
  call dspsv('u', n, 1, AI, ipiv, rhs, n, info)
  if (info.ne.0) then
     write(STDOUT, *) "  error in solving AI*x = rhs"
     if (info.gt.0) then
        write(STDERR, '(2x,a3,i1,a1,i1,a9)') "AI(", info, ",", info, ") is zero"
     else 
        write(STDERR, *) "  AI had illegal value"
     end if
     stop 2
  else
     if (verbose) write(STDOUT, *) "  DSPSV finished solving AI*x=rhs"
  end if

  theta(1 : (nvar +1)) = theta(1 : (nvar + 1)) + rhs(1 : (nvar + 1))

  if (verbose) write(STDOUT, *) "  solve returned successfully"
end subroutine updatetheta
