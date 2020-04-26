subroutine updatetheta(nvar, AI, rhs, theta, verbose)
  implicit none
  integer, intent(in)                                                 :: nvar
  logical, intent(in)                                                 :: verbose
  double precision, dimension(:), intent(in)                          :: AI, rhs
  double precision, dimension(:), intent(inout)                       :: theta
  integer                                                             :: n, info
  integer, dimension(:), allocatable, save                            :: ipiv

  if (verbose) write(6, *) "  In the subroutine solve"
  n = nvar + 1
  if (.not.allocated(ipiv)) allocate(ipiv(n))
  call dspsv('u', n, 1, AI, ipiv, rhs, n, info)
  if (info.ne.0) then
     write(6, *) "  error in solving AI*x = rhs"
     if (info.gt.0) then
        write(6, '(2x,a3,i1,a1,i1,a9)') "AI(", info, ",", info, ") is zero"
     else 
        write(6, *) "  AI had illegal value"
     end if
     stop 1
  else
     if (verbose) write(6, *) "  DSPSV finished solving AI*x=rhs"
  end if

  theta(1 : (nvar +1)) = theta(1 : (nvar + 1)) + rhs(1 : (nvar + 1))

  if (verbose) write(6, *) "  solve returned successfully"
end subroutine updatetheta
