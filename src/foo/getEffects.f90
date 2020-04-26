!This subroutine calculates the fixed and random effects based on the given variances (theta). 
!It then creates (or overwrites) three files with these information.
! In writing the random effects First comes slope effects; then intercept effects; then individual slope residual
!written by Masoud Ghaderi
!last modified 25 March 2020
subroutine getEffects(nobs, maxid, nfix, nvar, fixeffFile, raneffFile, varFile, theta, AmatFile, Vhat, Py, y, X, id, verbose)
  use constants
  implicit none
  logical                                                             :: verbose
  integer, intent(in)                                                 :: nobs, nfix, nvar, maxid
  character(len=*)                                                    :: fixEffFile, ranEffFile, varFile, AmatFile
  integer, dimension(:), intent(in)                                   :: id
  double precision, dimension(:), intent(in)                          :: theta, Py, y
  double precision, dimension(:,:), intent(in)                        :: Vhat, X

  integer                                                             :: iunfix, iunvar, iunran, iunAmat
  integer                                                             :: i, j
  character(len=60)                                                   :: formato
  double precision                                                    :: val1, s1, s2
  double precision, allocatable, dimension(:)                         :: fixeff
  type (doublePre_Array), dimension(:), allocatable, target           :: theZPy, raneff

  ! allocation
  allocate(theZPy(3), fixeff(nfix), raneff(3))
  allocate(theZPy(1)%level(maxid), raneff(1)%level(maxid)) ! slope effect (genetic)
  allocate(theZPy(2)%level(maxid), raneff(2)%level(maxid)) ! intercept effect (genetic)
  allocate(theZPy(3)%level(nobs), raneff(3)%level(nobs))   ! environment slope effect (diagonal)

  ! initialisation
  do i = 1, 3
     theZPy(i)%level(:) = 0.d0
     raneff(i)%level(:) = 0.d0
  end do


  ! fixed effects
  call dgemm('n', 'n', nfix, 1, nobs, 1.d0, Vhat, nfix, y, nobs, 0.d0, fixeff, nfix)
  if (verbose) write(6, *) 'fixed effects: ' , fixeff(1 : nfix)
268 format(a2, i1, a22)
  write(formato, 268) "((", (nfix-1), "(g24.15, 1x), g24.15))"
  open(newUnit = iunfix, file = fixEffFile)
  write(iunfix, '(2a24)') "slope", "intercept"
  write(iunfix, trim(formato)) fixeff(1 : nfix)
  close(iunfix)

  ! variances
!269 format(a10, i1, a20)
!  write(formato, 269) "(1x, a11, ", nvar, "(g24.15, 1x),g24.15)"
!  write(6, trim(formato)) " variance: ", theta(1 : (nvar + 1))
!  if (verbose) write(6, *) 
  open(newUnit = iunvar, file = varFile)
270 format(a2, i1, a21)
273 format(a1, i1, a4)
  write(formato, 273) "(", nvar, "a24)"
  if (nvar .eq. 4) then
     write(iunvar, formato) "var_A_slope","var_A_intercept", &
          "corr(A_int,A_slope)", "var_E_slope", "var_E_intercept"
     write(formato, 270) "(", nvar, "(f24.15, 1x), f24.15)"
     write(iunvar, formato)  theta(1:2), theta(4) / sqrt(theta(1) &
          * theta(2)) , theta(3), theta(5)
  elseif (nvar .eq. 3) then
     write(iunvar, formato) "var_A_slope","var_A_intercept", &
          "var_E_slope", "var_E_intercept"
     write(formato, 270) "(", nvar, "(f24.15, 1x), f24.15)"
     write(iunvar, formato)  theta(1:4)
  end if
  close(iunvar)

  ! random effects
  open(newUnit = iunAmat, file = amatFile)
  do i = 1, nobs
     j = id(i)
     val1 = Py(i) * X(i,1)
     theZPy(1)%level(j) = theZPy(1)%level(j) + val1
     theZPy(2)%level(j) = theZPy(2)%level(j) + Py(i)
     raneff(3)%level(i) = X(i,1) * Py(i)
  end do
  
  if (nvar == 4) then
     s1 = theta(4) / theta(1)
     s2 = theta(4) / theta(2)
     theZPy(1)%level(:) = theZPy(1)%level(:) + theZPy(2)%level(:) * s1
     theZPy(2)%level(:) = theZPy(2)%level(:) + theZPy(1)%level(:) * s2
  end if

  do 
     read(iunAmat, *, end = 271) i, j, val1
     raneff(1)%level(i) = raneff(1)%level(i) + theZPy(1)%level(j) * val1
     raneff(2)%level(i) = raneff(2)%level(i) + theZPy(2)%level(j) * val1
     if (i .ne. j) then
        raneff(1)%level(j) = raneff(1)%level(j) + theZPy(1)%level(i) * val1
        raneff(2)%level(j) = raneff(2)%level(j) + theZPy(2)%level(i) * val1
     else

     end if
  end do


271 continue
  close(iunAmat)
272 format(i12, 1x, g24.15)

  open(newUnit = iunran, file = raneffFile)
  do i = 1, maxid
     write(iunran, 272) i, raneff(1)%level(i)
  end do
  do i = 1, maxid
     write(iunran, 272) i, raneff(2)%level(i)
  end do
  do i = 1, nobs
     write(iunran, 272) i, raneff(3)%level(i)
  end do
  close(iunran)


end subroutine getEffects
