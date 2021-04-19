! order of matrices, and therefore variances shall be 
!   1     ZsGZs       genetic slope
!   2     ZiGZi       genetic intercept
!   3     ZsZs        environmental slope (stored as diagonal)
!   4     ZiGZs+ZsGZi genetic slope-intercept covariance
!   5     ZsZs        perm. env. effect slope (SAME AS NUMBER 3: but counted)
!   6     ZiZi        perm. env. effect intercept
!   7     ZsZi+ZiZs   perm. env. effect slope-intercept covariance
! (last)  Identity    env. intercept (NOT COUNTED IN theZGZ and it is LAST one)
subroutine reml(id, X, y, nfix, nobs, maxid, nelement, Gmatrix, nvar, nran,&
     theta, verbose, ipiv, Py, P, V, Vhat, temp, EmIterations, maxIters)

  use constants
  use global_module
  use blup_module
  implicit none
  !! ================ variable definitions  ================ !!
  logical, intent(in)                            :: verbose
  integer, intent(in) :: nobs, nvar, nfix, maxid, nran, nelement
  integer, dimension(:), intent(in) :: id ! real(KINDR) id of animals
  real(KINDR), dimension(:), intent(in) :: y ! phenotypes
  real(KINDR), dimension(:,:), intent(in) :: x ! incid. mat fixed effects
  real(KINDR), dimension(:), intent(inout) :: theta ! main output (variance component)
  real(KINDR), dimension((maxid*(maxid+1)/2)), intent(in) :: Gmatrix
  ! working arrays
  integer, dimension(1:nobs), intent(inout) :: ipiv
  real(KINDR), dimension(1:nobs), intent(inout) :: Py
  real(KINDR), dimension(1:nelement), intent(inout) :: P
  real(KINDR), dimension(1:nelement), intent(inout) :: V
  real(KINDR), dimension(1:nfix,1:nobs), intent(inout) :: Vhat
  real(KINDR), dimension(1:nobs,1:nfix), intent(inout) :: temp
  integer, intent(in), optional :: EmIterations, maxIters

  type(doublePre_Array), dimension(nvar) :: theZGZ
  type(JArr), dimension((nvar+1)) :: f
  real(KINDR), dimension((nvar+1)) :: oldtheta
  real(KINDR), dimension(((nvar+1)*(nvar+2)/2)) :: AI
  real(KINDR), dimension((nvar+1)) :: rhs
  real(KINDR), dimension(:), allocatable :: work
  real(KINDR) :: logl, epsilon = 1.d-6
  real(KINDR) :: val1, val2
  real(KINDR) :: detV, det_xt_vinv_x, yPy
  integer :: i, j, emIteration
  integer :: ifail, iter, maxIter
  real(KINDR), external :: dnrm2, ddot, dasum
  !! ================ No defintion after this line ================ !!
  I = nobs * nobs
  allocate(work(I))
  
  ! f is going to contain P*ZGZ_i
  do i = 1, nvar + 1
     allocate(f(i)%array(nobs))
  end do

  if (.not.present(EmIterations)) then
     EmIteration = 3
  else
     emIteration = emIterations
  end if
  if (.not.present(maxIters)) then
     maxIter = emIteration + 8
  else
     maxIter = maxIters
  end if

  ! order of variances: (As, Ai, Es, Cov, Ei) 
  oldtheta(1:(nvar + 1)) = theta(1:(nvar + 1))

  ! depending on the given correlation, we may need 3 or 4 ZGZ matrices.
  ! So better to check that because one matrix makes a lot of difference
  ! in using the amount of memory
  if (nran == 3) then
     if (nvar == 3) then
        if (any ( theta < 0 )) then
           write(STDERR, *) "n_var and initial guess not consistent"
           stop 2
        end if
        if (verbose) write(STDOUT, '(2x,a22)') "no correlation assumed"
     elseif (nvar > 3) then
        if (verbose) write(STDOUT, '(2x,a30)') "correlation taken into account"
     end if
  end if

  ! making G* matrices
!  allocate(theZGZ(nvar))
  i = nobs * (nobs + 1) / 2
  do j = 1, nvar 
     if (j .eq. 3) then
        allocate(theZGZ(j)%level(nobs))
     else        
        allocate(theZGZ(j)%level(i))
     end if
  end do
  if (nvar .eq. 3) then
     call getMatricesUncorrelated(verbose, nobs, nfix, maxid, X, Gmatrix, id, &
          theZGZ(1)%level, theZGZ(2)%level, theZGZ(3)%level)
  elseif (nvar .eq. 4) then
     call getMatricesCorrelated(verbose, nobs, nfix, maxid, X, Gmatrix, id, &
          theZGZ(1)%level, theZGZ(2)%level, theZGZ(3)%level, &
          theZGZ(4)%level)
  else
     call getMatrices(verbose, nobs, nfix, maxid, X, Gmatrix, id, theZGZ(1)%level)
  end if

69 format(a12, i3)
70 format(1x, a9)
71 format(1x, a10, 1x, f24.15, a10, 1x, f24.15)

  write(STDOUT, 69) "iteration: " ,0
  write(STDOUT, 70) " theta_0:"
  if (nvar .eq. 1) then
     write(STDOUT, '(a11,2x,a11)') "A", "E"
     write(STDOUT, *) theta(1:2)
  else
     write(STDOUT, '(a11,2x,a11)',advance = 'no') "A_slope", "A_intercept"
     if (nvar == 4) write(STDOUT, '(2x,a11)', advance = 'no') "covariance"
     write(STDOUT, '(2(2x,a11))') "E_slope", "E_intercept"
     write(STDOUT, '(f11.7,2x,f11.7)', advance = 'no') theta(1:2)
     if (nvar == 4) write(STDOUT, '(2x,f11.7)', advance = 'no') theta(4)
     write(STDOUT, '(2(2x,f11.7))') theta(3), theta(nvar + 1)
  end if

  iter = 0
  do 
     iter = iter + 1
     write(STDOUT, *)
     write(STDOUT, 69) "iteration: ", iter

     call calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose)
     if (verbose) write(STDOUT, *) " V is calculated"

     call detInv(nobs, V, detV, ipiv, Py, verbose) !Py is work array here
     if (verbose) write(STDOUT, *) " V is replaced by its inverse"

     call calculateP(nobs, nfix, V, X, P, det_xt_vinv_x, Vhat, work, temp, verbose)
     if (verbose) write(STDOUT, *) " P is calcuated"

     call calculateLogL(nobs, detV,det_xt_vinv_x,P, y, LogL,Py, yPy, verbose)
     if (verbose) write(STDOUT, *) " LogL, Py and yPy are calculated"
     write(STDOUT, '(1x,a8,g25.16)') " LogL = ", logl

     call calculaterhs(nobs, nvar, theZGZ, P, Py, rhs, f, verbose)
     if (verbose) write(STDOUT, *) " Right hand side is calculated"

     if (iter  <= emiteration) then
        if (verbose) write(STDOUT, *) 'em  iteration'
        do i = 1, nvar + 1
           theta(i) = theta(i) + 2 * (theta(i)**2) * rhs(i)  / dble(nobs)
        end do
     else
        if (verbose) write(STDOUT, *) 'ai iteration'
        call calculateAImatrix(nobs, nvar, P, AI, f, verbose)
        if (verbose) write(STDOUT, *) " AI matrix is calcuated"

        call updatetheta(nvar, AI, rhs, theta, verbose)
        if (verbose) write(STDOUT, *) " theta is updated"
     end if
     write(STDOUT, *) " variance vector:"
     if (nvar .eq. 1) then
        write(STDOUT, '(a11,2x,a11)') "A", "E"
        write(STDOUT, *) theta(1:2)
     else
        write(STDOUT, '(a11,2x,a11)',advance = 'no') "A_slope", "A_intercept"
        if (nvar == 4) write(STDOUT, '(2x,a11)', advance = 'no') "covariance"
        write(STDOUT, '(2(2x,a11))') "E_slope", "E_intercept"
        write(STDOUT, '(f11.7,2x,f11.7)', advance = 'no') theta(1:2)
        if (nvar == 4) write(STDOUT, '(2x,f11.7)', advance = 'no') theta(4)
        write(STDOUT, '(2(2x,f11.7))') theta(3), theta(nvar + 1)
     end if

!!!! Iteration completed !!!!

     val1 = dnrm2(nvar + 1, oldtheta, 1)
     oldtheta(1:(nvar + 1)) = oldtheta(1:(nvar + 1)) - theta(1:(nvar + 1))
     val2 = dnrm2(nvar + 1, oldtheta, 1) / val1
     val1 = dasum(nvar + 1, oldtheta, 1) / (nvar + 1)
     write(STDOUT, *) "Errors (iter): ",iter
     write(STDOUT, 71, advance='no') " l1 error:", val1 ,"; l2 error:", val2

     write(STDOUT, *) 

     if ((val1 < sqrt(epsilon)) .or. (val2 < epsilon)) then
        write(STDOUT, '(a10)') "converged!"
        exit
     elseif (iter > maxiter) then
        write(STDERR, '(a)') "reml did not converge"
        stop 1
     end if
     oldtheta(1 : (nvar + 1)) = theta(1 : (nvar + 1))
  end do
  deallocate(work)
end subroutine reml
