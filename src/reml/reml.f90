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
     theta, verbose, ipiv, Py, P, V, Vhat, temp, ifail, EmIterations, maxIters)

  use constants, only: KINDR, JArr, alloc1D, STDOUT, STDERR, ZERO
  use global_module, only: detInv
  use blup_module
  implicit none
  !! ================ variable definitions  ================ !!
  logical, intent(in) :: verbose
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
  integer, intent(inout) :: ifail
  integer, intent(in), optional :: EmIterations, maxIters

  type(Jarr), dimension(nvar) :: theZGZ
  type(JArr), dimension((nvar+1)) :: f
  real(KINDR), dimension((nvar+1)) :: oldtheta
  real(KINDR), dimension(((nvar+1)*(nvar+2)/2)) :: AI
  real(KINDR), dimension((nvar+1)) :: rhs
  real(KINDR), dimension(:), allocatable :: work
  real(KINDR) :: logl, oldlogl, epsilon = 1.d-6
  real(KINDR) :: val1, val2, val3
  real(KINDR) :: detV, det_xt_vinv_x, yPy
  integer :: i, j, emIteration
  integer :: iter, maxIter
  logical :: checkForLogL, exceptionPass
  real(KINDR), external :: dnrm2, ddot, dasum
  external :: updateTheta, calculateRHS, calculateLogL, calculateAIMatrix
  !! ================ No defintion after this line ================ !!
  I = nobs * nobs
  call alloc1D(work, I, "work","reml")
  checkForLogL = .true. ! default is true
  exceptionPass = .false. ! default is false (because exceptionPass /= checkForLogL)
  ifail = 0
  ! f is going to contain P*ZGZ_i
  do i = 1, nvar + 1
     call alloc1D(f(i)%array, nobs, "f(i)%array", "reml")
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
        call alloc1D(theZGZ(j)%array, nobs, "theZGZ(3)%array", "reml")
     else        
        call alloc1D(theZGZ(j)%array, i, "theZGZ(j)%array", "reml")
     end if
  end do
  if (nvar .eq. 3) then
     call getMatricesUncorrelated(verbose, nobs, nfix, maxid, X, Gmatrix, id, &
          theZGZ(1)%array, theZGZ(2)%array, theZGZ(3)%array)
  elseif (nvar .eq. 4) then
     call getMatricesCorrelated(verbose, nobs, nfix, maxid, X, Gmatrix, id, &
          theZGZ(1)%array, theZGZ(2)%array, theZGZ(3)%array, &
          theZGZ(4)%array)
  else
     call getMatrices(verbose, nobs, nfix, maxid, X, Gmatrix, id, theZGZ(1)%array)
  end if

69 format(a12, i3)
70 format(1x, a9)
71 format(3x, "-", a, ":", 1x, g24.15)
72 format("  warning: error in reml")
  write(STDOUT, 69) "iteration: " ,0
  write(STDOUT, 70) " theta_0:"
  if (nvar .eq. 1) then
     write(STDOUT, '(a24,1x,a24)') "A", "E"
     write(STDOUT, *) theta(1:2)
  else
     write(STDOUT, '(a24,1x,a24)',advance = 'no') "A_slope", "A_intercept"
     if (nvar == 4) write(STDOUT, '(1x,a24)', advance = 'no') "covariance"
     write(STDOUT, '(2(1x,a24))') "E_slope", "E_intercept"
     write(STDOUT, '(g24.15,1x,g24.15)', advance = 'no') theta(1:2)
     if (nvar == 4) write(STDOUT, '(1x,g24.15)', advance = 'no') theta(4)
     write(STDOUT, '(2(1x,g24.15))') theta(3), theta(nvar + 1)
  end if

  iter = 0
  do
     if (iter > 0) oldlogl = logl

     iter = iter + 1
     write(STDOUT, *)
     write(STDOUT, 69) "iteration: ", iter

     call calculateV(nobs, nvar, theta, theZGZ, V, verbose)
     if (verbose) write(STDOUT, *) " V is calculated"

     call detInv(nobs, V, detV, ipiv, Py, verbose, ifail) !Py is work array here
     if (ifail /= 0) then
        write(STDOUT, 72) 
        return
     end if
     if (verbose) write(STDOUT, *) " V is replaced by its inverse"

     call calculateP(nobs, nfix, V, X, P, det_xt_vinv_x, Vhat, work, temp,&
         verbose, ifail)
     if (ifail /= 0) then
        write(STDOUT, 72)
        return
     end if
     if (verbose) write(STDOUT, *) " P is calcuated"

     call calculateLogL(nobs, detV,det_xt_vinv_x,P, y, LogL,Py, yPy, verbose)
     if (verbose) write(STDOUT, *) " LogL, Py and yPy are calculated"
     write(STDOUT, '(1x,a8,g24.15)') " LogL = ", logl

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

        call updatetheta(nvar, AI, rhs, theta, verbose, ifail)
        if (ifail /= 0) then
           write(STDOUT, 72)
           return
        end if
        if (verbose) write(STDOUT, *) " theta is updated"
     end if

     write(STDOUT, *) " variance vector:"
     if (nran .eq. 1) then
        write(STDOUT, '(a24,1x,a24)') "A", "E"
        if ((theta(1) * theta(2) < 0).and.(theta(1) < 0)) then
           ! fixing negative genetic variance component (only for nran = 1)
           theta(1) = 0.001 * theta(2)
           exceptionPass = .true. ! switch for not checking LogL in next iteration (on)
           write(STDOUT, '(g24.15,a)', advance = 'no') theta(1), "* "
           write(STDOUT, '(g24.15)') theta(2)
        else
           exceptionPass = .false.! switch for not checking LogL in next iter (off => check)
           write(STDOUT, '(g24.15,1x,g24.15)') theta(1:2)
        end if
     else
        write(STDOUT, '(a24,1x,a24)',advance = 'no') "A_slope", "A_intercept"
        if (nvar == 4) write(STDOUT, '(1x,a24)', advance = 'no') "covariance"
        write(STDOUT, '(2(1x,a24))') "E_slope", "E_intercept"
        write(STDOUT, '(g24.15,1x,g24.15)', advance = 'no') theta(1:2)
        if (nvar == 4) write(STDOUT, '(1x,g24.15)', advance = 'no') theta(4)
        write(STDOUT, '(2(1x,g24.15))') theta(3), theta(nvar + 1)
     end if

     ! finally, if the plan is to logl check 
     if (.not. exceptionPass) then
        ! but this is the first ai iteration after em
        if ( (emiteration > 0) .and. (iter  == (emiteration+1)) ) then
           exceptionPass = .true. ! logl check should be ignored for next iteration
        end if
        !otherwise leave the switch off
     end if
     
     !!!! Iteration completed !!!!

     val1 = dnrm2(nvar + 1, oldtheta, 1) ! l2 norm of oldtheta
     oldtheta(1:(nvar + 1)) = oldtheta(1:(nvar + 1)) - theta(1:(nvar + 1))
     val2 = dnrm2(nvar + 1, oldtheta, 1) / val1 ! l2 norm of delta_theta / l2 norm of oldtheta
     val1 = dasum(nvar + 1, oldtheta, 1) / (nvar + 1) ! l1 norm of delta_theta
     write(STDOUT, '(a, i0)') "Errors for iteration ",iter
     write(STDOUT, 71, advance='no') " variance (l1)", val1 
     write(STDOUT, 71) " variance (l2)", val2
     if (iter > 1) then
      val3 = dabs((LogL - oldLogL)/oldLogL) ! l1 norm of delta_logl
      write(STDOUT, 71) " LogL (relative)", val3
      if (checkForLogL .and. (logL < oldLogL)) then
         write(STDOUT, '(A)') " Invalid variance comp. because LogL decreased"
         ifail = 1
         return
      end if
     end if
     ! next iteration checkForLogL is determined here
     checkForLogL = .not. exceptionPass 

     if ((val1 < sqrt(epsilon)) .or. (val2 < epsilon) .and. (val3 < epsilon)) then
        write(STDOUT, '(a10)') "converged!"
        exit
     elseif (iter > maxiter) then
        if ((val1 < sqrt(epsilon)).or.(val2 < epsilon) .and. (val3 > epsilon)) then
           write(STDOUT, '(a)') " warning: variance vector converged but not LogL"
        elseif ((val1 > sqrt(epsilon)).and.(val2 > epsilon)) then
           write(STDOUT, '(a)') " warning: variance vector did not converge"
           ifail = 1
           return
        end if
     end if
     oldtheta(1 : (nvar + 1)) = theta(1 : (nvar + 1))
  end do
  deallocate(work)

  ! reml may fail however (for nran = 3, because for nran = 1 it is modified)
  if (nran==1) then
     if (theta(2) < ZERO) then
        write(STDOUT, '(a)') "Invalid variance component in Reml"
        ifail = 1
        return
     end if
  elseif (nran==3) then
     if (any(theta(1:3) < ZERO) .or. (theta(nvar + 1) < ZERO) ) then
        write(STDOUT, '(a)') "Invalid variance component in Reml"
        ifail = 1
        return
     end if
  end if
  ifail = 0

end subroutine reml
