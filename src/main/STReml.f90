module ST_Reml
contains
  subroutine STReml(id, X, y, nfix, nobs, maxid, Gmatrix, nvar, theta, &
       fixEffects, ranEffects, verbose, EmIterations, maxIters)

    use constants
    use global_module
    use blup_module
    use reml_module
    implicit none
    !! ================ variable definitions  ================ !!
    logical, intent(in)                            :: verbose
    integer, intent(in)                            :: nobs, nvar, nfix, maxid
    integer, dimension(:), intent(in)              :: id
    double precision, dimension(:), intent(in)     :: y
    double precision, dimension(:,:), intent(in)   :: x
    double precision, dimension(:), intent(inout)  :: theta
    double precision, dimension(:), intent(in)     :: Gmatrix
    integer, intent(in), optional                  :: EmIterations, maxIters

    double precision, dimension(:), intent(out)    :: fixEffects
    type(doublePre_Array),dimension(:),intent(out) :: ranEffects

    type(doublePre_Array),dimension(:),allocatable :: theZGZ
    type(ArrOfArr), dimension(:), allocatable      :: f
    double precision, dimension(:), allocatable    :: oldtheta, Py, P, V, AI, rhs, work
    double precision                               :: logl, epsilon = 1.d-6
    double precision                               :: val1, val2
    double precision                               :: detV, det_xt_vinv_x, yPy
    integer, dimension(:), allocatable             :: ipiv
    double precision, dimension(:,:), allocatable  :: Vhat
    integer                                        :: i, k, j, emIteration
    integer                                        :: ifail, iter, maxIter
    double precision, external                     :: dnrm2, ddot, dasum
    !! ================ No defintion after this line ================ !!

    allocate(oldtheta(nvar + 1))
    allocate(Vhat(nfix,nobs), Py(nobs)) 
    I = nobs * (nobs + 1) / 2
    allocate(P(I),V(I))
    I = (nvar + 1) * (nvar + 2) / 2
    allocate(AI(I), rhs(nvar + 1))
    I = nobs * nobs
    allocate(work(I),ipiv(nobs))
    allocate(f(nvar+1))
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

    oldtheta(1 : (nvar + 1)) = theta(1 : (nvar + 1))


    allocate(theZGZ(nvar))
    i = nobs * (nobs + 1) / 2
    do j = 1, nvar
       allocate(theZGZ(j)%level(i))
    end do

    call getMatrices(verbose, nobs, X, Gmatrix, id, theZGZ(1)%level)

69  format(a12, i3)
70  format(1x, a9)
71  format(1x, a10, 1x, f24.15, a10, 1x, f24.15)

    write(stdout, 69) "iteration: " ,0
    write(stdout, 70) " theta_0:"
    write(stdout, '(a11,2x,a11)') "A", "E"
    write(stdout, *) theta(1:2)

    iter = 0
    do 
       iter = iter + 1
       write(stdout, *)
       write(stdout, 69) "iteration: ", iter
       theta(:) = oldtheta(:)
       write(6, *) 'Theta::: ', theta
       
       call calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose)
       if (verbose) write(stdout, *) " V is calculated"

       call detInv(nobs, V, detV, ipiv, work, verbose)
       if (verbose) write(stdout, *) " V is replaced by its inverse"

       call calculateP(nobs, nfix, V, X, P, det_xt_vinv_x, Vhat, verbose)
       if (verbose) write(stdout, *) " P is calcuated"

       call calculateLogL(nobs, detV,det_xt_vinv_x,P, y, LogL,Py, yPy, verbose)
       if (verbose) write(stdout, *) " LogL is calculated"
       write(stdout, '(1x,a8,g25.16)') " LogL = ", logl

       call calculaterhs(nobs, nvar, theZGZ, P, Py, rhs, f, verbose)
       if (verbose) write(stdout, *) " Right hand side is calculated"

       if (iter  <= emiteration) then
          if (verbose) write(stdout, *) 'em  iteration'
          do i = 1, nvar + 1
             theta(i) = theta(i) + 2 * (theta(i)**2) * rhs(i)  / dble(nobs)
          end do
       else
          if (verbose) write(stdout, *) 'ai iteration'
          call calculateAImatrix(nobs, nvar, P, AI, f, verbose)
          if (verbose) write(stdout, *) " AI matrix is calcuated"

          call updatetheta(nvar, AI, rhs, theta, verbose)
          if (verbose) write(stdout, *) " theta is updated"
       end if

       write(stdout, *) " variance vector:"
       write(stdout, '(a11,2x,a11)') "A", "E"
       write(stdout, *) theta(1:2)

       !!!! Iteration completed !!!!

       val1 = dnrm2(nvar + 1, oldtheta, 1)
       oldtheta(1 : (nvar + 1)) = oldtheta(1 : (nvar + 1)) - theta(1 : (nvar + 1))
       val2 = dnrm2(nvar + 1, oldtheta, 1) / val1
       val1 = dasum(nvar + 1, oldtheta, 1) / (nvar + 1)
       write(stdout, *) "Errors (iter): ",iter
       write(stdout, 71, advance='no') " l1 error:", val1 ,"; l2 error:", val2

       write(stdout, *) 

       if ((val1 < sqrt(epsilon)) .or. (val2 < epsilon)) then
          write(stdout, '(a10)') "converged!"
          exit
       elseif (iter > maxiter) then
          write(stdout, '(a16)') "did not converge"
          stop 1
       end if
       oldtheta(1 : (nvar + 1)) = theta(1 : (nvar + 1))
    end do

  end subroutine STReml
end module ST_Reml

