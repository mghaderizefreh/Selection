module iteration
contains
  subroutine iterate(nobs, nvar, nfix, theZGZ, y, x, logl, theta, Py, Vhat, iter, emiteration, verbose)
    use constants
    use global_module
    use blup_module
    use reml_module
    implicit none
    logical, intent(in)                                                 :: verbose
    integer, intent(in)                                                 :: nvar, nobs, nfix, emiteration, iter
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
    type (ArrOfArr), dimension(:), allocatable, save                    :: f

    if (.not.allocated(P)) then
       I = nobs * (nobs + 1) / 2
       allocate(P(I),V(I))
       I = (nvar + 1) * (nvar + 2) / 2
       allocate(AI(I), rhs(nvar + 1))
       I = nobs * nobs
       allocate(work(I),ipiv(nobs))
       ! f is going to contain P*ZGZ_i
       allocate(f(nvar+1))
       do i = 1, nvar + 1
          allocate(f(i)%array(nobs))
       end do
    end if

    call calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose)
    if (verbose) write(stdout, *) " V is calculated"

    call detInv(nobs, V, detV, ipiv, work, verbose)
    if (verbose) write(stdout, *) " V is replaced by its inverse"

    call calculateP(nobs, nfix, V, X, P, det_xt_vinv_x, Vhat, verbose)
    if (verbose) write(stdout, *) " P is calcuated"

    call calculateLogL(nobs, detV, det_xt_vinv_x, P, y, LogL,Py, yPy, verbose)
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
    write(stdout, '(a11,2x,a11)',advance = 'no') "A_slope", "A_intercept"
    if (nvar == 4) write(stdout, '(2x,a11)', advance = 'no') "covariance"
    write(stdout, '(2(2x,a11))') "E_slope", "E_intercept"
    write(stdout, '(f11.7,2x,f11.7)', advance = 'no') theta(1:2)
    if (nvar == 4) write(stdout, '(2x,f11.7)', advance = 'no') theta(4)
    write(stdout, '(2(2x,f11.7))') theta(3), theta(nvar + 1)

  end subroutine iterate

end module iteration
