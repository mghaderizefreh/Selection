module ST_Blup
contains
  subroutine STBlup(id, X, y, nfix, nobs, maxid, Gmatrix, nvar, theta, &
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
    double precision, dimension(:), allocatable    :: Py, P, V, work
    double precision                               :: detV, det_xt_vinv_x
    integer, dimension(:), allocatable             :: ipiv
    double precision, dimension(:,:), allocatable  :: Vhat
    integer                                        :: i, j
    integer                                        :: ifail
    double precision, external                     :: dnrm2, ddot, dasum
    !! ================ No defintion after this line ================ !!
    allocate(Vhat(nfix,nobs), Py(nobs)) 
    I = nobs * (nobs + 1) / 2
    allocate(P(I),V(I))
    allocate(work(I),ipiv(nobs))

    if (present(EmIterations)) I = EmIterations
    if (present(maxIters)) I = maxIters

    allocate(theZGZ(nvar))
    i = nobs * (nobs + 1) / 2
    do j = 1, nvar
       allocate(theZGZ(j)%level(i))
    end do

    call getMatrices(verbose, nobs, X, Gmatrix, id, theZGZ(1)%level)

    call calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose)
    if (verbose) write(stdout, *) " V is calculated"

    call detInv(nobs, V, detV, ipiv, work, verbose)
    if (verbose) write(stdout, *) " V is replaced by its inverse"

    call calculateP(nobs, nfix, V, X, P, det_xt_vinv_x, Vhat, verbose)
    if (verbose) write(stdout, *) " P is calcuated"

    call dspmv('u', nobs, 1.d0, P, y, 1, 0.d0, Py, 1)
    if (verbose) write(6, *) "  DSPMV finished calculating Py (=P * y)"

    call getEffects(nobs, maxid, nfix, nvar, theta, Gmatrix, Vhat, Py, y, X,&
         id, fixeffects, raneffects, verbose)

  end subroutine STBlup
end module ST_Blup

