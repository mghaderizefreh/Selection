subroutine iterate(nobs, nvar, nfix, theZGZ, y, x, logl, theta, Py, Vhat, verbose)
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


subroutine calculateV(nobs, nvar, theta, theZGZ, ifail, V, verbose) 
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs, nvar
  double precision, dimension(:), intent(in)                          :: theta
  type (doublePre_array), dimension(:), intent(in)                    :: theZGZ
  integer, intent(out)                                                :: ifail
  double precision, dimension(:), intent(out)                         :: V

  integer                                                             :: i, k, ipos, irow, isize
  double precision                                                    :: Val1

  if (verbose) write(6,*) "  In the subroutine CalculateV"
  ifail = 1
  isize = (nobs + 1) * nobs/2
  v(1 : isize) = 0.d0
  val1 = theta(nvar + 1)
  ipos = 0
  k = 3 ! the index of the diagonal matrix ZsZs
  do irow = 1, nobs
     ipos    = ipos + irow   ! the position of this diagonal
     V(ipos) = val1 + theta(k) * theZGZ(k)%level(irow)
  end do

  do i = 1, nvar
     if (i .ne. k) then
        v(1 : isize) = v(1 : isize) + theZGZ(i)%level(1 : isize) * theta(i)
     end if
  end do

  ifail=0
  if (verbose) write(6,*) "  calculateV returned succesfully"
end subroutine calculatev

subroutine detInv(nobs, V, detV, ipiv, work, verbose)
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs
  integer, dimension(:), intent(in)                                   :: ipiv
  double precision, dimension(:), intent(in)                          :: work
  double precision, dimension(:), intent(inout)                       :: V
  double precision, intent(out)                                       :: detV

  integer                                                             :: ifail, ineg, info

  if (verbose) write(6,*) "  In the subroutine detinv"
  call dsptrf('u', nobs, V, ipiv, info)
  if (verbose) write ( *, * ) "  info after DSPTRF", info
  if ( info == 0 ) then
     call dspdrf_Ldet ( 'u', nobs, V, ipiv, detV, ineg, ifail )
     if (verbose) write(6,*) "  log detV ineg", detV, ineg
     call dsptri('u', nobs, V, ipiv, work, info)
     if (verbose) write(6,*) "  info after DSPTRI", info
  else if ( info < 0 ) then
     write(6, *) "  error with input variables (inside detinv)", info
     stop 1
  end if
  if (verbose) write(6,*) "  detinv return successfully"
end subroutine detinv


subroutine calculateP(nobs, nfix, Vinv, X, P, det_xt_vinv_x, Vhat, verbose)
  implicit none
  logical, intent(in)                                                   :: verbose
  double precision, dimension(:), intent(in)                            :: Vinv
  double precision, dimension(:,:), intent(in)                          :: X
  integer, intent(in)                                                   :: nobs, nfix
  double precision, dimension(:), intent(out)                           :: P
  double precision, intent(out)                                         :: det_xt_vinv_x
  double precision, dimension(:,:), intent(out)                         :: Vhat

  integer                                                               :: info, I
  integer, dimension(:), allocatable, save                              :: ipiv2
  double precision, dimension(:,:), allocatable, save                   :: Vinvfull, mat, temp
  double precision, dimension(:), allocatable, save                     :: vec, work2

  double precision, external                                            :: ddot   !blas function

  if (verbose) write(6,*) "  In the subroutine CalculateP"
  if (.not.allocated(Vinvfull)) then
     I = nfix * nfix
     allocate(Vinvfull(nobs,nobs), ipiv2(nfix), work2(I), temp(nfix,nobs))
  end if

  if (verbose) write(6,'(3x,a6,i1,a34)') "using ", nfix, "-dimensional matrix for X"
  if (.not.allocated(mat)) allocate(mat(nfix, nfix), vec((nfix+1)*nfix/2))
  call dunpack('u', nobs, Vinv, Vinvfull, nobs, info)
  if (verbose) write(6, *) "  info after DUNPACK", info
  if (info .ne. 0) then
     write(6,*) "  error. DUNPACK failed for some reason. Error inside calculateP"
     stop 1
  else
     if (verbose) write(6, *) "  DUNPACK finished unpacking vinv to Vinvfull"
  end if

  call dgemm('t', 'n', nfix, nobs, nobs, 1.d0, X, nobs, VinvFull, nobs, 0.d0, temp, nfix)
  if (verbose) write(6, *) "  DGEMM finished calculating temp ( = X' * Vinv )"

  call dgemm('n', 'n', nfix, nfix, nobs, 1.d0, temp, nfix, X, nobs, 0.d0, mat, nfix)
  if (verbose) write(6, *) "  DGEMM finished calculating mat (= X'Vinv * X)"

  call dtrttp('u', nfix, mat, nfix, vec, info)
  if (verbose) write(6, *) "  info after DTRTTP", info
  if (info .ne. 0) then
     write(6, *) "  error. DTRTTP failed to pack mat for some reason. Error inside calculateP"
     stop 1
  else
     if (verbose) write(6, *) "  DTRTTP finished packing mat to vec"
  end if

  call detinv(nfix, vec, det_xt_vinv_x, ipiv2, work2,  verbose)
  if (verbose) write(6, *) "  DETINV finished calculating inverse and determinant of vec (=[X'VinvX]inv)"

  call dunpack('u', nfix, vec, mat, nfix, info)
  if (info .ne. 0) then
     write(6, *) "  error. DUNPACK failed to unpack vec for some reasons. Error inside calculateP"
     stop 1
  else
     if (verbose) write(6, *) "  DUNPACK finished unpacking vec to mat"
  end if

  call dgemm('n', 'n', nfix, nobs, nfix, 1.d0, mat, nfix, temp, nfix, 0.d0, Vhat, nfix)
  if (verbose) write(6, *) "  DGEMM finished calculating Vhat (= [X'VinvX]inv * X'Vinv )"

  call dgemm('t', 'n', nobs, nobs, nfix, -1.d0, temp, nfix, Vhat, nfix, 0.d0, Vinvfull, nobs)
  if (verbose) write(6, *) "  DGEMM finished calculating Vinvfull (= - VinvX[X'VinvX]inv * X'Vinv )"

  call dtrttp('u', nobs, Vinvfull, nobs, P, info)
  if (verbose) write(6, *) "  info after DTRTTP", info
  if (info .ne. 0) then
     write(6, *) "  error. DTRTTP failed to pack Vinvfull to P for some reason. Error inside calculateP"
     stop 1
  else
     if (verbose) write(6, *) "  Vinvfull packed to P"
  end if
  I = nobs * (nobs + 1) / 2
  call daxpy(I, 1.d0, vinv, 1, P, 1)
  if (verbose) write(6, *) "  DAXPY finished calculating P (=vinv + P_old)"

  if (verbose) write(6,*) "  calculateP returned successfully"
end subroutine calculatep


subroutine calculateLogL(nobs, detV, det_xt_vinv_x, P, y, LogL,Py, yPy, verbose)
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs
  double precision, intent(in)                                        :: detV, det_xt_vinv_x
  double precision, dimension(:), intent(in)                          :: P, y
  double precision, dimension(:), intent(out)                         :: Py
  double precision, intent(out)                                       :: yPy, LogL
  double precision, external                                          :: ddot

  if (verbose) write(6, *) "  In the subroutine calculateLogL"

  call dspmv('u', nobs, 1.d0, P, y, 1, 0.d0, Py, 1)
  if (verbose) write(6, *) "  DSPMV finished calculating Py (=P * y)"

  yPy = ddot(nobs, y, 1, Py, 1)
  if (verbose) write(6, *) "  DDOT finished calculating yPy"

  LogL = -.5 * (detV + det_xt_vinv_x + yPy)
  if (verbose) write(6, *) "  calculateLogL returned successfully"

end subroutine calculateLogL

subroutine calculateAImatrix(nobs, nvar, theZGZ, P, Py, AI, verbose)
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs, nvar
  type (doublePre_array), dimension(:), intent(in)                    :: theZGZ
  double precision, dimension(:), intent(in)                          :: P, Py
  double precision, dimension(:), intent(out)                         :: AI
  double precision, external                                          :: ddot
  integer                                                             :: i, k, j
  type ArrOfArr
     double precision, dimension(:), allocatable                      :: array
  end type ArrOfArr
  type (ArrOfArr), dimension(:), allocatable, save                    :: f
  double precision, dimension(:), allocatable, save                   :: temp

  if (verbose) write(6, *) "  In the subroutine calculateAImatrix"
  if (.not. allocated(f)) allocate(f(nvar+1))
  if (.not. allocated(f(nvar+1)%array)) allocate(f(nvar+1)%array(nobs))
  if (.not. allocated(temp))  allocate(temp(nobs))

  f(nvar+1)%array(1:nobs) = Py(1:nobs)
  k = 3 ! index of the diagonal matrix ZsZs
  do i = 1, nvar 
     if (.not. allocated(f(i)%array)) allocate(f(i)%array(nobs))
     if (i .ne. k) then
        call dspmv('u', nobs, 1.d0, theZGZ(i)%level, Py, 1, 0.d0, f(i)%array, 1)
     else
        f(i)%array(1:nobs) = theZGZ(i)%level(1:nobs) * Py(1:nobs)
     end if
     if (verbose) write(6, '(2x,a30, i1, a1)') "DSPMV finished calculating f(", i, ")"
  end do

  k = 1
  do i = 1, (nvar + 1)
     call dspmv('u', nobs, 1.d0, P, f(i)%array, 1, 0.d0, temp, 1)

     if (verbose) write(6, '(2x,a32, i1, a1)') "DSPMV finished calculating P*f(", i, ")"
     do j = 1, i
        AI(k) = 0.5d0 * ddot(nobs, temp, 1, f(j)%array, 1)
        k = k + 1 
     end do
  end do

  if (verbose) write(6, *) "  calculateAImatrix returend successfully"
end subroutine calculateAImatrix


subroutine calculaterhs(nobs, nvar, theZGZ, P, Py, rhs, work, verbose)
  implicit none
  logical, intent(in)                                                 :: verbose
  integer, intent(in)                                                 :: nobs, nvar
  type (doublePre_array), dimension(:), intent(in)                    :: theZGZ
  double precision, dimension(:), intent(in)                          :: P, Py
  double precision, dimension(:), intent(out)                         :: rhs
  double precision, external                                          :: ddot !, traceA, traceAxB, traceAxBdiag
  double precision, dimension(:), intent(inout)                       :: work
  integer                                                             :: i, k

  if (verbose) write(6, *) "  In the subroutine calculateRHS"

  rhs(nvar + 1) = -0.5d0 *(traceA(P, nobs) - ddot(nobs, Py, 1, Py, 1)) 
  k = 3
  do i = 1, nvar
     if (i .ne. k) then
        call dspmv('u', nobs, 1.d0, theZGZ(i)%level, Py, 1, 0.d0, work(1:nobs), 1)
        rhs(i) = -0.5d0 * (traceAxB(P, theZGZ(i)%level, nobs) - ddot(nobs, Py, 1, work, 1))
     else
        work(1 : nobs) = theZGZ(i)%level * Py
        rhs(i) = -0.5d0 * (traceAxBdiag(P, theZGZ(i)%level, nobs) - ddot(nobs, Py, 1, work, 1))
     end if
  end do

  if (verbose) write(6, *) "  calculateRHS returned successfully"
end subroutine calculaterhs



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

function traceA(A, nobs) result(trace)
  ! simply calcuates the trace of A with size nobsxnobs
  implicit none
  double precision, dimension(:), intent(in)                          :: A
  integer, intent(in)                                                 :: nobs
  double precision                                                    :: trace
  integer                                                             :: i, j
  trace = 0.d0
  j = 0
  do i = 1, nobs
     j = j + i
     trace = trace + A(j)
  end do
end function traceA

function traceAxB(A, B, nobs) result(trace)
  implicit none
  ! This function calculates the trace of AxB when A and both
  ! are symmetric and of size nobsxnobs. It uses ddot of BLAS
  double precision, dimension(:), intent(in)                          :: A, B
  integer, intent(in)                                                 :: nobs
  double precision                                                    :: trace
  double precision, external                                          :: ddot
  integer                                                             :: i, k

  trace = 2 * ddot(nobs * (nobs + 1) / 2, A, 1, B, 1)
  k = 0
  do i = 1, nobs
     k = k + i
     trace = trace - A(k) * B(k)
  end do
end function traceAxB

function traceAxBdiag(A, B, nobs) result(trace)
  implicit none
  ! This function calculates the trace of AxB when A is symmetric
  ! and packed of size nobsxnobs and B is diagonal.
  double precision, dimension(:), intent(in)                          :: A, B
  integer, intent(in)                                                 :: nobs
  double precision                                                    :: trace
  double precision, external                                          :: ddot
  integer                                                             :: i, k
  trace = 0.d0
  k = 0
  do i = 1, nobs
     k = k + i
     trace = trace + A(k) * B(i)
  end do
end function traceAxBdiag


subroutine dunpack(uplo, n, ap, a, lda, info)
  implicit none
  ! unpacks a double precision vector to a full matrix. The difference with the lapack method is
  ! that here I pack both upper and lower part of the matrix
  character(len = 1), intent(in)                                      :: uplo ! not used just kept for consistency
  integer, intent(in)                                                 :: n, lda
  double precision, dimension(:), intent(in)                          :: ap
  double precision, dimension(:,:), intent(inout)                     :: a
  integer, intent(out)                                                :: info
  integer                                                             :: i, j, k

  if (n .ne. lda) then
     info = 2
     return
  end if
  info = 1
  if ((size(a, dim = 1) .ne. n) .or. (size(a, dim = 2) .ne. n)) then
     return
  end if
  k = 1
  do i = 1, n
     do j = 1, i
        a(i,j) = ap(k)
        a(j,i) = ap(k)
        k = k + 1
     end do
  end do
  info = 0

end subroutine dunpack

