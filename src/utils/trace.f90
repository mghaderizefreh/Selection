function traceA(A, nobs) result(trace)
  ! simply calcuates the trace of A with size nobsxnobs
  use constants
  implicit none
  real(KINDR), dimension(:), intent(in)                          :: A
  integer, intent(in)                                                 :: nobs
  real(KINDR)                                                    :: trace
  integer                                                             :: i, j
  trace = 0.d0
  j = 0
  do i = 1, nobs
     j = j + i
     trace = trace + A(j)
  end do
end function traceA

function traceAxB(A, B, nobs) result(trace)
  use constants
  implicit none
  ! This function calculates the trace of AxB when A and both
  ! are symmetric and of size nobsxnobs. It uses ddot of BLAS
  real(KINDR), dimension(:), intent(in)                          :: A, B
  integer, intent(in)                                                 :: nobs
  real(KINDR)                                                    :: trace
  real(KINDR), external                                          :: ddot
  integer                                                             :: i, k

  trace = 2 * ddot(nobs * (nobs + 1) / 2, A, 1, B, 1)
  k = 0
  do i = 1, nobs
     k = k + i
     trace = trace - A(k) * B(k)
  end do
end function traceAxB

function traceAxBdiag(A, B, nobs) result(trace)
  use constants
  implicit none
  ! This function calculates the trace of AxB when A is symmetric
  ! and packed of size nobsxnobs and B is diagonal.
  real(KINDR), dimension(:), intent(in)                          :: A, B
  integer, intent(in)                                                 :: nobs
  real(KINDR)                                                    :: trace
  real(KINDR), external                                          :: ddot
  integer                                                             :: i, k
  trace = 0.d0
  k = 0
  do i = 1, nobs
     k = k + i
     trace = trace + A(k) * B(i)
  end do
end function traceAxBdiag


