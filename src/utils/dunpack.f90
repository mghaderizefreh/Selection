subroutine dunpack(uplo, n, ap, a, lda, info)
  use constants
  implicit none
  ! unpacks a real(KINDR) vector to a full matrix. The difference with the lapack method is
  ! that here I pack both upper and lower part of the matrix
  character(len = 1), intent(in) :: uplo ! not used just kept for consistency
  integer, intent(in) :: n, lda
  real(KINDR), dimension(1:(n*(n+1)/2)), intent(in) :: ap
  real(KINDR), dimension(1:n,1:n), intent(out) :: a
  integer, intent(out) :: info
  integer :: i, j, k

  if (uplo == 'L') then
  end if

  if (n .ne. lda) then
     info = 2
     write(STDOUT, '(a,i3)') "  warning: error in dunpack. info:", info
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

