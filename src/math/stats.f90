function variance(X, n) result(out)
  use constants
  implicit none
  integer, intent(in) :: n
  real(KINDR), dimension(:), intent(in) :: X
  real(KINDR) :: out
  real(KINDR) :: v
  v = sum(X(1:n)) / n
  out = sum((X(1:n) - v)**2) / (n - 1)
end function variance

function covariance(X, Y, n) result(out)
  use constants
  implicit none
  integer, intent(in) :: n
  real(KINDR), dimension(n), intent(in) :: X, Y
  real(KINDR) :: out
  out = sum(X(1:n) * Y(1:n)) / n - &
       sum(X(1:n)) * sum(Y(1:n)) / (n**2)
end function covariance

function correlation(X, Y, n) result(out)
  use constants
  implicit none
  integer, intent(in) :: n
  real(KINDR), dimension(n), intent(in) :: X, Y
  real(KINDR) :: out
  real(KINDR) :: v1, v2
  v1 = sum(X(1:n)) / n
  v2 = sum(Y(1:n)) / n
  out = sum((X(1:n) - v1) * (Y(1:n) - v2)) /&
       sqrt(sum((X(1:n) - v1)**2) * sum((Y(1:n) - v2)**2))
end function correlation




