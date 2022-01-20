subroutine getGen0Variance(nvar, nran, nanim, nobs, interval, chalvl, vars, y,&
     theta)
  use constants, only : KINDR, variances, STDERR, ONE
  use math, only : variance
  implicit none
  !! ================ variable definitions  ================ !!
  integer, intent(in) :: nvar, nran, nanim, nobs
  real(KINDR), dimension(1:2), intent(in) :: interval, chalvl
  type(variances), intent(in) :: vars
  real(KINDR), dimension(1:nobs) :: y
  real(KINDR), dimension(1:(nvar+1)) :: theta

  real(KINDR) :: val1, val2

  if (nanim > 0) continue
  if (nran == 1) then
     if (nvar .ne. 1) then
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) "with nran = 1, one variance component is expected&
             &(nvar = 1)"
        stop 2
     end if
  elseif (nran == 2) then
     write(STDERR, '(a)') "Error:"
     write(STDERR, *) "This is not implented"
     stop 2
  end if

  if (nran == 1) then
     val1 = (interval(1) + interval(2)) / 2 ! x=x_middle
     ! nominator
     val2 = val1 * val1 * vars%A(1) + vars%A(2) + 2 * val1 * vars%cov(1,2)
     ! denominator
     val1 = val2 + val1 * val1 * vars%E(1) + vars%E(2)
     val1 = val2 / val1
     theta(1) = variance(y, nobs) * val1
     theta(2) = variance(y, nobs) * (ONE - val1)
  elseif (nran == 3) then
     val1 = (interval(2) - interval(1)) / (chalvl(2) - chalvl(1))
     theta(1) = vars%A(1) * (val1 ** 2)
     theta(2) = vars%A(2)
     theta(3) = vars%E(1) * (val1 ** 2)
     theta(nvar+1) = vars%E(2)
     if (nvar == 3) then
     elseif (nvar == 4) then
        theta(4) = vars%cov(1,2)
     else
        write(STDERR, '(a)') "Error:"
        write(STDERR, *) "This is not implemented(2)"
     end if
  end if
end subroutine getGen0Variance

subroutine getTrueVariance(nvar, nran, nanim, ncomp, nobs, interval, vars,&
     tbv, y, theta)
  use constants, only : KINDR, variances, ONE
  use math, only : variance, covariance
  implicit none
  !! ================ variable definitions  ================ !!
  integer, intent(in) :: nvar, nran, nanim, ncomp, nobs
  real(KINDR), dimension(1:2), intent(in) :: interval
  type(variances), intent(in) :: vars
  real(KINDR), dimension(1:ncomp,1:nanim), intent(in) :: tbv
  real(KINDR), dimension(1:nanim), intent(in) :: y
  real(KINDR), dimension(1:(nvar+1)), intent(out) :: theta

  real(KINDR) :: val1, val2
 

  if (nran .eq. 1) then
     val1 = (interval(1) + interval(2)) / 2 ! x=x_middle
     ! nominator
     val2 = val1 * val1 * variance(tbv(1:nanim,1), nanim) + &
          variance(tbv(1:nanim,2), nanim) + 2 * val1 * &
          covariance(tbv(1:nanim,1), tbv(1:nanim,2), nanim)
     ! denominator
     val1 = val2 + val1 * val1 * vars%E(1) + vars%E(2)
     val1 = val2 / val1 ! heritability
     theta(1) = variance(y, nobs) * val1
     theta(2) = variance(y, nobs) * (ONE - val1)
  elseif (nran .eq. 3) then
     val1 = sum(tbv(1:nanim,1))/nanim
     val2 = sum(tbv(1:nanim,2))/nanim
     theta(1) = variance(tbv(1:nanim,1), nanim)
     theta(2) = variance(tbv(1:nanim,2), nanim)
     theta(3) = vars%E(1)
     theta((nvar+1)) = vars%E(2)
     if (nvar .eq. 3) then
     elseif (nvar .eq. 4) then
        theta(4) = covariance(tbv(1:nanim,1), tbv(1:nanim, 2), nanim)
     end if
  end if
  
end subroutine getTrueVariance
