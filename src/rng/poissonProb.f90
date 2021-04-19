!     Last change:  RPW  12 Oct 2010   10:36 am
!
!=================================================================================
! subroutine poissonProb(mu,k,prob)
!
! calculate the probability of k successes from a poisson distribution with mean mu
!
! prob(k|mu) = mu^k * exp(-mu) / (k!)
!
!              k!= 1*2* ... *k
!                  by definition for k=0 ==> k! = 1)
!=================================================================================

subroutine poissonProb(mu,k,prob)

  use constants
  implicit none
  real(KINDR), intent(in):: mu
  integer, intent(in):: k
  real(KINDR), intent(out):: prob

  real(KINDR) :: lkp,lp
  integer :: i

  !-------------------------
  !if mu is exact zero
  !-------------------------
  ! if mu=0 the the only number of events that it could happend is no event
  !  (i.e. prob =1 if no event and zero oyther cases)
  if (mu <= ZERO) then
     if (k==0) then
        prob = ONE
     else
        prob = ZERO
     end if
     return
  end if

  !-------------------------
  ! if mu is gretaer than zero
  !-------------------------
  prob = ZERO
  if (k == 0) prob = exp(-mu)

  if (k > 0 ) then
     lkp = ZERO
     if (k > 1) then
        do i = 2, k
           lkp = lkp + log(i + ZERO) ! changed from log(dble(i))
        end do
     end if
     lp = -mu + (k * log(mu)) - lkp
     prob = exp(lp)
  end if

end subroutine poissonProb
