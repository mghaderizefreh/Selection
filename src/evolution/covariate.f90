subroutine covariate(nComp, nObs, nAnim, nLox, TBV, E, phen, locations, means)
  use constants
  implicit none

  integer, intent(in) :: nComp, nObs, nAnim, nLox
  real(KINDR), dimension(nanim, ncomp), intent(in) :: TBV
  real(KINDR), dimension(ncomp) :: means
  real(KINDR), dimension(nAnim, nLox), intent(in) :: locations
  real(KINDR), dimension(nobs, ncomp), intent(in) :: E
  real(KINDR), dimension(nobs), intent(out) :: phen

  real(KINDR), dimension(:), allocatable, save :: cte
  integer :: i, j, k
  
  if (nComp .ne. 2) then
     write(STDERR, '(a)') "Error"
     write(STDERR, *) " 'covariate' only works for ncomp = 2 at the moment"
     stop 3
  end if

  if (.not.allocated(cte)) then
     call alloc1D(cte, nComp, "cte", "covariate")
     cte(1:nComp) = ZERO
     do i = 1, nComp
        cte(i) = means(i) - sum(TBV(1:nAnim, i))/nAnim&
             - sum(E(1:nObs, i))/nObs
     end do
  end if

  ! all individuals have phenotype at a few locations
  k = 0
  do i = 1, nAnim
     do j = 1, nlox
        k = 1 + k
        phen(k) = TBV(i, 2) + E(k, 2) + cte(2) + locations(i, j) * &
             (TBV(i, 1) + E(k, 1) + cte(1))
     end do
  end do

end subroutine covariate

