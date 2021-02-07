subroutine covariate(nComp, nobs, nAnim, TBV, E, phen, locations, row, col, means)
  use constants
  implicit none

  integer, intent(in) :: nAnim, nComp, nobs
  integer, intent(in) :: row, col
  !  integer, dimension(nanim), intent(in) :: indiv
  real(KINDR), dimension(nanim, ncomp), intent(in) :: TBV
  real(KINDR), dimension(ncomp) :: means
  real(KINDR), dimension(row,col), intent(in) :: locations
  real(KINDR), dimension(nobs, ncomp), intent(in) :: E
  !  integer, dimension(:), allocatable, intent(out) :: ids
  real(KINDR), dimension(nobs), intent(out) :: phen

  real(KINDR), dimension(:), allocatable :: cte
  logical , save :: gen1 = .true.
  integer :: i, j, k
  allocate(cte(ncomp))
  cte(1:nComp) = ZERO
  if (gen1) then
     do i = 1, nComp
        cte(i) = means(i) - (sum(TBV(1:nAnim,i)) + sum(E(1:nAnim, i)))/nAnim
     end do
     gen1 = .false.
  end if

  j = row * col
  if (row == 1) then
     i = nAnim
     ! first case: all individual are phenotyped at one location
     phen(1:i) = TBV(1:i,1) + E(1:i,1) + cte(1) + locations(1,1) * &
          (TBV(1:i,2) + E(1:i,2) + cte(2))
  elseif ((j < nAnim) .and. (row .eq. 1)) then
     ! case 2: all individuals have phenotype at a few locations
     k = 0
     do i = 1, nAnim
        do j = 1, col
           k = 1 + k
           phen(k) = TBV(i, 1) + E(k, 1) + cte(1) + locations(1, j) * &
                (TBV(i, 2) + E(k, 2) + cte(2))
        end do
     end do
  elseif (row .eq. nanim) then
     ! case 2: all individuals have phenotype at a few locations
     k = 0
     do i = 1, nAnim
        do j = 1, col
           k = 1 + k
           phen(k) = TBV(i, 1) + E(k, 1) + cte(1) + locations(i, j) * &
                (TBV(i, 2) + E(k, 2) + cte(1))
        end do
     end do
  end if

end subroutine covariate

