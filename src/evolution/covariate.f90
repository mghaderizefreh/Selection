subroutine covariate(nComp, nAnim, TBV, E, phen, locations, means, cte)
  use constants
  implicit none

  integer, intent(in) :: nAnim, nComp
!  integer, dimension(nanim), intent(in) :: indiv
  real(KINDR), dimension(nanim, ncomp), intent(in) :: TBV
  real(KINDR), dimension(ncomp) :: means
  real(KINDR), dimension(:,:), intent(in) :: locations, E
!  integer, dimension(:), allocatable, intent(out) :: ids
  real(KINDR), dimension(:), intent(out) :: phen
  real(KINDR), intent(out) :: cte
  
  integer :: i, j, k

  j = size(locations)
  
  if (size(locations) == 1) then
     i = nAnim
     ! first case: all individual are phenotyped at one location
     phen(1:i) = TBV(1:i,1) + E(1:i,1) + locations(1,1) * &
          (TBV(1:i,2) + E(1:i,2))
     cte = sum(phen) / nAnim - (means(1) + means(2) * locations(1,1))
  elseif ((j < nAnim) .and. (size(locations, 1) .eq. 1)) then
     ! case 2: all individuals have phenotype at a few locations
     k = 0
     do i = 1, nAnim
        do j = 1, size(locations, 2)
           k = 1 + k
           phen(k) = TBV(i, 1) + E(k, 1) + locations(1, j) * &
                (TBV(i, 2) + E(k, 2))
        end do
     end do
     cte = sum(phen) / k - means(1)
     cte = cte - means(2) * sum(locations(1,:)) / size(locations)
  elseif (size(locations, 1) .eq. nanim) then
     ! case 2: all individuals have phenotype at a few locations
     k = 0
     do i = 1, nAnim
        do j = 1, size(locations, 2)
           k = 1 + k
           phen(k) = TBV(i, 1) + E(k, 1) + locations(i, j) * &
                (TBV(i, 2) + E(k, 2))
        end do
     end do
     cte = sum(phen) / k - means(1)
     k = size(locations)
     do i = 1, size(locations, 2)
        cte = cte - means(2) * sum(locations(:,i)) / k
     end do
  end if


  
end subroutine covariate

