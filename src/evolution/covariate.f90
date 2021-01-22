subroutine covariate(nComp, nAnim, indiv, TBV, E, phen, locations, means, cte)
  use constants
  implicit none

  integer, intent(in) :: nAnim, nComp
  integer, dimension(nanim), intent(in) :: indiv
  double precision, dimension(nanim, ncomp), intent(in) :: TBV
  double precision, dimension(ncomp) :: means
  double precision, dimension(:,:), intent(in) :: locations
  integer, dimension(:), allocatable, intent(out) :: ids
  double precision, dimension(:,:), allocatable, intent(out) :: phen

  integer :: i, j

  j = size(locations)
  if (size(locations) == 1) then
     ! first case: all individual are phenotyped at one location
     phen(1:j) = TBV(1:j,1) + E(1:j,1) + locations(1:1) * (TBV(1:j,2) + E(1:j,2))
     cte = sum(phen) / j - (means(1) + means(2) * locations(1,1))
  else!if ((j < nAnim) .and. (size(locations, 1) .eq. 1)) then
     ! case 2: all individuals have phenotype at a few locations
     write(STDERR, *) "error"
     write(STDERR, *) "not implemented for size(locations) > 1"
     write(STDERR, *) "This needs to be discussed"
     stop 2
  end if
  
end subroutine covariate

