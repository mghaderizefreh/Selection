subroutine SimulatePhenotype(nAnim, nComp, indiv, TBV, variances, means, &
     locations, ids, phen, proc, cte)

  use constants
  use rng_module
  implicit none

  integer, intent(in) :: nAnim, nComp
  integer, dimension(nanim), intent(in) :: indiv
  real(KINDR), dimension(nanim, ncomp), intent(in) :: TBV
  type(variance), intent(in) :: variances
  real(KINDR), dimension(ncomp) :: means
  real(KINDR), dimension(:,:), intent(in) :: locations
  integer, dimension(:), allocatable, intent(out) :: ids
  real(KINDR), dimension(:,:), allocatable, intent(out) :: phen
  character(len = *) :: proc
  real(KINDR), intent(out) :: cte

  integer :: i, j, k
  real(KINDR), dimension(:),allocatable :: temp
  real(KINDR), dimension(:,:), allocatable :: temp2
  real(KINDR), dimension(:), allocatable :: tempr
  real(KINDR), allocatable, dimension(:,:) :: A, E!, PE

  allocate(temp2(nComp, nComp), tempr(nComp))
  temp2(1:nComp, 1:nComp) = ZERO

  do i = 1, nComp
     temp2(i,i) = variances%E(i)
  end do
  tempr(1:nComp) = ZERO

  j = size(locations)
  write(6, *) j
  if (j == 1) then
     i = nAnim
     allocate(phen(i,1), ids(i), E(i,nComp))
     ids = indiv
  elseif ((j < nAnim) .and. (size(locations, 1) .eq. 1)) then
     i = j * nAnim
     allocate(phen(i,2), ids(i))
     do k = 1, j
        phen(k:i:j, 1) = locations(1, k)
        ids(k:i:j) = indiv(1:nAnim)
     end do
  elseif (size(locations, 1) == nAnim) then
     i = size(locations, 2)
     allocate(phen(j,2), ids(j))
     write(6, *) "i", i
     write(6, *) "phen", shape(phen)
     do k = 1, i
        write(6, *) 'k, phen to hold', k, shape(phen(k:j:i,1))
        phen(k:j:i, 1) = locations(1:nAnim, k)
        ids( ((k - 1) * nAnim + 1) : (k * nAnim) ) = indiv(1:nAnim)
     end do
  else
     write(STDERR, *) "Error: The format of 'locations' is wrong"
     write(STDERR, *) "Use either:"
     write(STDERR, *) " - [size:    1 x 1   ] all animals in the same location"
     write(STDERR, *) " - [size:    1 x nLox] all have phenotypes in nlox places"
     write(STDERR, *) " - [size:nAnim x nLox] phenotypes at different places"
     stop 2
  end if

  call gnormal(tempr, temp2, nComp, i, E)
   ! todo: implementation for PE if really it is required

  select case (proc(1:3))
  case ("COV", "COv", "CoV", "Cov", "cOV", "cOv", "coV", "cov")
     call covariate(nComp, nAnim, indiv, TBV, E, phen, locations, means, cte)
  case default
     write(STDERR, *) "error:"
     write(STDERR, *) "case '", proc, "' not implemented"
     stop 2
  end select

end subroutine SimulatePhenotype
