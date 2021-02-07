subroutine selectbyIntercept(nanim, indiv, sex, n_m, n_fpm, male, female, ranEff)
  use constants
  use quickSort
  implicit none

  integer, intent(in) :: nanim
  integer, dimension(nanim), intent(in) :: indiv
  logical, dimension(:), intent(inout) :: sex
  type(doublePre_Array), dimension(:), intent(in) :: raneff
  integer, intent(in) :: n_m, n_fpm
  integer, dimension(:), intent(out) :: male, female

  integer :: nmale, nfemale, index
  integer :: i
  real(KINDR), dimension(:), allocatable :: tempR
  integer, dimension(:), allocatable :: tempI

  if (size(raneff) > 3) then
     write(STDERR, *) " ERROR:"
     write(STDERR, *) " selectByIntercept not implemented for raneff > 3"
     stop 2
  end if
  index = merge(1, 2, size(raneff) == 1)
  nmale = count(sex)
  nfemale = size(sex) - nmale
  i = max(nmale, nfemale) !using one array for both female and male per type
  allocate(tempI(i))
  allocate(tempR(i))
  
  tempI(1:nmale) = pack(indiv, sex) ! true is male
  tempR(1:nmale) = raneff(index)%level(tempI(1:nmale))
  call sortrx(nmale, tempR, tempI)
  i = nmale - n_m + 1
  male(1:n_m) = tempI(i:nmale) ! sort is ascending

  tempI(1:nfemale) = pack(indiv, .not.sex)
  tempR(1:nfemale) = raneff(index)%level(tempI(1:nfemale))
  call sortrx(nfemale, tempR, tempI)
  i = nfemale - n_fpm * n_m + 1
  female(1:(n_fpm * n_m)) = tempI(i:nfemale) + nmale !offset indices by nmale
  
end subroutine selectbyIntercept



subroutine SelectBySlope(nanim, indiv, sex, n_m, n_fpm, male, female, ranEff)
  use constants
  use quickSort
  implicit none
  
  integer, intent(in) :: nanim
  integer, dimension(nanim), intent(in) :: indiv
  logical, dimension(:), intent(inout) :: sex
  type(doublePre_Array), dimension(:), intent(in) :: raneff
  integer, intent(in) :: n_m, n_fpm
  integer, dimension(:), intent(out) :: male, female

  integer :: nmale, nfemale, index
  integer :: i
  real(KINDR), dimension(:), allocatable :: tempR
  integer, dimension(:), allocatable :: tempI

  if (size(raneff) .ne. 3) then
     write(STDERR, *) " ERROR:"
     write(STDERR, *) " selectBySlope must have size(raneff) = 3"
     stop 2
  end if
  index = 1
  nmale = count(sex)
  nfemale = size(sex) - nmale
  i = max(nmale, nfemale)
  allocate(tempI(i))
  allocate(tempR(i))
  
  tempI(1:nmale) = pack(indiv, sex)
  tempR(1:nmale) = raneff(index)%level(tempI(1:nmale))
  call sortrx(nmale, tempR, tempI)
  i = nmale - n_m + 1
  male(1:n_m) = tempI(i:nmale)

  tempI(1:nfemale) = pack(indiv, .not.sex)
  tempR(1:nfemale) = raneff(index)%level(tempI(1:nfemale))
  call sortrx(nfemale, tempR, tempI)
  i = nfemale - n_fpm * n_m + 1
  female(1:(n_fpm * n_m)) = tempI(i:nfemale) + nmale

end subroutine SelectBySlope





function makePedigree(n_m, n_fpm, n_opf, male, female, start) result(pedigree)
  use constants
  implicit none
  integer, intent(in) :: n_m, n_fpm, n_opf, start
  integer, dimension(:), intent(in) :: male, female
  integer, dimension((n_m*n_fpm*n_opf),3) :: pedigree
  integer :: i, no, N, nf

  no = n_fpm * n_opf
  nf = n_m * n_fpm
  N = n_m * no
  do i = 1, n_m
     pedigree(((i-1)*no+1):(i*no), 2) = male(i)
  end do
  do i = 1, nf
     pedigree(((i-1)*n_opf+1):(i*n_opf), 3) = female(i)
  end do
  do i = 1, N
     pedigree(i, 1) = start + i
  end do
end function makePedigree
