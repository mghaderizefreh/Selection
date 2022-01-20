subroutine selectParents(nanim, indiv, sex, n_m, n_fpm, male, female,&
     effects)
  use constants, only : KINDR, alloc1I, alloc1D
  use quickSort, only : sortrx
  implicit none

  integer, intent(in) :: nanim
  integer, dimension(nanim), intent(in) :: indiv
  logical, dimension(1:nAnim), intent(inout) :: sex
  real(KINDR), dimension(1:nAnim), intent(in) :: effects ! breeding values
  integer, intent(in) :: n_m, n_fpm
  integer, dimension(1:(n_m*n_fpm)), intent(out) :: female
  integer, dimension(1:n_m), intent(out) :: male

  integer :: nmale, nfemale
  integer :: i, j
  real(KINDR), dimension(:), allocatable, save :: mR, fR ! real array for (fe)male
  integer, dimension(:), allocatable, save :: mI, fI, tempMI, tempFI

  nmale = count(sex)
  nfemale = nAnim - nmale
  if (.not.allocated(mI)) then
     call alloc1I(mI, nmale, "mI", "select")
     call alloc1I(fI, nfemale, "fI", "select")
     call alloc1D(mR, nmale, "mR", "select")
     call alloc1D(fR, nfemale, "fR", "select")
     call alloc1I(tempMI, nmale, "tempMI", "select")
     call alloc1I(tempFI, nfemale, "tempFI", "select")
  end if
  
  mI(1:nmale) = pack(indiv, sex) ! true is male, mI holds index of males
  mR(1:nmale) = effects(mI(1:nmale))! mR holds the bv of males
  call sortrx(nmale, mR, tempMI) ! bv are sorted (ascending)
  i = nmale - n_m + 1
  ! male holds the index of best males
  male(1:n_m) = mI(tempMI(i:nmale)) ! sort is ascending
!  if(verbose) write(STDOUT, *) "best male (ind,val)", male(n_m), &
!       mR(tempMI(nmale)), sex(male(n_m))
!  if(verbose) write(STDOUT, *) "2nd best male      ", male(n_m-1),&
!       mR(tempMI(nmale-1)), sex(male(n_m-1))

  fI(1:nfemale) = pack(indiv, (.not.sex))
  fR(1:nfemale) = effects(fI(1:nfemale))
  call sortrx(nfemale, fR, tempFI)
  j = n_fpm * n_m
  i = nfemale - j + 1
  female(1:j) = fI(tempFI(i:nfemale))
!  if(verbose) write(STDOUT, *) "best female (ind,val)", female(j), &
!       fR(tempFI(nfemale)), sex(female(j))
!  if(verbose) write(STDOUT, *) "2nd best female      ", female(j-1), &
!       fR(tempFI(nfemale-1)), sex(female(j-1))

end subroutine selectParents

function makePedigree(n_m, n_fpm, n_opf, male, female, start) result(pedigree)
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
