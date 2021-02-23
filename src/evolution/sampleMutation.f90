subroutine sampleMutation(id, igam, genotypes, nloci, nmut, mutCum0,&
     mutCumP, istore)
  use constants
  implicit none

  integer, intent(in) :: istore  ! how SNP are stored

  integer, intent(in) :: id,igam
  integer, intent(in) :: nloci !nblock = nloci if storage is as integer
  integer, dimension(:,:,:), intent(inout) :: genotypes

  ! max # of mutations (after that the prob is so low that assumed zero)
  integer, intent(in) :: nmut
  real(KINDR), intent(in) :: mutcum0  ! probability of no mutation

  ! cumulative probability of number of mutations
  real(KINDR), dimension(:), intent(in) :: mutcumP

  REAL(KINDR)   :: rand
  integer       :: imut
!  real(KINDR)   :: val1
  integer, dimension(nmut) :: mutation_points
  integer :: ipos, g, i, array_pos, bit_pos !,j, k, ihap, frompos, topos

  external :: getPosition

  ! sampling the number of recombination events
  call random_number(rand)
  if (rand <= mutcum0) then
     imut = 0
  else
     imut = nmut
     do i = 1, (nmut - 1)
        if (rand <= mutcumP(i)) then
           imut= i
           exit
        end if
     end do
  end if

!  imu = imut

  ! find the recombination_points (do we mean mutation points?)
  mutation_points(1) = 0
  if (imut > 0) then

     !sampling the positions
     i = 1
     do while(i <= imut)
        call random_number(rand)

        ipos = int(rand * nloci) + 1
        if (ipos >= nloci) ipos = nloci - 1 ! just in case the rand was exactly 1.0

        mutation_points(i) = ipos
        i = i + 1

     end do

  end if
  !================================
  ! now mutate  the gamete to the appropriate place in the offspring

  ! allele stored as integer
  !---------------------------------------------
  if (istore == 0) then
     i = 0
     do i = 1, imut
        ipos = mutation_points(i)
        genotypes(id, igam, ipos) = 3 - genotypes(id, igam, ipos)
     end do

     !---------------------------------------------
     ! allele stored as binary
  else
     do i = 1, imut
        ipos = mutation_points(i)
        call getPosition(ipos, NBITS, array_pos, bit_pos)
        g = genotypes(id, igam, array_pos)
        if (Btest(g,bit_pos)) then
           g = ibclr(g,bit_pos) !if 1 => 0
        else
           g = ibset(g,bit_pos) !if 0 => 1
        end if
        genotypes(id, igam, array_pos) = g
     end do
  end if

  !============================================

end subroutine sampleMutation

