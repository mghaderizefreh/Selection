! sample gammete
! assumptions:
!
! equidistance loci (no position is needed to be passed)
! given nloci  => nsegments = nloci-1
!
! no chromosome length needed (as probability of number of chiasma is passed)
!
! no mutation (this is done independently in a different subroutine)
!
!============================================================================
subroutine sampleGamete(parent, id, igam, parentGen, nloci, nblock, &
     nchiasma, chiasmaCum0, chiasmaCumP, istore, &
     offGenInp, positions, samepos)
  use constants
  implicit none

  integer, intent(in)   :: istore  ! how SNP are stored (1=binary/0=integer)

  integer, intent(in) :: parent, id, igam
  integer, intent(in) :: nloci, nblock !nblock = loci if storage is as integer
  integer, dimension(:,:,:), intent(inout), target :: parentGen

  ! maximum number of chiasma (after that the prob is so low that assumed zero)
  integer, intent(in) :: nchiasma

  real(KINDR), intent(in) :: chiasmacum0  ! probability of no recombination
  ! cumulative prob of i recombination
  real(KINDR), dimension(:), intent(in) :: chiasmacumP  

  ! if the genotyps of offspring are store in a different array
  integer, dimension(:,:,:), intent(inout), optional, target :: offGenInp

  ! if the loci are not equidistance, then the position (Morgan or cM).
  real(KINDR), dimension(:), intent(in), optional :: positions
  integer, intent(in), optional :: samepos   ! flag showing if there 

  real(KINDR)  :: rand
  integer :: ichiasma
  !  real(KINDR)   :: val1
  integer, dimension(nchiasma) :: recombination_points, recombination_block,&
       recombination_bit
  integer :: ipos, i, j, k, g, ihap, frompos, frombit, topos, array_pos, &
       bit_pos
  integer, dimension(:,:,:), pointer :: OffGen

  real(KINDR) :: startpos, endpos, interval, rpos
  logical :: accepted

  external :: getPosition

  if (nloci <= 1) then
     write(STDERR, *) " error. subroutine need at least two loci"
     stop 2
  end if

  ! Check if the array containing the parents and offspring genotypes are 
  ! the same or different
  if (present(OffGenInp)) then
     OffGen => OffGenInp   ! offspring genotype are stored in different array
  else
     OffGen => ParentGen   ! offs. gen. are stored in same array as parents' gen.
     if (id == parent) then
        write(STDERR, *) "ERROR:"
        write(STDERR, *) &
             " genotypes for parent and offspring are stored in the same array"
        write(STDERR, *) &
             " Then parent and offspring cannot have the same identification"
        stop 2
     end if
  end if

  recombination_points(1) = 0
  recombination_bit(1) = -1
  recombination_block(1) = 0

  ichiasma = 0
  ! sampling the number of recombination events
  ! i.e., assigning {0,...,nchiasma} randomly to ichiasma
  call random_number(rand)
  if (rand <= chiasmacum0) then
     ichiasma = 0
  else
     ichiasma = nchiasma
     do i = 1, (nchiasma-1)
        if (rand <= chiasmacumP(i)) then
           ichiasma = i
           exit
        end if
     end do
  end if

  !  ich = ichiasma  !storing in the global variable (to check progreees)

  ! if NO recombination event then sample one unrecombined haplotype
  ! but it should be randomly taken from one haploid (ihap\in{1,2})
  if (ichiasma == 0) then
     ipos = 1
     call random_number(rand)
     ihap=1                  ! start from paternal haplotype
     if (rand > HALF) ihap = 2  ! start from maternal haplotype
     if (istore == 0) topos = nloci
     if (istore == 1) topos = nblock

     do ipos = 1, topos
        ! the individual inherited the entire haplotype without recombination
        offGen(id, igam, ipos) = parentGen(parent, ihap, ipos)
     end do
     return
  end if

  recombination_points(1:nchiasma)=0

  !=================================
  ! if the code reaches here, then there is at least one recombination event,
  ! sample the recombination points and later the recombined haplotype
  !=================================

  !-------------------------------------------
  ! sample the recombination_points
  !   sampling algorithm depends on whether the loci positions are given or not.
  !-------------------------------------------
  chiasGT0: if (ichiasma > 0) then

     posPresent: if (present(positions)) then
        ! position or loci are given so assumed that loci are not equidistance.
        ! Needs to sample recombination position first, and then find the bracket
        ! containing the position
        ! it is assumed that all loci are ordered in term of their position
        startpos = positions(1)
        endpos = positions(nloci)
        interval = endpos - startpos
        do i = 1, ichiasma
           ! first: sampling the position where the recombination happens (it
           !   must be between ]startpos-endpos]
           ! the recombination cannot happen before the position of the first loci
           accepted = .false.
           do while (.not.accepted)! find a non zero random number
              call random_number(rand)
              !this to avoid getting an exactly zero value (no point to have 
              ! the recombination in the first loci)
              if (rand > ZERO ) accepted = .TRUE.  
           end do
           ! the position where the recombination will be
           rpos = startpos + interval * rand 

           ! second: finding the bracket where the recombination happens

           ! the bracket where the recombination happend (set initially to zero 
           !  before it is assigned)
           ipos = 0 
           accepted = .false.
           j = 1      ! searching from
           k = nloci  ! searching until
           if (rpos == positions(nloci)) then
              ! the recombination happens exactly in the last loci. Then recom 
              !  is in last bracket
              ipos = nloci-1 
              accepted = .true.
           end if

           ! finding the bracket where the recombination happens
           accepting: do while(.not.accepted)  
              ipos = int((j + k) / 2)

              ! locus j and k are adjacent. Recombination happens at bracket 
              !  ipos which is between the loci j and k.
              if (ipos == j) then 
                 accepted = .true.
                 exit
              else
                 if (rpos > positions(ipos)) j = ipos
                 if (rpos < positions(ipos)) k = ipos
                 if (rpos == positions(ipos)) then 
                    ! the recombination is exactly at this locus

                    ipos = ipos - 1  ! the recombination is in the bracket before
                    accepted = .true.
                    exit
                 end if
              end if

           end do accepting

           ! third: if there is at least a pair of loci sharing the same position
           !   (i.e. flag samepos > 0) , check if this position where the 
           !   recombination happen is shared by several loci.
           ! if so, make the recombination just before the first locus.
           !  (in this way it is ensured that all loci in the same position will
           !  alwyas be inherited together)
           if (samepos > 0 .and. ipos > 1 ) then 
              !if the recom is in the first bracket there is no need
              ! to check lower bracket

              k = 0
              do while(k == 0)
                 if (positions(ipos) == positions(ipos - 1)) then 
                    !the prev loci has the same positions. lower the recomb

                    ipos = ipos - 1
                 else
                    k = 1
                 end if
              end do
           end if
           recombination_points(i) = ipos

        end do
     else
        ! no position is given. So assuming that all loci are equidistance and
        !  no loci share the same position
        ! there is no need to sample the recombination position, but only the
        !  bracket where the recombination happens
        do i = 1, ichiasma
           call random_number(rand)
           ipos = int(rand * (nloci-1)) + 1
           if (ipos >= nloci) ipos = nloci - 1 
           ! just in case the random number was exactly 1.0

           recombination_points(i) = ipos
        end do
     end if posPresent

     !-------------------------------------------
     !   now sorting the recombination  (buble sort ok: as the number of 
     !    recombination events will be very few )
     !-------------------------------------------
     if (ichiasma > 1) then
        do i = 1, (ichiasma-1)
           do j = (i+1), ichiasma
              if (recombination_points(i) > recombination_points(j))then
                 k = recombination_points(i)
                 recombination_points(i) = recombination_points(j)
                 recombination_points(j) = k
              end if
           end do
        end do
     end if

     !-------------------------------------------
     !   if storage is in bit then map the recombinations points
     !-------------------------------------------
     if (istore == 1) then
        do i = 1, ichiasma
           topos = recombination_points(i)
           call getposition(topos,NBITS,array_pos,bit_pos)
           recombination_block(i) = array_pos
           recombination_bit(i) = bit_pos
        end do
     end if
  end if chiasGT0

  !===========================================================================
  ! now copying the gamete to the appropriate place in the offspring
  !-------------------------------------------

  frompos = 1
  call random_number(rand)
  ihap = 1                 ! start from paternal haplotype
  if(rand > HALF) ihap = 2  ! start from maternal haplotype

  ! allele stored as integer
  !---------------------------------------------
  if (istore == 0) then
     frompos = 1
     call random_number(rand)
     ihap = 1 ! start from paternal haplotype
     if(rand > HALF) ihap = 2  ! start from maternal haplotype
     i = 0
     do while (frompos <= nloci)

        i = i + 1
        if(i <= ichiasma) topos = recombination_points(i)
        if(i  > ichiasma) topos = nloci

        do ipos = frompos, topos
           ! the individual inherited the haplotype interval [<frompos>: <topos>]
           offGen(id, igam, ipos) = parentGen(parent, ihap, ipos) 
        end do

        ! now resetting indices for next segment
        ! choosing the alternative haplotype as a recombinatio event happen 
        !  at this point (1=> 2 or 2=> 1)
        ihap = 3 - ihap
        frompos = topos + 1
        !     check if next recombinatiopn was in the same bracket
        if (i < ichiasma ) then 
           do while (recombination_points(i) == recombination_points(i+1) )
              ! the recombination in the same bracket (change the other haplotype)
              ihap = 3 - ihap
              i = i + 1
              if (i == ichiasma) exit
           end do
        end if

     end do

     !---------------------------------------------
     ! allele stored as binary
  else

     frompos = 1  !next segment starts in first block (as it is the first loci)
     frombit = NBITS-1 !next segment starts in first bit (as it is the first loci)
     call random_number(rand)
     ihap = 1                   ! start from paternal haplotype
     if (rand > HALF) ihap = 2  ! start from maternal haplotype
     i = 0

     do while (frompos <= nblock)
        i = i + 1              !the next recombination

        if(i <= ichiasma) then
           topos = recombination_points(i)           !where the segment ends
           array_pos = recombination_block(i)
           bit_pos = recombination_bit(i)
        else
           topos = nloci
           call getPosition(topos,NBITS,array_pos,bit_pos)
        end if

        !    copying the block where the recombination happen
        k = parentGen(parent, ihap, frompos)   ! the haplotype from parent

        g = offGen(id,igam,frompos)   ! the haplotype from offspring

        ipos = 0
        j = frombit + 1

        ! insert the bits (starting in ipos:(j-1)) from variable <k> to 
        !  variable <g> starting at ipos
        call mvbits(k,ipos,j, g,ipos)  

        offGen(id,igam,frompos) = g

        ! copy the remaining blocks (if the next recombination is in a 
        !  different block)
        if (array_pos > frompos) then
           j = frompos + 1   !unrecombined block starts here

           do ipos = j, array_pos, 1
              ! the indiv. inherited the haplotype interval [<frompos>: <topos>]
              offGen(id, igam, ipos) = parentGen(parent, ihap, ipos)
           end do
        end if

        !     now resetting indices for next segment
        ihap = 3 - ihap              ! swaping the haplotype inherited from parents

        bit_pos = bit_pos - 1

        ! if the next loci is in a different block (reset the indices)
        if (bit_pos < 0) then
           array_pos = array_pos + 1 ! next posi is in next block
           bit_pos = NBITS - 1     ! next posi starts in first bit of next block
        end if
        frompos = array_pos
        frombit = bit_pos

        ! check if next recombination is in the same bracket
        if (i < ichiasma) then 
           ! the recombination in the same bracket (change to the other haplotype
           do while (recombination_points(i) == recombination_points(i+1)) 
              ihap = 3 - ihap
              i = i + 1
              if (i == ichiasma) exit
           end do
        end if

        if (topos == nloci) exit

     end do

  end if

  !============================================
  return

end subroutine sampleGamete

