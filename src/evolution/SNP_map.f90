!     Last change:  RPW  15 Jan 2010    2:59 pm
! subroutine getPosition(position,nbits,array_pos,bit_pos)
! boolean    function getBit(the_integer, bit_pos)
! subroutine putAllele(allele_array,id,bit_pos,allele)
! subroutine putGenotype(alelle_array1,allele_array2,id,bit_pos,allele1,allele2)
! subroutine putAlleleLocus(allele_array,locus,allele)
! subroutine putGenotypeLocus(alelle_array1,allele_array2,SNP,allele1,allele2)
! subroutine genotypeBin2Int(genotype,allele1,allele2)
! subroutine genotypeInt2Bin(genotype,allele1,allele2)
! subroutine getAllGenotypesInt(genotypes,alleles1,alleles2,nbit,locus)
!==============================================================================

subroutine getPosition(pos,nbits,array_pos,bit_pos)
  implicit none
  integer, intent(in) :: pos         ! the n bit to be found
  integer, intent(in) :: nbits       ! number of bit per integer
  integer, intent(out):: array_pos   ! the integer in the array
  integer, intent(out):: bit_pos     ! the bit within the integer 
  integer :: remaining

  array_pos = int(pos/nbits)
  remaining = mod(pos,nbits)
  if(remaining == 0) then
     bit_pos=0
  else
     array_pos = array_pos+1
     bit_pos   = nbits-remaining
  endif
end subroutine getPosition

!=======================================================================
!function getBit(the_integer, bit_pos) result(gb)
!  implicit none
!  integer, intent(in) :: the_integer  ! the integer in the array
!  integer, intent(in) :: bit_pos      ! the bit within the integer
!  logical :: gb
!  gb = btest(the_integer,bit_pos)
!end function getBit

!=======================================================================
!subroutine putAllele(allele_array,id,bit_pos,allele)
!  implicit none
!  integer, dimension(:),intent(inout) :: allele_array
!  integer              ,intent(in)    :: id
!  integer              ,intent(in)    :: bit_pos
!  logical              ,intent(in)    :: allele
!
!  if(.not. allele) then  !allele is 0
!     allele_array(id)=ibclr(allele_array(id), bit_pos)
!  else                   !allele is 1
!     allele_array(id)=ibset(allele_array(id), bit_pos)
!  endif
!end subroutine putAllele
!
!!=======================================================================
!subroutine putGenotype(allele_array1,allele_array2,id,bit_pos,allele1,allele2)
!  implicit none
!  integer, dimension(:),intent(inout) :: allele_array1
!  integer, dimension(:),intent(inout) :: allele_array2
!  integer          ,intent(in)    :: id
!  integer          ,intent(in)    :: bit_pos
!  logical          ,intent(in)    :: allele1,allele2
!
!  call putAllele(allele_array1,id,bit_pos,allele1)
!  call putAllele(allele_array2,id,bit_pos,allele2)
!end subroutine putGenotype
!
!!=======================================================================
!subroutine putAlleleLocus(allele_matrix,id,locus,allele)
!  !--------------------------------------------------------------------
!  ! same as putAllele but the input is the locus number. Hence, this subroutine 
!  ! needs to find the specific bit which will store the allele whilst this
!  ! subroutine seems cleaner it needs to do two other operations hence more
!  ! computationally intensive. This could seriously affect the performance of
!  ! the programme if a lot of loci are considered per animal
!  !--------------------------------------------------------------------
!  implicit none
!  integer, dimension(:,:),intent(inout),target :: allele_matrix
!  integer          ,intent(in)    :: id
!  integer          ,intent(in)    :: locus
!  logical          ,intent(in)    :: allele
!
!  integer             :: array_pos
!  integer             :: bit_pos
!  integer             :: nbits
!
!  nbits = bit_size(allele_matrix(1,1))                     !bits per integer
!  !finding the column and position where locus is stored
!  call getPosition(locus,nbits,array_pos,bit_pos)
!  call putAllele(allele_matrix(:,array_pos),id,bit_pos,allele) !storing the allele
!end subroutine putAlleleLocus
!
!!=======================================================================
!subroutine putGenotypeLocus(allele_matrix1,allele_matrix2,id,locus,allele1,allele2)
!  implicit none
!  integer, dimension(:,:),intent(inout) :: allele_matrix1
!  integer, dimension(:,:),intent(inout) :: allele_matrix2
!  integer          ,intent(in)    :: id
!  integer          ,intent(in)    :: locus
!  logical          ,intent(in)    :: allele1,allele2
!
!  integer             :: array_pos
!  integer             :: bit_pos
!  integer             :: nbits
!
!  nbits=bit_size(allele_matrix1(1,1))                 !bits per integer
!  !finding the column and position where locus is stored
!  call getPosition(locus,nbits,array_pos,bit_pos)     
!  call putGenotype(allele_matrix1(:,array_pos),allele_matrix2(:,array_pos), &
!       id,bit_pos,allele1,allele2)
!
!end subroutine putGenotypeLocus
!
!!=======================================================================  
!subroutine genotypeBin2Int(genotype,allele1,allele2)
!  integer, intent(out) :: genotype
!  logical, intent(in)  :: allele1
!  logical, intent(in)  :: allele2
!
!  !    genotype 1=AA,2=AB, 3=BB
!  genotype=1
!  if(allele1) genotype=genotype+1
!  if(allele2) genotype=genotype+1
!end subroutine genotypeBin2Int
!
!!=======================================================================  
!subroutine genotypeInt2Bin(genotype,allele1,allele2)
!  integer, intent(in)    :: genotype
!  logical, intent(out)   :: allele1
!  logical, intent(out)   :: allele2
!
!  !    genotype 1=AA,2=AB, 3=BB
!  if(genotype == 1) then
!     allele1=.FALSE.
!     allele2=.FALSE.
!  elseif(genotype == 2) then
!     allele1=.TRUE.
!     allele2=.FALSE.
!  elseif(genotype == 3) then
!     allele1=.TRUE.
!     allele2=.TRUE.
!  endif
!end subroutine genotypeInt2Bin
!!======================================================================
!subroutine getAllGenotypesInd(genotypes,genotyped,alleles1,alleles2,nbit,nSNP)
!  !  it returns all genotypes for a given animal
!  integer, dimension(:), intent(out)   :: genotypes
!  !the vector row stating if with individual is genotyped for this locus
!  integer, dimension(:), intent(in )   :: genotyped
!  !the vector row with individual's allele 1
!  integer, dimension(:), intent(in )   :: alleles1
!  !the vector row with individual's allele 2
!  integer, dimension(:), intent(in )   :: alleles2
!  integer              , intent(in)    :: nSNP,nbit
!  integer :: i,blocks,bit_pos,gen,ilocus
!  logical :: all1,all2,accepted
!
!  blocks=int(nSNP/nbit)+1
!  ilocus=0
!  Allloci: do i=1,blocks
!     do bit_pos=(nbit-1),0, -1
!        ilocus=ilocus+1
!        if(ilocus.le.nSNP) then
!           accepted= btest(genotyped(i),bit_pos)
!           if  ( accepted ) then
!              all1 = btest(alleles1(i),bit_pos)
!              all2 = btest(alleles2(i),bit_pos)
!              call genotypeBin2Int(gen,all1,all2)
!              genotypes(ilocus)= gen
!           else
!              genotypes(ilocus)= 0     !genotype missing for this locus
!           endif
!        else
!           exit Allloci
!        endif
!     end do
!  end do Allloci
!end subroutine getAllGenotypesInd
!
!!=====================================================================
!
!
!!======================================================================
!subroutine getAllGenotypesLocus(genotypes,alleles1,alleles2,bit_pos,nanim,genotyped, ngen)
!  !  it returns all genotypes coded as integer for a given SNP
!  !  if the dummy variable <genotyped> is present, then it check if the genotype
!  !   is missing (return 0)
!
!  integer, dimension(:), intent(out)          :: genotypes
!  !the vector col with individual's allele 1
!  integer, dimension(:), intent(in) :: alleles1
!  !the vector col with individual's allele 2
!  integer, dimension(:), intent(in) :: alleles2
!  integer              , intent(in)           :: nanim,bit_pos
!  !the vector col stating if individual is genotyped for this locus
!  integer, dimension ( : ), intent ( in ), optional :: genotyped
!  !vector returning the number of individual for each genotype
!  integer, dimension ( : ), intent ( out ), optional :: ngen
!
!  integer, dimension(3) :: totalgen
!  integer :: id,gen
!  logical :: all1,all2,accepted
!
!  totalgen(:) = 0
!  if (present(genotyped)) then
!     do id=1,nanim
!        accepted=btest(genotyped(id),bit_pos)
!        if  ( accepted ) then
!           all1 = btest(alleles1(id),bit_pos)
!           all2 = btest(alleles2(id),bit_pos)
!           call genotypeBin2Int(gen,all1,all2)
!           genotypes(id) = gen
!           totalgen( gen ) = totalgen( gen )+1
!        else
!           genotypes(id) = 0   !genotype missing for this locus
!        endif
!     end do
!  else
!     do id=1,nanim
!        all1 = btest(alleles1(id),bit_pos)
!        all2 = btest(alleles2(id),bit_pos)
!        call genotypeBin2Int(gen,all1,all2)
!        genotypes(id) = gen
!        totalgen( gen ) = totalgen( gen ) + 1
!     end do
!
!  endif
!  if(present(ngen)) ngen=totalgen
!
!end subroutine getAllGenotypesLocus
!
!!======================================================================
!
!subroutine getAGenotype(genotype,allele1,allele2,bit_pos,genotyped)
!  !  it returns genotypes coded as integer for a given SNP
!  !  if the dummy variable <genotyped> is present, then it check if the genotype
!  !   is missing (return 0)
!
!  integer, intent(out)          :: genotype
!  integer, intent(in) :: allele1       !the vector col with individual's allele 1
!  integer, intent(in) :: allele2       !the vector col with individual's allele 2
!  integer, intent(in)           :: bit_pos
!  !the vector col stating if with individual is genotyped for this locus
!  integer, intent(in ),optional :: genotyped
!
!  integer :: gen
!  logical :: all1,all2,accepted
!
!  if (present(genotyped)) then
!     accepted=btest(genotyped,bit_pos)
!     if  ( accepted ) then
!        all1 = btest(allele1,bit_pos)
!        all2 = btest(allele2,bit_pos)
!        call genotypeBin2Int(gen,all1,all2)
!        genotype = gen
!     else
!        genotype = 0   !genotype missing for this locus
!     endif
!  else
!     all1 = btest(allele1,bit_pos)
!     all2 = btest(allele2,bit_pos)
!     call genotypeBin2Int(gen,all1,all2)
!     genotype = gen
!
!  endif
!end subroutine getAGenotype
!
!
!
!
