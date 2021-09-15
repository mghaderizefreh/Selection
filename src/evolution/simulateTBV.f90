! Written by Masoud Ghaderi. Documented on 4 Nov, 2019
subroutine SimulateTBV(nAnim, nChr, nComp, indiv, genome, chr_nlocibefore,&
     QTLlist, TBV, verbose)

  use constants
  implicit none

  integer, intent(in) :: nAnim, nChr, nComp
  integer, dimension(nanim), intent(in) :: indiv
  type(chromosome), dimension(nChr), intent(in) :: genome
  integer, dimension(nchr), intent(inout) :: chr_nlocibefore
  type(QTL_Array), intent(in) :: QTLlist
  real(KINDR), dimension(1:nAnim, 1:nComp), intent(out) :: TBV
  logical, intent(in) :: verbose

  integer, dimension(:), allocatable :: haplotype1, haplotype2
  real(KINDR), dimension(:,:), allocatable :: effect
  integer, save :: totLoci = 0, maxBlock = 0
  integer :: i, j, k, id, iloci, a1, a2, nloci, nblock
  integer :: iblck1, ibit1
  integer :: ichr

  if (totLoci .eq. 0) then
     ! eval totLoci, totQTL, chr_nlocibefore
     do ichr = 1, nchr
        if (genome(ichr)%nblock > maxBlock) maxBlock = genome(ichr)%nblock
        chr_nlocibefore(ichr) = totLoci
        totLoci               = totLoci + genome(ichr)%nloci
     end do
  end if

  call alloc1I(haplotype1, maxBlock, "haplotype1", "simulateTBV")
  call alloc1I(haplotype2, maxBlock, "haplotype2", "simulateTBV")
  call alloc2D(effect, totLoci, nComp, "effect", "simulateTBV")
  effect(1:totLoci, 1:nComp) = ZERO

  ! reading effects
  do iChr = 1, nChr
     effect(QTLlist%indices(iChr, 1:QTLlist%nQTL) +&
          chr_nlocibefore(ichr), 1:nComp) = QTLlist%values(iChr, &
          1:QTLlist%nQTL, 1:nComp)
  end do

  tbv(1:nAnim, 1:nComp) = ZERO

  !================================================
  !now reading the genotype and and calc tbv
  !================================================
  id = 0
  do j = 1, nanim
     id = indiv(j)
     do ichr = 1, nchr
        nLoci = genome(iChr)%nLoci
        nBlock = genome(iChr)%nBlock

        haplotype1(1:nblock) = genome(iChr)%genotypes(id, 1, 1:nblock)
        haplotype2(1:nblock) = genome(iChr)%genotypes(id, 2, 1:nblock)

        iblck1 = 0        !block for loci i
        ibit1 = 0        !bit for loci i
        do iloci = 1, nloci
           ibit1 = ibit1 - 1       ! position of the new loci i to be check
           if ( ibit1 < 0 ) then
              iblck1 = iblck1 + 1          ! a new block
              ibit1 = NBITS - 1          ! start from the first bit
           end if
           a1 = 1 ! first allele is clear bit
           a2 = 1 ! second allele is clear bit
           if (Btest(haplotype1(iblck1), ibit1)) a1 = 2    ! the
           if (Btest(haplotype2(iblck1), ibit1)) a2 = 2

           !the number of SNP cumulative across chrs
           k = chr_nlocibefore(ichr) + iloci  

           i = a1 + a2 - 3  !genotypes score -1,0,1
           tbv(id, 1:nComp) = tbv(id, 1:nComp) + dble(i) * effect(k, 1:nComp)
        end do
     end do
  end do
  if (verbose) then
  end if

  deallocate(haplotype1)
  deallocate(haplotype2)
  deallocate(effect)

end subroutine SimulateTBV
