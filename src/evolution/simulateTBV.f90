! Written by Masoud Ghaderi. Documented on 4 Nov, 2019
subroutine SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, genome, QTLlist, SNPlist,&
     TBV, verbose)

  use constants
  implicit none

  integer, intent(in) :: nAnim, nChr, nComp, nSNP
  integer, dimension(nanim), intent(in) :: indiv
  type(chromosome), dimension(nChr), intent(in) :: genome
  type(QTL_Array), intent(in) :: QTLlist
  integer, dimension(nChr, nSNP) :: SNPlist
  real(KINDR), dimension(1:nAnim, 1:nComp), intent(out) :: TBV
  logical, intent(in) :: verbose

  integer, dimension (:), allocatable, save :: haplotype1, haplotype2
  integer, dimension (:), allocatable, save :: chr_nlocibefore, pruningSNP
  real(KINDR), dimension(:,:), allocatable, save :: effect
  integer, save :: totLoci, totQTL
  integer :: i, j, k, id, iloci, a1, a2, nloci, nblock
  integer :: iblck1, ibit1
  integer :: ichr

  j = 0
  k = 0
  
  ! eval totLoci, totQTL, chr_nlocibefore, pruningSNP and effect
  if (.not.allocated(chr_nlocibefore)) then
     allocate(chr_nlocibefore(nChr))
     i = 0
     totLoci = 0
     do ichr = 1, nchr
        if ( genome(ichr)%nblock > i ) i = genome(ichr)%nblock !i=maxblock
        chr_nlocibefore(ichr) = totLoci
        totLoci               = totLoci + genome(ichr)%nloci
     end do
     allocate(effect(totLoci, nComp), pruningSNP(totLoci))
     effect(1:totLoci, 1:nComp) = ZERO
     pruningSNP(1:totLoci) = 0
     ! reading SNP effects and list of SNP chip (different lists)
     do iChr = 1, nChr
        pruningSNP(SNPlist(iChr,1:nSNP) + chr_nlocibefore(ichr)) = 1
     end do
     totQTL = nChr * QTLlist%nQTL
  end if

  do iChr = 1, nChr
     effect(QTLlist%indices(iChr, 1:QTLlist%nQTL) +&
          chr_nlocibefore(ichr), 1:nComp) = QTLlist%values(iChr, &
          1:QTLlist%nQTL, 1:nComp)
  end do

  if (.not.allocated(haplotype1)) allocate(haplotype1(i), haplotype2(i))

  tbv(1:nAnim, 1:nComp) = ZERO

  !================================================
  !now reading the genotype and and calc tbv (if needed)
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

end subroutine SimulateTBV

