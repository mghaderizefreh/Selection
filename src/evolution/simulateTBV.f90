! Written by Masoud Ghaderi. Documented on 4 Nov, 2019
subroutine SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, genome, QTLlist, SNPlist,&
     TBV, verbose)

  use constants
  implicit none

  integer, intent(in) :: nAnim, nChr, nComp, nSNP
  integer, dimension(nanim), intent(in) :: indiv
  !  type(variance), intent(in) :: varainces
  type(chromosome), dimension(nChr), intent(in) :: genome
  type(QTL_Array), intent(in) :: QTLlist
  integer, dimension(nChr, nSNP) :: SNPlist
  real(KINDR), dimension(:,:), allocatable, intent(out) :: TBV
  logical, intent(in) :: verbose

  integer, dimension (:), allocatable :: haplotype1, haplotype2
  integer, dimension (:), allocatable :: chr_nlocibefore, pruningSNP !, pruningSNP2
  real(KINDR), dimension(:,:), allocatable :: effect
  integer :: i, j, k, id, iloci, a1, a2, nloci, nblock
  integer :: iblck1, ibit1, nbits
  integer :: ichr, totLoci, totQTL
!  real(KINDR) :: val1

  nbits = 32
  allocate(chr_nlocibefore(nChr))
  totQTL=0

  i = 0
  j=0
  k=0
  totLoci=0
  
  do ichr = 1, nchr
     if ( genome(ichr)%nblock > i ) i = genome(ichr)%nblock
     if ( genome(ichr)%nloci  > j ) j = genome(ichr)%nloci
     chr_nlocibefore(ichr) = totLoci
     totLoci               = totLoci + genome(ichr)%nloci
  end do
  allocate(haplotype1(i), haplotype2(i), pruningSNP(totLoci))

  ! allocating the other array 
  allocate(effect(totLoci, nComp))
  allocate(tbv(nanim, nComp))

  effect(1:totLoci, 1:nComp) = 0.D0
  tbv(1:nAnim, 1:nComp) = 0.D0

  pruningSNP(1:totLoci)=0
  ! reading SNP effects and list of SNP chip (different lists)
  do iChr = 1, nChr
     effect(&
          QTLlist%indices(iChr, 1:QTLlist%nQTL) + chr_nlocibefore(ichr),&
          1:nComp) = QTLlist%values(iChr, 1:QTLlist%nQTL, 1:nComp)
     pruningSNP(SNPlist(iChr,1:nSNP) + chr_nlocibefore(ichr)) = 1
  end do
  totQTL = nChr *  QTLlist%nQTL

  ! reading list of inviduals to phenotype
  ! allocate(pruningIND(nanim))
  ! pruningIND(:) = 0
  ! purningIND(i) = 1, if i (index of animal) is to be phenotyped otherwise 0
  ! TODO: get logic from the input indiv

  !================================================
  !now reading the genotype and and calc tbv (if needed)
  !================================================
  id = 0
  do j = 1, nanim
     id = indiv(j)
     do ichr = 1, nchr
        nLoci = genome(iChr)%nLoci
        nBlock = genome(iChr)%nBlock

        ! TODO X: this is not appropriate here as the genome is already read
        !idmissing = chr_idmissing( ichr )

        ! TODO (X): remove or keep based on the X
        !        if ( idmissing == 0 ) then
        !        read ( iungen, *, end = 10 ) haplotype1( 1 : nblock )
        !        else
        !           read ( iungen, *, end = 10 ) i, haplotype1( 1 : nblock )
        !        end if
        !        read ( iungen, *, end = 11 ) haplotype2( 1 : nblock )
        haplotype1(1:nblock) = genome(iChr)%genotypes(id, 1, 1:nblock)
        haplotype2(1:nblock) = genome(iChr)%genotypes(id, 2, 1:nblock)

        iblck1 = 0        !block for loci i
        ibit1 = 0        !bit for loci i
        do iloci = 1, nloci
           ibit1 = ibit1 - 1       ! position of the new loci i to be check
           if ( ibit1 < 0 ) then
              iblck1 = iblck1 + 1          ! a new block
              ibit1 = nbits - 1          ! start from the first bit
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

  if (verbose) write(STDOUT, *) "breeding values simulated"

end subroutine SimulateTBV

