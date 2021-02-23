subroutine getGmatrix(nanim, nChr, nSNP, ident, genome, SNPlist, iscaled, &
     ivar, imiss, addDom, Amat, verbose)
  use constants
  implicit none

  logical, intent(in) :: verbose
  integer, intent(in) :: nanim, nChr, nSNP !(nSNP per chromosome)
  integer, dimension(nanim), intent(in) :: ident
  type(chromosome), dimension(nChr), intent(in) :: genome
  integer, dimension(nChr,nSNP), intent(in) :: SNPlist
  integer, intent(in) :: iscaled !(0:no, 1:yes)
  integer, intent(in) :: ivar !(0:sample, 1:2pq, 2:2p'q')
  integer, intent(in) :: imiss !(0:mean, 1:ignore)
  integer, intent(in) :: addDom !(1:additive, 2:dominance)
  real(KINDR), dimension(:), allocatable :: Amat

  integer :: totalSNP!, maxLength = 1
  integer :: i, k
  integer :: a1, a2
  integer :: totLoci, ibit1, iblck1, iloci, id, iChr, ipos, isnp
  integer, dimension(:,:), allocatable :: genotype
  integer, dimension(:), allocatable :: chr_nlocibefore, pruningSNP, usedSNPMat
  real(KINDR), dimension(:), allocatable :: cumvarMat
  character(len=3) :: effect
  integer, dimension(2,2) :: genocode
  integer, dimension(0:3) :: ngeno
  integer :: ifail
!  character(len = 100) :: formato

  external :: bsribscalc1, bsribscalc1a

  genocode(1,1) = 1
  genocode(1,2) = 2
  genocode(2,1) = 2
  genocode(2,2) = 3

  if (ident(1) > 0) then
  end if
  !--------------------------------------------
  ! reading genome to create the genotype
  !--------------------------------------------
  allocate(chr_nlocibefore(nChr))
  totLoci = 0
  do iChr = 1, nChr
!     if (genome(iChr)%nblock > i) i = genome(iChr)%nblock
!     if (genome(iChr)%nloci > j) j = genome(iChr)%nloci
     chr_nlocibefore(iChr) = totLoci
     totLoci = totLoci + genome(iChr)%nloci
  end do
  allocate(pruningSNP(totLoci))
  pruningSNP(1:totLoci) = 0
  totalSNP = nSNP * nChr
  allocate(genotype(nanim, totalSNP))

  do id = 1, nanim
     k = 0
     do iChr = 1, nChr

        pruningSNP(chr_nlocibefore(iChr) + SNPlist(iChr,:)) = 1

        i = 1
        iblck1 = 0 !block for loci i
        ibit1 = 0  ! bit for loci i
        do iloci = 1, genome(iChr)%nloci
           ibit1 = ibit1 - 1 ! position of the new loci i to be checked

           if (ibit1 < 0) then
              iblck1 = iblck1 + 1 ! a new block
              ibit1 = NBITS - 1 ! start from the 1st bit
           end if

           if (pruningSNP(chr_nlocibefore(iChr) + iloci) .eq. 0) cycle

           k = k + 1
           a1 = 1
           a2 = 1
           if (btest(genome(iChr)%genotypes(id, 1, iblck1), ibit1)) a1 = 2
           if (btest(genome(iChr)%genotypes(id, 2, iblck1), ibit1)) a2 = 2
           genotype(id, k) = genocode(a1, a2)
        end do
     end do
  end do
  !====================================
  !------------------------------------
  ! counting non-segregating snp?
  !------------------------------------
  ipos=0    ! counting non-segregating snp?
  do isnp = 1, totalSNP
     ngeno(0:3) = 0
     do i=1, nanim
        ngeno(genotype(i,isnp)) = ngeno(genotype(i,isnp)) + 1
     end do
     i = 0
     if (ngeno(1) > 0) i = i + 1
     if (ngeno(2) > 0) i = i + 1
     if (ngeno(3) > 0) i = i + 1
     if (i <= 1) then
        if (verbose) write(STDOUT, *) " snp to eliminate", isnp, ngeno
        ipos = ipos+1
     endif
  end do
  if (verbose) write(STDOUT, *) " number of SNP non usable", ipos, nsnp
  !====================================
  effect = "add"
  if (addDom == 2) effect = "dom"
  if (verbose) write(STDOUT, '(a11,a3)') "effect is ", effect

  i = nAnim * (nAnim + 1) / 2   
  if (.not.allocated(amat)) allocate(amat(i))

  if (imiss == 1) then
     allocate(usedSNPMAT(i))
     allocate( cumvarMAT(i))
  endif

  ifail=0  
  if (verbose) write(STDOUT, '(a)') " starting ibs cal"

  if (imiss == 0) then
     call bsribscalc1(genotype, amat, nanim, i, totalSNP, effect, iscaled,&
          ivar, ifail, verbose)
  else
     call bsribscalc1a(genotype, amat, nanim, i, totalSNP, effect, iscaled,&
          ivar, usedSNPMAT, cumvarMAT, ifail, verbose)
  end if

end subroutine getGmatrix
