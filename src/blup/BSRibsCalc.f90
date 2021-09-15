subroutine BSRibsCalc1(genotypes, amat, nanim, nobs, nSNP, effect, iscaled, ivar,&
     ifail, verbose)
  !=======================================================================
  !
  ! it calculate the relations IBS matrix for additive or dominance effect
  ! given segregating genotypes
  !
  !
  ! it calculates the IBS assuming that mssing values has a genotype score equal to 
  ! the population mean
  !
  !
  !
  !
  ! ifail  = 0    matrix was calculated without problems
  ! ifail /= 0    matrix cannot be calculated
  !        = 1  wrong effect pased as input; it has to be: "add" or "dom"
  !        = 2  no snp has the genotypes required to calculate the effect
  !             additive effect needs two genotypes
  !             dominance       needs the heterozygote and at least one homozygote
  !========================================================================
  use constants
  use omp_lib
  implicit none

  character(len=*), intent(in) :: effect
  integer, intent(in) :: nanim
  integer, intent(in) :: nobs
  integer, intent(in) :: nSNP
  integer, intent(in) :: iscaled
  integer, intent(in) :: ivar
  integer, dimension(nanim, nSNP), intent(in) :: genotypes
  !put as inout so amat is not force to be adjusted if failed
  real(KINDR), dimension(nobs), intent(inout) :: amat  
  integer, intent(out) :: ifail
  logical, intent(in) :: verbose

  integer :: i, j, id1, id2, ipos, isnp
  real(KINDR) :: val1, val4
  integer :: k

  real(KINDR), dimension(:), allocatable :: temparray
  character(len=3) :: theeffect
  integer :: n, usedSNP, igen, ieffect
  !isHomozygote(0:3)
  integer, dimension(0:3) :: ngen
  real(KINDR), dimension(0:3) :: genscore
  real(KINDR), dimension(0:3, 0:3) :: IBSstatus
  real(KINDR) :: sumvar, mean, vari

  !------------------------------------------------------------------------------
  ! calculating the genotype score
  ! depending if matrix is additive or dominance
  !     additive  score = 1,2,3 for aa,ab,bb
  !     dominance score = 0,1,0 for aa,ab,bb
  !------------------------------------------------------------------------------
  ! tranfering effect to theeffect in case the input character variable had more
  ! than 3 characters

  theeffect = effect
  ieffect = 0
  genscore(0) = ZERO
  select case (theeffect)
  case ("ADD", "ADd", "AdD", "aDD", "Add", "aDd", "adD", "add")  
     ! if matrix is additive effect
     genscore(1) = ONE
     genscore(2) = TWO
     genscore(3) = 3._KINDR
     ieffect     = 1
  case ("DOM", "DOm", "DoM", "dOM", "Dom", "dOm", "doM", "dom") 
     ! if matrix is dominance effect
     genscore(1) = ZERO
     genscore(2) = ONE
     genscore(3) = ZERO
     ieffect     = 2
  case default
     write(STDERR, *) " wrong selection of effect"
     write(STDERR, *) &
          " matris is either for addtive (add) or dominance (dom) effect"
     ifail=1
     return
  end select
  if (verbose) write(STDOUT, *) " calculating ", theeffect, ieffect
  if (verbose) write(STDOUT, *) " scaled  ", iscaled
  !------------------------------------------------------------------------------
  ! calculating the IBS matrix
  !------------------------------------------------------------------------------

  ifail = 0

  usedSNP = 0
  sumvar = ZERO
  amat(1:nobs) = ZERO

  !$OMP PARALLEL DEFAULT(NONE)&
  !$OMP PRIVATE(val4, k, temparray)&
  !$OMP SHARED(nobs, nSNP, ieffect, genscore, ivar, iscaled, nanim, genotypes)&
  !$OMP SHARED(sumvar, usedSNP, AMAT)
  call alloc1D(temparray, nobs, "tempArray", "BSRibsCalc1")
  val4 = ZERO
  k = 0
  temparray(1:nobs) = ZERO

  !$OMP DO LASTPRIVATE(i, IBSstatus, ngen, n, mean, vari, val1, j, igen)&
  !$OMP LASTPRIVATE(ipos, id1, id2)
  do isnp = 1, nSNP
     !   counting number of individuals in each genotype class
     do i = 0, 3
        ngen(i) = count(genotypes(1:nanim, isnp) == i)
     end do

     n = ngen(1) + ngen(2) + ngen(3)
     if(n == 0) then
        cycle      !SNP has all genotypes missing (SNP not used in IBS calc)
     end if

     !now checking if SNP can be used in IBS calculation
     i = 0
     if( ngen(2) > 0) i = i + 1
     if(ieffect == 1) then 
        ! if additive matrix it is needed at least two genotype present
        if(ngen(1) > 0) i = i + 1
        if(ngen(3) > 0) i = i + 1
     else
        ! if dominance matrix, it is needed the heterozygote and at least one
        ! homozygote onbserved
        if(ngen(1) > 0 .or. ngen(3) > 0) i = i + 1
     end if
     ! IF(i <= 1) cycle ! the SNP has all genotypes missing or they are fixed to
     !  a given genotype. NO possible to be used in IBS calc
     if(i <= 1) then
        cycle ! the SNP has all genotypes missing or they are fixed to a given 
        ! genotype. NO possible to be used in IBS calc
     end if

     !------------------------------------------------------------------------------
     !   if reaching here means that the SNP can be used
     !  calculating mean and variance of genotype score  for SNP isnp
     !  and ibsstatus
     !------------------------------------------------------------------------------
     mean = sum(dble(ngen(1:3)) * genscore(1:3)) ! sum of genscores = 
     ! nAA*efAA + nAb*efAB + nBB*efBB
     vari = sum(dble(ngen(1:3)) * genscore(1:3) * genscore(1:3)) ! sum of squared
     vari = (vari - mean * mean / dble(n)) / dble(n) ! variance
     mean = mean / dble(n) !mean
     if(ivar == 1) then

        ! variance of additive  effect assuming HWE = 2pq
        vari = (TWO * ((mean - ONE) / TWO) * (ONE - (mean - ONE) / TWO))
        ! variance of dominance effect assuming HWE = (2pq)*(2pq)
        if (ieffect == 2) then
           vari = vari * vari 
        end if
     end if
     
     val4 = val4 + vari
     k = k + 1    ! SNP to be used in calculation

     !------------------------------------------------------------------------------
     ! calculating IBS status for pairs with all possible genotype score of ind in
     !   pair because there are only 9 possible combination (aa,ab,bb X aa,ab,bb)
     !   calculating the IBS status once and later reading from table should
     !   speed-up calculation
     !------------------------------------------------------------------------------
     val1 = ONE
     if(iscaled == 1) val1 = vari   !genotype score to be scaled so variance is 1
     IBSstatus(0,0) = ZERO
     do i = 1, 3
        IBSstatus(0,i) = ZERO
        IBSstatus(i,0) = ZERO
        do j = 1, 3
           IBSstatus(i,j) = (genscore(i) - mean) * (genscore(j) - mean) / val1
        end do
     end do
     !------------------------------------------------------------------------------
     ! now do the ibs calculation for all ind with genotype
     ! calculation is sum of all IBSstatus given genotype score acroos all valid
     ! SNP
     !------------------------------------------------------------------------------

     do id1 = 1, nanim
        igen = genotypes(id1, isnp)
        if(igen == 0) cycle ! genotype missing,  not need to estimate IBS
        do id2 = 1, id1
           ipos = (id1 - 1) * id1 / 2 + id2
           temparray(ipos) = temparray(ipos) + IBSstatus(igen, genotypes(id2,isnp))
        end do
     end do
     !------------------------------------------------------------------------------
  end do
  !$OMP END DO
  !$OMP CRITICAL
  usedSNP = usedSNP + k
  sumvar = sumvar + val4
  amat(1:nobs) = amat(1:nobs) + temparray(1:nobs)
  deallocate(temparray)
  !$OMP END CRITICAL
  !$OMP END PARALLEL
  !------------------------------------------------------------------------------
  ! now divide the matrix by a denominator
  ! denominator is sum(genscore var)
  !       if genscore was scaled, then this sum(var) is the number of SNP used in 
  !       calculation
  !------------------------------------------------------------------------------

  if (usedSNP > 1) then

     ! if genscore was scaled then denominator is number fo snp used
     val4 = dble(usedSNP) 
     if(iscaled == 0) val4 = sumvar ! if not,scaled then denominator is sum variances
     i = nanim * (nanim + 1) / 2
     amat(1:i) = amat(1:i) / val4
     if (verbose) WRITE(STDOUT, *) " number of SNP used ", usedSNP,nSNP
  else
     ifail = 2   !no SNP was segregating to get a matrix calculated
  end if

end subroutine BSRibsCalc1

!=======================================================================
subroutine BSRibsCalc1a(genotypes,amat, nanim, nobs, nSNP,effect,iscaled, ivar,&
     usedSNPMAT, cumvarMAT, ifail, verbose)
  !=======================================================================
  !
  ! it calculates the relations IBS matrix for additive or dominance effect
  ! given segregating genotypes
  !
  !
  !  IBS for a given pair are calculated only with know genotyped 
  !
  !
  !
  !
  ! ifail  = 0    matrix was calculated without problems
  ! ifail /= 0    matrix cannot be calculated
  !        = 1  wrong effect pased as input; it has to be: "add" or "dom"
  !        = 2  no snp has the genotypes required to calculate the effect
  !             additive effect needs two genotypes
  !             dominance       needs the heterozygote and at least one homozygote
  !========================================================================
  use constants
  implicit none

  character(len=*), intent(in) :: effect
  integer, intent(in) :: nanim
  integer, intent(in) :: nobs
  integer, intent(in) :: nSNP
  integer, intent(in) :: iscaled
  integer, intent(in) :: ivar
  integer, dimension(nanim,nSNP), intent(in) :: genotypes
  !put as inout so amat is not forced to be adjusted if failed
  real(KINDR), dimension(nobs), intent(inout) :: amat
  integer, intent(out):: ifail
  logical, intent(in) :: verbose

  !put as inout so cumvarMAT is not forced to be adjusted if failed
  real(KINDR), dimension(:), intent(inout) :: cumvarMAT  
  !put as inout so usedSNPMAT is not forced to be adjusted if failed
  integer, dimension(:), intent(inout) :: usedSNPMAT 
  integer :: i, j, id1, id2, ipos, isnp
  real(KINDR):: val1, val4
  integer :: k

  character(len=3) :: theeffect
  integer :: n, usedSNP, igen, ieffect
  integer, dimension(0:3) :: ngen
  real(KINDR), dimension(0:3) :: genscore
  real(KINDR), dimension(0:3,0:3) :: IBSstatus
  real(KINDR) :: sumvar, mean, vari

  real(KINDR), dimension(0:3) :: IBSvar
  integer, dimension(0:3) :: IBSnSNP

  integer, dimension(nanim), automatic :: idused

  !------------------------------------------------------------------------------
  ! calculating the genotype score
  ! depending if matrix is additive or dominance
  !     additive  score = 1,2,3 for aa,ab,bb
  !     dominance score = 0,1,0 for aa,ab,bb
  !------------------------------------------------------------------------------
  ! tranfering effect to theeffect in case the input character variable had more 
  ! than 3 characters
  theeffect = effect  
  ieffect = 0
  genscore(0) = ZERO
  select case (theeffect)
  case ("ADD", "ADd", "AdD", "aDD", "Add", "aDd", "adD", "add")
     ! if matrix is additive effect
     genscore(1) = ONE
     genscore(2) = TWO
     genscore(3) = 3._KINDR
     ieffect     = 1
  case ("DOM", "DOm", "DoM", "dOM", "Dom", "dOm", "doM", "dom")
     ! if matrix is dominance effect
     genscore(1) = ZERO
     genscore(2) = ONE
     genscore(3) = ZERO
     ieffect     = 2
  case default
     write(STDERR, *) " wrong selection of effect"
     write(STDERR, *) &
          " matris is either for addtive (add) or dominance (dom) effect"
     ifail = 1
     return
  end select
  if (verbose) write(STDOUT, *) " calculating ", theeffect, ieffect
  if (verbose) write(STDOUT, *) " scaled  ", iscaled

  !------------------------------------------------------------------------------
  ! setting to zero the IBS status for a pair where at least one has missing genotype
  ! genscore for a missing genotype is the pop mean (or zero after centering).
  ! hence since IBD is the crosproduct of genscore, the IBS of this type of pair is 0
  !------------------------------------------------------------------------------
  IBSvar(0:3) = ZERO
  IBSnSNP(0) = 0
  IBSnSNP(1:3) = 1

  IBSstatus(0, 0) = ZERO
  do i = 1, 3
     IBSstatus(0,i) = ZERO
     IBSstatus(i,0) = ZERO
  end do
  !------------------------------------------------------------------------------
  ! calculating the IBS matrix
  !------------------------------------------------------------------------------
  k = nSNP / 20
  if(k == 0) k = 1
  ifail = 0
  ipos = nanim * (nanim + 1) / 2
  j = size(amat, dim = 1)

  amat(1:ipos) = ZERO
  cumvarMAT(1:ipos) = ZERO
  usedSNPMAT(1:ipos) = 0

  sumvar = ZERO
  usedSNP = 0

  do isnp = 1, nSNP
!     if(mod(isnp, k ) == 1) write(*, '(a4,4i9,f25.8 )') 'at ', isnp, usedSNP, nSNP,&
!          k, sumvar

     !   counting number of individuals in each genotype class
     do i = 0, 3
        ngen(i) = count(genotypes(1:nanim, isnp) == i)
     end do
     n = ngen(1) + ngen(2) + ngen(3)
     if(n == 0) then
        if (verbose) write(STDOUT, '(a,6i9)') " snp not used", isnp, ngen(0:3), n
        cycle   !SNP has all genotypes missing (SNP not used in IBS calc)
     end if
     !now checking if SNP can be used in IBS calculation
     i = 0
     if(ngen(2) > 0) i = i + 1

     ! if additive ma2trix it is needed at least two genotype present
     if(ieffect == 1) then
        if(ngen(1) > 0) i = i + 1
        if(ngen(3) > 0) i = i + 1

        ! if dominance matrix, it is needed the heterozygote and at least one 
        ! homozygote onbserved
     else
        if(ngen(1) > 0 .or. ngen(3) > 0) i = i + 1
     end if
     if(i <= 1) then
        if (verbose) write(STDOUT,'(a,6i9)')' snp not used', isnp, ngen(0:3), n
        cycle ! the SNP has all genotypes missing or they are fixed to a given 
        ! genotype. NO possible to be used in IBS calc
     end if

     !------------------------------------------------------------------------------
     !   if reaching here means that the SNP can be used 
     !------------------------------------------------------------------------------

     !------------------------------------------------------------------------------
     !  calculating mean and variance of genotype score  for SNP isnp
     !  and ibsstatus
     !------------------------------------------------------------------------------
     mean = sum(dble(ngen(1:3)) * genscore(1:3)) ! sum of genscores 
     ! = nAA*efAA + nAb*efAB + nBB*efBB
     vari = sum(dble(ngen(1:3)) * genscore(1:3) * genscore(1:3))  ! sum of squared
     vari = (vari - mean * mean / dble(n)) / dble(n) !variance
     mean = mean / dble(n) !mean (assuming genescore 1,2,3 for aa,ab, bb)
     if (ivar == 1) then
        ! variance of additive  effect assuming HWE = 2pq
        vari = (TWO * ((mean - ONE) / TWO) * (ONE - (mean - ONE) / TWO)) 
        ! variance of dominance effect assuming HWE = (2pq)*(2pq)
        if (ieffect == 2) vari = vari * vari 
     end if

     sumvar = sumvar + vari
     usedSNP = usedSNP + 1    ! SNP to be used in calculation

     !------------------------------------------------------------------------------
     ! calculating IBS status for pairs with all possible genotype score of ind in 
     !   pair because there are only 9 possible combination (aa,ab,bb X aa,ab,bb)
     !   calculating the IBS status once and later reading from table should
     !   speed-up calculation
     !------------------------------------------------------------------------------
     val1 = ONE
     if(iscaled==1) val1 = vari   !genotype score to be scaled so variance is 1
     do i = 1, 3
        do j = 1, 3
           IBSstatus(i,j) = (genscore(i) - mean) * (genscore(j) - mean) / val1
        end do
     end do

     IBSvar(1:3) = vari
     !------------------------------------------------------------------------------
     ! now do the ibs calculation for all ind with genotype
     ! calculation is just sum of all IBSstatus given genotype score acroos all valid
     ! SNP
     !
     ! matrix is symmetric so only half store (lower diag row)
     !------------------------------------------------------------------------------
     do id1 = 1, nanim
        idused(id1) = genotypes(id1, isnp)
     end do
     ipos = 0
     do id1 = 1, nanim
        igen = idused(id1)

        if(igen == 0) cycle ! genotype missing,  not need to estimate IBS

        do id2 = 1, id1
           ! the calculation of next position needs to done this way as 
           ! there is a cyle event
           ipos = (id1 - 1) * id1 / 2 + id2  
           amat(ipos) = amat(ipos) + IBSstatus(igen, idused(id2))

           cumvarMAT(ipos) = cumvarMAT(ipos) + IBSvar(idused(id2))
           ! if genotype of id2 is missing then cumulation is zero, otherwise 1
           usedSNPMAT(ipos) = usedSNPMAT(ipos) + IBSnSNP(idused(id2))
        end do
     end do
     !------------------------------------------------------------------------------
  end do

  !------------------------------------------------------------------------------
  ! now divide the matrix by a denominator
  ! denominator is sum(genscore var)
  !       if genscore was scaled, then this sum(var) is the number of SNP used in 
  !       calculation
  !------------------------------------------------------------------------------

  if(usedSNP > 1) then
     ! if genscore was scaled then denominator is number fo snp used
     val4 = dble(usedSNP)
     if (iscaled == 0) val4 = sumvar!if not,scaled then denominator is sum variances
     i = nanim * (nanim + 1) / 2
     ! total normalised by sum of variances
     if (iscaled == 0) amat(1:i) = amat(1:i) / cumvarMAT(1:i)
     ! total normalised by number of SNPused
     if (iscaled /= 0) amat(1:i) = amat(1:i) / usedSNPMAT(1:i)

     if (verbose) write(STDOUT, *) " number of SNP used ", usedSNP, nSNP
  else
     ifail = 2   !no SNP was segregating to get a matrix calculated
  end if

  return

end subroutine BSRibsCalc1a

!=======================================================================
