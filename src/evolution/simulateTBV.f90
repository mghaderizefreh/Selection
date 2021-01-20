!    Creates the phenotype file for intercept and slope given variances
!
!    Inputs are genetic variance for slope and intercept, 
!                    variance for permanent env effect for slope and intercetp (can be zero)
!                    environment variance for slope and intercept
!                    number of equidistance intervals in the environment
!                    first environment (double precision)
!                    last environemnt (double precision)
!                    filename of the id of the inidividuals which are going to be phenotyped
!                    filename of the id of the founders (variances are based on the founders)
!                    name of a file for storing true breeding values for intercept __ see makePhenotype.f90 for output type
!                    name of a file for stroing true breeding values for slope __ see makePhenotype.f90 for output type
!
!    The assumptions are:
!                    individuals are genotyped in files with name `pedigree.ch???` where ??? is chromosome number
!                    there are 26 chromosomes
!                    QTL values (for slope and intercept) are stored in QTLlist.txt
!                    SNP list are stored in SNPlist.txt
!
!    Outputs are stored as:
!                    phen.slo as the phenotype of slope with format ``i, Phen_slo(i), TBV_slo(i), PE_slo(i), E_slo(i,j)''
!                                    where i is individual and j is the environment
!                    phen.inc as phenotype of intercept similar to phen.slo
!                    phen.is as phenotype of individuals as ``i, env(j), Phen(i,j)'' with j being nested in i
!
! Written by Masoud Ghaderi. Documented on 4 Nov, 2019

subroutine SimulateTBV(nAnim, nChr, nComp, nSNP, indiv, genome, QTLlist, SNPlist,&
     TBV)

  use constants
  implicit none

  integer, intent(in) :: nAnim, nChr, nComp, nSNP
  integer, dimension(nanim), intent(in) :: indiv
  !  type(variance), intent(in) :: varainces
  type(chromosome), dimension(nChr), intent(in) :: genome
  type(QTL_Array), intent(in) :: QTLlist
  integer, dimension(nChr, nSNP) :: SNPlist
  double precision, dimension(:,:), allocatable, intent(out) :: TBV
  !  double precision, dimension(nAnim), intent(out) :: phen


  integer, dimension (:), allocatable :: haplotype1, haplotype2
  integer, dimension (:), allocatable :: pruningSNP, chr_nlocibefore
  double precision, dimension(:,:), allocatable :: effect
  integer :: i, j, k, id, iloci, a1, a2, nloci, nblock
  integer :: iblck1, ibit1, nbits
  integer :: ichr, totLoci, totQTL

!  real :: rand
!  real, dimension (:,:), allocatable :: frequency  
!  double precision, dimension (:), allocatable :: slopeeffect, inteffect, envarray
  !double precision, dimension (:), allocatable :: tbvslo, tbvint
 ! CHARACTER ( LEN = 256 ) :: fileeffect, filetbvs, filetbvi,  filepruningIND
!  double precision :: val1,val2,val3, val4, val5, val6, phenotype, phenotype_slo
!  double precision ::  phenotype_int
!  double precision :: h2int, h2sl, sd_A_slo, sd_A_int, sd_E_slo, sd_E_int, sd_PE_slo
!  double precision :: sd_PE_int, sd_P_int, sd_P_sl , env_first, env_last
!  double precision :: muint, muslo
!  double precision :: var_A_slo, var_A_int, var_E_slo, var_E_int, var_PE_slo
!  double precision :: var_PE_int, var_P_slo, var_P_int, Eij_slo, Eij_int
!  CHARACTER ( LEN = 30 ) :: baseName, startfile
!  INTEGER :: i, j,k, id, iloci, a1, a2, ihap, offfileun
!  integer :: iblck1, iblck2, ibit1, ibit2, nbits, nseed
!  integer :: nQTL, totQTL, nsteps, ifile, sfile, nparent, nOff, nFounder
!  integer :: iungen ,iuneff,iuntbvs, iuntbvi,iunpruneSNP,iunpruneIND, iunps,iunpi
!  integer :: iunpo, iunpOff
!  integer :: idmissing, nloci, nblock, ichr, totLoci,totsavedLoci, totSavedInd
!  integer, dimension(:), allocatable :: chr_nblock, chr_nloci, chr_nlocibefore
!  character (len=256) :: filepruneSNP, fileQTL
!  real, DIMENSION ( :, :, : ), allocatable :: chr_frequency
!  integer, dimension ( : ), allocatable :: pruningIND, pruningSNP
!  double precision, allocatable, dimension(:) :: A, Ps, PEs, PEi, temp
!  logical :: allexist, L_EXISTS, forOff

  !  call askInteger(nFounder, " input the number of founders")
  !  call askInteger(nsteps, ' Enter the number of equidistance intervals')
  !  read *, env_first
  !  read *, env_last
  !  var_A_slo = variances%A(1)
  !  var_A_int = variances%A(2)
  !  var_PE_slo = variances%PE(1)
  !  var_PE_int = variances%PE(2)
  !  var_E_slo = variances%E(1)
  !  var_E_int = variances%E(2)

  !  var_P_int = var_A_int + var_PE_int + var_E_int
  !  var_P_slo = var_A_slo + var_PE_slo + var_E_slo

  !  muslo = means(1)
  !  muint = means(2)

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
  pruningSNP( : ) = 1   !all SNP has been selected to be saved

  ! allocating the other array 
  allocate(effect(totLoci, nComp))
  write(6, *) 'here'
  allocate(tbv(nanim, nComp))

  effect(1:totLoci, 1:nComp) = 0.D0
  tbv(1:nAnim, 1:nComp) = 0.D0

  ! reading SNP effects
  ! TODO : improve this loop
  do iChr = 1, nChr
     do j = 1, QTLlist%nQTL
        ! add the cumulative value
        i = QTLlist%indices(iChr, j) + chr_nlocibefore(ichr)  
        effect(i, 1:nComp) = QTLlist%values(iChr, j, 1:nComp)
        !intslope(i) = val1
        !inteffect(i) = val2
     end do
  end do
  totQTL = nChr *  QTLlist%nQTL


  ! reading list of SNP
  ! TODO: improve this loop
  pruningSNP(:)=0
  do iChr = 1, nChr
     do j = 1, nSNP
        i = SNPlist(iChr,j) + chr_nlocibefore(ichr)  ! add the cumulative value
        pruningSNP(i) = 1
     end do
  end do

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
        ! TODO X: do I need to implement missing id (ask Ricardo)
        !        idmissing = chr_idmissing( ichr )

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

  ! making phenotypes
  !  write(*,*)' writing tbv(s)' ! Literally duplicating Ricardo's code 
  !  !(i.e. centralising and scaling)
  !  val1=0.d0  
  !  val2=0.d0
  !  val4 = 0.d0
  !  val5 = 0.d0
  !  do i=1,nanim
  !     val1 = val1 + tbvslo(i)
  !     val4 = val4 + tbvint(i) 
  !     val2 = val2 + tbvslo(i) * tbvslo(i)
  !     val5 = val5 + tbvint(i) * tbvint(i)
  !  end do
  !  val3 = val1 / dble(nanim)
  !  val6 = val4 / dble(nanim)
  !  val2 = (val2 - (val1 * val1) / dble(nanim)) / dble(nanim - 1)
  !  val5 = (val5 - (val4 * val4) / dble(nanim)) / dble(nanim - 1)
  !  write ( iuntbvs, '(a82,3g25.16)') &
  !       'id tbv_standarised_to_var1 tbv_meanzero tbv_raw# mean var sd=', &
  !       val3, val2, sqrt(val2)
  !  write ( iuntbvi, '(a82,3g25.16)') &
  !       'id tbv_standarised_to_var4 tbv_meanzero tbv_raw# mean var sd=', &
  !       val6, val5, sqrt(val5)
  !  val2 = sqrt(val2)
  !  val5 = sqrt(val5)
  !  do i = 1, nanim
  !     write(iuntbvs, '(i12,3g25.16)') i, (tbvslo(i) - val3) / val2, &
  !          tbvslo(i) - val3, tbvslo(i)
  !     write(iuntbvi, "(i12,3g25.16)") i, (tbvint(i) - val6) / val5, &
  !          tbvint(i) - val6, tbvint(i)
  !  end do
  !  close(iuntbvs)
  !  close(iuntbvi)
  !
  !  allocate(temp(nanim))
  !
  !  nparent = nfounder
  !
  !  !! what I do here is replacing the original tbvslope by its centralised version
  !  ! scaled by parents sd. this is to make sure tvb's have an expectation zero while
  !  ! they give a desired value for h2
  !
  !  val1 = sum(tbvslo(1 : nparent))
  !  temp = tbvslo * tbvslo
  !  val2 = sum(temp(1 : nparent))
  !  sd_A_slo = sqrt( (val2 - val1 * val1 / dble(nparent)) / dble(nparent - 1))
  !
  !  tbvslo = (tbvslo - val1 / dble(nparent)) / sd_A_slo * sqrt(var_A_slo)
  !  sd_E_slo = sqrt(var_E_slo)!  sd_E_slo = sqrt(1 - h2sl) * sd_P_sl
  !
  !  val4 = sum(tbvint(1 : nparent))
  !  temp = tbvint * tbvint
  !  val5 = sum(temp(1 : nparent))
  !  sd_A_int = sqrt( (val5 - val4 * val4 / dble(nparent)) / dble(nparent - 1))
  !  tbvint = (tbvint - val4 / dble(nparent)) / sd_A_int * sqrt(var_A_int)
  !  sd_E_int = sqrt(var_E_int) ! sqrt(1 - h2int) * sd_P_int
  !
  !  ! permanent and random effect to be added later in the loop
  !  ! Ps = musl + tbvslope  
  !  ! permanent and random effect to be added later in the loop
  !  ! Pi = muint + tbvint   
  !
  !  allocate(envarray(nsteps))
  !
  !  envarray(1:nsteps-1) = (/ (env_first + i * (env_last - env_first)/&
  !       (nsteps - 1), i = 0, nsteps - 2) /)
  !  envarray(nsteps) = env_last
  !
  !  sd_E_slo = sqrt(var_E_slo)
  !  sd_E_int = sqrt(var_E_int)
  !  sd_PE_slo = sqrt(var_PE_slo)
  !  sd_PE_int = sqrt(var_PE_int)
  !
  !  do i = 1, nanim
  !     forOff = (pruningInd(i) == 1)
  !     call normdev(nseed, rand)
  !     PEs(i) = dble(rand) * sd_PE_slo
  !     call normdev(nseed, rand)
  !     PEi(i) = dble(rand) * sd_PE_int
  !     do j = 1, nsteps
  !        call normdev(nseed, rand)
  !        Eij_int = dble(rand) * sd_E_int
  !        call normdev(nseed, rand)
  !        Eij_slo = dble(rand) * sd_E_slo
  !
  !        phenotype_int = muint + tbvint(i) + PEi(i) + Eij_int
  !        phenotype_slo = muslo + tbvslo(i) + PEs(i) + Eij_slo
  !        phenotype = phenotype_int + phenotype_slo * EnvArray(j)
  !
  !        write(iunps, 110) i, phenotype_slo, tbvslo(i), PEs(i), Eij_slo
  !        write(iunpi, 110) i, phenotype_int, tbvint(i), PEi(i), Eij_int
  !        if (forOff) then
  !           write(iunpoff, *) i, EnvArray(j), phenotype 
  !        else
  !           write(iunpo, 120) i, EnvArray(j), phenotype 
  !        end if
  !     end do
  !  end do
  !
  !  close(iunps)
  !  close(iunpi)
  !  close(iunpo)
  !  close(iunpoff)
  !
  !110 format(i12, 1x, 4g24.15)
  !120 format(i12, 1x, 2g24.15)

end subroutine SimulateTBV






