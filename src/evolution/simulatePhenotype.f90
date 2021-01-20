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

subroutine SimulatePhenotype(nAnim, nComp, indiv, TBV, locations,&
     ids, phen)

  use constants
  use rng_module
  implicit none

  integer, intent(in) :: nAnim, nComp
  integer, dimension(nanim), intent(in) :: indiv
  double precision, dimension(nanim, ncomp), intent(in) :: TBV
!  type(variance), intent(in) :: varainces
!  double precision, dimension(ncomp) :: means
  double precision, dimension(:,:), intent(in) :: locations
  integer, dimension(:), allocatable, intent(out) :: ids
  double precision, dimension(:), allocatable, intent(out) :: phen

  integer :: i
!  real :: rand
!  real, dimension (:,:), allocatable :: frequency  
!  double precision, dimension (:), allocatable :: tempval, envarray
  !double precision, dimension (:), allocatable :: tbvslo, tbvint
 ! CHARACTER ( LEN = 256 ) :: fileeffect, filetbvs, filetbvi,  filepruningIND
!  double precision :: val1,val2,val3, val4, val5, val6, phenotype, phenotype_slo
!  double precision :: phenotype_int
!  double precision :: h2int, h2sl, sd_A_slo, sd_A_int, sd_E_slo, sd_E_int, sd_PE_slo
!  double precision :: sd_PE_int, sd_P_int, sd_P_sl , env_first, env_last
!  double precision :: Eij_slo, Eij_int
!  double precision :: var_A_slo, var_A_int, var_E_slo, var_E_int, var_PE_slo
!  double precision :: var_PE_int, var_P_slo, var_P_int
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
!  double precision, allocatable, dimension(:) :: temp!, A, Ps, PEs, PEi
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
!
!  var_P_int = var_A_int + var_PE_int + var_E_int
!  var_P_slo = var_A_slo + var_PE_slo + var_E_slo

!  allocate(temp(nanim))
  
!  nparent = nfounder

  if (size(locations) == 1) then
     i = nAnim
  else
     i = size(locations)
  end if
  allocate(phen(i), ids(i))
  phen(1:i) = 0.d0
  ids(1:i) = 0
  !
  !! what I do here is replacing the original tbvslope by its centralised version
  ! scaled by parents sd. this is to make sure tvb's have an expectation zero while
  ! they give a desired value for h2
  
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

end subroutine SimulatePhenotype






