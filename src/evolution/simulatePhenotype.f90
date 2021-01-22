subroutine SimulatePhenotype(nAnim, nComp, indiv, TBV, variances, means, &
     locations, ids, phen, proc, cte)

  use constants
  use rng_module
  implicit none

  integer, intent(in) :: nAnim, nComp
  integer, dimension(nanim), intent(in) :: indiv
  double precision, dimension(nanim, ncomp), intent(in) :: TBV
  type(variance), intent(in) :: variances
  double precision, dimension(ncomp) :: means
  double precision, dimension(:,:), intent(in) :: locations
  integer, dimension(:), allocatable, intent(out) :: ids
  double precision, dimension(:,:), allocatable, intent(out) :: phen
  character(len = *) :: proc
  double precision, intent(out) :: cte

  integer :: i, j, k
  double precision, dimension (:),allocatable :: temp
  real, dimension(:,:), allocatable :: temp2
  real, dimension (:), allocatable :: tempr
  double precision, allocatable, dimension(:,:) :: A, E!, PE

  allocate(temp2(nComp, nComp))
  temp2(1:nComp, 1:nComp) = 0.0
  do i = 1, nComp
     temp2(i,i) = real(variances%E(i))
  end do
  tempr(1:nComp) = 0.0

  j = size(locations)
  if (j == 1) then
     i = nAnim
     allocate(phen(i,1), ids(i), E(i,nComp))
     ids = indiv
  elseif ((j < nAnim) .and. (size(locations, 1) .eq. 1)) then
     i = j * nAnim
     allocate(phen(i,2), ids(i))
     do k = 1, j
        phen(k:i:j, 1) = locations(1, k)
        ids(k:i:j) = indiv(1:nAnim)
     end do
  elseif (size(locations, 1) == nAnim) then
     i = size(locations, 2)
     allocate(phen(j,2), ids(j))
     write(6, *) "i", i
     write(6, *) "phen", shape(phen)
     do k = 1, i
        write(6, *) 'k, phen to hold', k, shape(phen(k:j:i,1))
        phen(k:j:i, 1) = locations(1:nAnim, k)
        ids( ((k - 1) * nAnim + 1) : (k * nAnim) ) = indiv(1:nAnim)
     end do
  else
     write(STDERR, *) "Error: The format of 'locations' is wrong"
     write(STDERR, *) "Use either:"
     write(STDERR, *) " - [size:    1 x 1   ] all animals in the same location"
     write(STDERR, *) " - [size:    1 x nLox] all have phenotypes in nlox places"
     write(STDERR, *) " - [size:nAnim x nLox] phenotypes at different places"
  end if

  call gnormal(tempr, temp2, nComp, i, E)
  ! todo: implementation for PE if really it is required
  
  select case (proc(1:3))
  case ("COV", "COv", "CoV", "Cov", "cOV", "cOv", "coV", "cov")
     call covariate(nComp, nAnim, indiv, TBV, E, phen, locations, means, cte)
  case default
     write(STDERR, *) "error:"
     write(STDERR, *) "case '", proc, "' not implemented"
     stop 2
  end select
     

  open(1, file = 'phentest')
  do i = 1, size(phen,1)
     write(1, *) ids(i), phen(i,:)
  end do
  close(1)
  stop
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






