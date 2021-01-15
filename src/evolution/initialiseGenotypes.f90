!==============================================================================================
! this subroutine initialise the genotype for the first generation of the historical generations
!
! there are tree options for setting genotypes (genstart parameter)
! genstart=1  =  all loci set to be fixed to 
! genstart=2  =  all loci are with frequ 0.5, in HWE and in linkage equilibrium (i.e. genotypes are sampled independently to be in HWE )
! genstart=3  =  genotypes are read from files
!
!if option 1 and 2 are used, all chromosomes have similar 
!

subroutine initialiseGenotypes(nchr, nanim, genstart, nloci, nblock, istore, genome,&
     maxloci, maxblock, ifail, prefixfilename)
  use constants
  implicit none

  integer, intent(in) :: genstart, nanim, nchr, istore
  character(len=*), optional, intent(in) :: prefixfilename
  integer, intent(inout) :: nloci
  type(chromosome), DIMENSION(:), intent(out) :: genome
  integer, intent(out) :: ifail, nblock

  integer, dimension(:,:,:), pointer :: genotypes
  integer :: iun, i, k, id, j, ichr, igam, a1, iblock
  real :: rand
  integer :: maxblock, maxloci
  character(len=256) :: filename

  !nloci is input if genstart is 1 or 2
  write(stdout, '(a24, i1)') " initialising genotypes",genstart
  !=====================================================
  ! IF GENSTART=1 or 2
  ! all chrmosmoe have the same number of SNP
  !
  ! the number of loci is passed as input in dummy variable <nloci>
  !=====================================================
  !
  ! This block of code initialise the array of genome TYPE so that
  ! each element has its 'genotypes' pointer array property
  ! allocated to nanim x 2 x nblock
  if (genstart .eq. 1 .or. genstart .eq. 2) then
     if(istore .eq. 1) then
        i = 32
        nblock = nloci/i
        if(mod(nloci,i) > 0) nblock = nblock + 1
     else
        nblock = nloci
     endif
     maxblock = nblock
     maxloci = nloci
     do ichr = 1, nchr
        genome(ichr)%nloci  = nloci
        genome(ichr)%nblock = nblock
        allocate(genome(ichr)%genotypes(nanim, 2, nblock))
     end do
  end if
  ifail = 0

  !=====================================================
  ! GENSTART=1
  ! ALL SNP in ALL CHR SET TO BE  FIXED 
  !
  !=====================================================
  if (genstart .eq. 1) then
     k = 0
     if (istore .eq. 1) then
        do j = 31, 0, - 1
           k = ibclr(k, j ) !ibclr(k,j) returns the value of k  with the bit at position j set to zero.
        end do
     else
        k = 1  !genotype are 1
     end if
     do ichr = 1, nchr
        genome(ichr)%genotypes(:,:,:) = k !all animalss are homozygous for all genes
     end do

     !=====================================================
     ! GENSTART=2
     !  SNP ARE SEGREGATING WITH FREQ =0.5 IN HWE AND IN LINKAGE EQUILIBRIUM
     !
     !=====================================================
  elseif (genstart .eq. 2) then
     ! if the initial genotype are segregating at freq 0.5 and in LE across SNPs
     ! when genotypes are store in bits
     if (istore .eq. 1) then
        do ichr = 1, nchr
           do iblock = 1, nblock
              do igam = 1, 2
                 do id = 1, nanim
                    a1 = k
                    do j = 31, 0, - 1
                       call random_number(rand)
                       !ibset(i,j) returns the value of i with 
                       !the bit at position j set to one.
                       if ( rand > 0.5 ) a1 = ibset( a1, j ) 
                    end do
                    genome(ichr)%genotypes(id, igam, iblock ) = a1
                 end do
              end do
           end do
        end do
     else
        !  when genotypes are store as integer
        a1 = 1
        do ichr = 1, nchr
           do iblock = 1, nblock
              do igam = 1, 2
                 do id = 1, nanim
                    call random_number(rand)
                    !change the allele (it is ok as frequency
                    ! is the same of both alleles)
                    if ( rand > 0.5 ) a1 = 3 - a1  
                    genome(ichr)%genotypes(id, igam, iblock ) = a1
                 end do
              end do
           end do
        end do
     end if

     !=====================================================
     ! GENSTART=3
     ! SNP GENOTYPE ARE READ FROM FILES
     ! CHROMOSOMES CAN HAVE DIFFERENT NUMBER OF SNP/LOCI
     !
     !=====================================================
  elseif (genstart .eq. 3) then

     maxblock = 0
     maxloci = 0
     do ichr = 1, nChr
        write(filename,'(a,i3.3)') trim(prefixfilename), ichr

        write(stdout,*)' reading ', trim(filename)
        open(newunit = iun, file=filename, status = 'old')
        read(iun,*) i, nloci, nblock, k

        write(*,'(a,8i6)') 'ich i nloci nbloc k nanim', ichr, i,&
             nloci, nblock, k, nanim

        if (i .lt. nanim) then
           ifail = 1
           write(stderr,*) &
                ' Error. Genotype file does not have enough individuals',&
                trim(filename)
           write(stderr,*) 'i nanim', i, nanim
           stop 2
        end if
        if (nloci > maxloci) maxloci = nloci
        if (nblock > maxblock) maxblock = nblock

        genome(ichr)%nloci  = nloci
        genome(ichr)%nblock = nblock
        allocate (genome(ichr)%genotypes(nanim, 2, nblock))

        if (k .eq. 0) then
           do id = 1, nanim
              read(iun, *) (genome(ichr)%genotypes(id, 1, j), j = 1, nblock)
              read(iun, *) (genome(ichr)%genotypes(id, 2, j), j = 1, nblock)
           end do
        else
           do id = 1, nanim
              read(iun, *) i, (genome(ichr)%genotypes(id, 1, j), j = 1, nblock)
              read(iun, *)    (genome(ichr)%genotypes(id, 2, j), j = 1, nblock)
           end do
        end if
        close(iun)
     end do
  end if

  write(stdout,'(a)')' end initialising genotypes'

end subroutine initialiseGenotypes


