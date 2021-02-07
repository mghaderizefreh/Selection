subroutine getQTLandSNP(verbose, nChr, nQTL, nSNP, nComp, randomMAF, genome, &
     QTLlist, SNPlist, covMat, baseNameFreq, MAF)
  use constants
  use rng_module
  use quickSort
  implicit none

  logical, intent(in) :: verbose
  integer, intent(in) :: nChr, nQTL, nSNP, nComp
  logical, intent(in) :: randomMAF
  type(chromosome), dimension(:), intent(in) :: genome
  type(QTL_Array), intent(out) :: QTLlist
  integer, dimension(:,:), allocatable, intent(out) :: SNPlist
  real(KINDR), dimension(nComp, nComp), intent(in) :: covMat
  character (len=20), intent(in), optional :: baseNameFreq
  real(KINDR), intent(in), optional :: MAF

  integer :: iChr, i, k, j, iun
  integer :: nReq, nAvail
  logical :: l_exists
  character(len = 256):: fileName
  integer, dimension(:), allocatable :: iReq, temp
  real(KINDR), dimension(:,:), allocatable  :: values_1D
  real(KINDR), dimension(:), allocatable :: means
  type(JArr) , dimension(:) , allocatable:: MAFArray
  real(KINDR) :: rand, freq
  if (.not.randomMAF) then
     if (.not.present(baseNameFreq).or..not.present(MAF)) then
        write(STDERR, *) "Error"
        write(STDERR, *) "baseNameFreq and/or MAF not present"
        write(STDERR, *) "and non-random selection instructed"
        stop 2
     else
        do ichr = 1, nchr
           write(filename, '(a,i3.3)') trim(baseNameFreq), ichr
           inquire(file=filename, exist=l_exists)
           if (.not. l_exists) then
              write(STDERR,'(a4,1x,a,i3.3,1x,a)') "File", trim(baseNameFreq),&
                   ichr,"does not exist"
              write(STDERR, '(a)') "exiting..."
              stop 2
           end if
        end do
     end if
  end if

  allocate(SNPlist(nChr, nSNP))
  allocate(QTLlist%indices(nChr, nQTL), QTLlist%values(nChr, nQTL, nComp))
  allocate(means(ncomp))
  means(1:ncomp) = 0.d0
  k = nQTL * nChr
  allocate(values_1D(nComp, k))
  QTLlist%nComp = nComp
  QTLList%nQTL = nQTL
  nReq = nQTL + nSNP
  allocate(iReq(nReq))
  call gnormal(means, covMat, nComp, k, values_1D)
  if (verbose) write(STDOUT, *) " random values for tbv created"
  if (.not.randomMAF) then
     allocate(MAFArray(nChr))

     ICHRLOOP: do iChr = 1, nChr 
        allocate(MAFArray(iChr)%array(genome(iChr)%nloci))
        write(fileName,'(a,i3.3)') trim(baseNameFreq), ichr
        open(newUnit = iun, file = fileName, status = 'old')
        read(iun, *) 
        do i = 1, genome(iChr)%nloci
           read(iun, *) j, rand, k, freq
           MAFArray(iChr)%array(i) = min(freq, 1- freq)
        end do
        close(iun)

        allocate(temp(genome(iChr)%nloci))
        call sortrx(genome(iChr)%nloci, MAFArray(ichr)%array, temp)

        i = 0
        do while(i < genome(iChr)%nLoci)
           i = i + 1
           if (MAFarray(ichr)%array(temp(i)) .gt. maf) then
              i = i - 1
              exit
           end if
        end do

        nAvail = genome(iChr)%nloci - i

        if (nAvail < nReq) then
681        format(" Number of QTL is less than nLoci for Chromosome ", &
                i2," for MAF > ", f10.8)
682        format(" iChr, NReq, Nloci, available, < maf", i2, 4x, i4, 3x,&
                i4, 2x, i4, 2x, i4)
           write(STDERR, 681) iChr, maf
           write(STDERR, 682) iChr, NReq, genome(iChr)%nLoci, nAvail , i
           write(STDERR, '(a)') "exiting..."
           stop 2
        end if

        temp(1:nAvail) = (/(i, i = 1, nAvail)/)
        call choice(temp, nAvail, nReq, iReq)

        QTLlist%indices(iChr, 1:nQTL) = iReq(1:nQTL)
        do i = 1, nComp
           k = (iChr - 1) * nQTL
           QTLlist%values(iChr, 1:nQTL, i) = values_1D(i, (k + 1) : (k + nQTL))
        end do
        SNPlist(ichr, 1:nSNP) = iReq((nQTL+1):nReq)

        deallocate(temp)
     end do ICHRLOOP

  else ! i.e., when the selection is random

     ICHRLOOP2: do iChr = 1, nChr
        k = genome(iChr)%nLoci
        allocate(temp(k))
        temp = (/(i, i = 1, k)/)
        call choice(temp, k, nReq, iReq)
        QTLlist%indices(iChr, 1:nQTL) = iReq(1:nQTL)
        SNPlist(iChr, 1:nSNP) = iReq((nQTL+1):nReq)
        deallocate(temp)

        do i = 1, nComp
           k = (iChr - 1) * nQTL
           QTLlist%values(iChr, 1:nQTL, i) = values_1D(i, (k+1):(k+nQTL))
        end do

     end do ICHRLOOP2
  end if
  if (verbose) write(STDOUT, *) "QTL and SNP list simulated"
end subroutine getQTLandSNP

