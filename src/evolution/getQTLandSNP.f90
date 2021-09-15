subroutine getQTLandSNP(verbose, nChr, nQTL, nSNP, nComp, randomMAF, genome, &
     QTLlist, SNPlist, covMat, baseNameFreq, MAF)
  use constants
  use rng_module
  use quickSort
  implicit none

  logical, intent(in) :: verbose
  integer, intent(in) :: nChr, nQTL, nSNP, nComp
  logical, intent(in) :: randomMAF
  type(chromosome), dimension(1:nChr), intent(in) :: genome
  type(QTL_Array), intent(inout) :: QTLlist
  integer, dimension(nChr,nSNP), intent(out) :: SNPlist
  real(KINDR), dimension(nComp, nComp), intent(in) :: covMat
  character (len=10), intent(in), optional :: baseNameFreq
  real(KINDR), intent(in), optional :: MAF

  integer :: iChr, i, k, j, iun
  integer :: nReq, nAvail, maxloci
  logical :: l_exists
  character(len = 256):: fileName
  integer, dimension(:), allocatable :: iReq, temp
  real(KINDR), dimension(:,:), allocatable  :: values_1D
  real(KINDR), dimension(ncomp) :: means
  type(JArr) , dimension(:) , allocatable:: MAFArray
  real(KINDR) :: rand, freq
  if (.not.randomMAF) then
     if (.not.present(baseNameFreq).or..not.present(MAF)) then
        write(STDERR, '(a)') "Error"
        write(STDERR, *) "baseNameFreq and/or MAF not present"
        write(STDERR, *) "and non-random selection instructed"
        stop 2
     else
        do ichr = 1, nchr
           write(filename, '(a,i3.3)') trim(baseNameFreq), ichr
           inquire(file=filename, exist=l_exists)
           if (.not. l_exists) then
              write(STDERR, '(a)') "Error"
              write(STDERR,'(a5,1x,a,i3.3,1x,a)') " File", trim(baseNameFreq),&
                   ichr,"does not exist"
              write(STDERR, '(a)') " exiting..."
              stop 2
           end if
        end do
     end if
  end if

  means(1:ncomp) = ZERO
  k = nQTL * nChr
  call alloc2D(values_1D, k, nComp, "values_1D", "getQTLandSNP")
  QTLlist%nComp = nComp
  QTLList%nQTL = nQTL
  nReq = nQTL + nSNP
  call alloc1I(iReq, nReq, "iReq", "getQTLandSNP")
  call gnormal(means, covMat, nComp, k, values_1D)
  if (verbose) write(STDOUT, *) " random values for tbv created"
  maxLoci = 0
  do iChr = 1, nChr
     if (genome(iChr)%nLoci > maxLoci) maxLoci = genome(iChr)%nLoci
  end do
  call alloc1I(temp, maxLoci, "temp", "getQTLandSNP")
  temp = (/(i, i = 1, maxLoci)/)

  if (.not.randomMAF) then
     allocate(MAFArray(nChr))
     ICHRLOOP: do iChr = 1, nChr
        call alloc1D(MAFArray(iChr)%array, genome(iChr)%nloci, &
             "MAFArray(iChr)%array", "getQTLandSNP")
        write(fileName,'(a,i3.3)') trim(baseNameFreq), ichr
        open(newUnit = iun, file = fileName, status = 'old')
        read(iun, *)
        do i = 1, genome(iChr)%nloci
           read(iun, *) j, rand, k, freq
           MAFArray(iChr)%array(i) = min(freq, 1 - freq)
        end do
        close(iun)
        
        ! sortrx does not need to know the number of elements in temp
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
           write(STDERR, '(a)') "Erorr"
681        format(" Number of QTL is less than nLoci for Chromosome ", &
                i2," for MAF > ", f10.8)
682        format(" iChr, NReq, Nloci, available, < maf", i2, 4x, i4, 3x,&
                i4, 2x, i4, 2x, i4)
           write(STDERR, 681) iChr, maf
           write(STDERR, 682) iChr, NReq, genome(iChr)%nLoci, nAvail , i
           write(STDERR, '(a)') " exiting..."
           stop 2
        end if

        temp(1:nAvail) = (/(i, i = 1, nAvail)/)
        call choice(temp, maxloci, nAvail, nReq, iReq, nReq)

        QTLlist%indices(iChr, 1:nQTL) = iReq(1:nQTL)
        do i = 1, nComp
           k = (iChr - 1) * nQTL
           QTLlist%values(iChr, 1:nQTL, i) = values_1D((k + 1) : (k + nQTL), i)
        end do
        SNPlist(ichr, 1:nSNP) = iReq((nQTL+1):nReq)

     end do ICHRLOOP
     do ichr = 1, nchr
        deallocate(MAFArray(iChr)%array)
     end do
     deallocate(MAFArray)
     
  else ! i.e., when the selection is random

     ICHRLOOP2: do iChr = 1, nChr
        !          src,   dim,     size,              n,  output, dim
        call choice(temp, maxLoci, genome(iChr)%nLoci, nReq, iReq, nReq)
        QTLlist%indices(iChr, 1:nQTL) = iReq(1:nQTL)
        SNPlist(iChr, 1:nSNP) = iReq((nQTL+1):nReq)
        do i = 1, nComp
           j = (iChr - 1) * nQTL
           QTLlist%values(iChr, 1:nQTL, i) = values_1D((j+1):(j+nQTL), i)
        end do
     end do ICHRLOOP2
  end if

deallocate(iReq)
deallocate(temp)
deallocate(values_1D)

end subroutine getQTLandSNP

