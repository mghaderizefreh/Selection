!    Creates QTLlist and SNPlist (for intercept and slope) based on genotype files
!    The list is made using SNPs with highest MAF or above a threshold. If latter is chosen,
!    frequency files should be provided and a threshold.
!    Inputs are:
!                number of chromosomes
!                base name for genotype of individuals for each chromosomes (it should be pedigree.ch???)
!                threshold value and frequency files in the form freq.txt??? if lists are to be chosen using a threshold
!                number of QTLs on each chromosomes (The sum of QTLs and SNPs should be less than total SNP for each chromosome)
!                number of SNPs on each chromosomes (The sum of QTLs and SNPs should be less than total SNP for each chromosome)
!                name of files for SNP and QTL list (normall SNPlist, QTLlist)
!    Outputs are
!                QTLlist.txt in the form `iChr, loc(i,j), val1, val2'
!                            where iChr is chromosome number and loc(i,j) is location of jth QTL on ith chromosome and val1 and 
!                            val2 are samples from two independent distribution with mean 0 and std 1
!                    SNPlist.txt if the form `iChr, loc(i,j)'
!                            where loc(i,j) in SNPlist MUST BE DIFFERENT than the one in QTLlist
!
! Written by Masoud Ghaderi. Modified and document on 4 Nov 2019
program makeListExt
  use constants
  use global_module
  use rng_module
  use quickSort
  implicit none

  character (len=20) :: startfile, baseNameGenotype, QTLfile, SNPfile, baseNameFreq
  integer :: nChr, iChr, iun, iun2, iun3, i, iostat, iid, nanim, randomSelection,k, j, ii
  integer :: nQTL, nSNP, nReq, nAvail, nComp
  double precision :: maf
  logical :: l_exists
  character(len = 30):: fileName, status, estatus
  character(len = 256):: corrStructFile
  integer, dimension(:), allocatable :: nLoci, nBlocks, iQTL,iSNP,iReq, temp, seed
  integer, dimension(:,:), allocatable :: QTLlist, SNPlist
  real, dimension(:,:,:), allocatable :: values
  real, dimension(:,:), allocatable  :: covMat, values_1D
  real, dimension(:), allocatable ::  means
  type JaggedArray
     real, dimension(:), allocatable :: ROW
  end type JaggedArray
  type(JaggedArray) , dimension(:) , allocatable:: MAFArray
  real :: rand, freq

  write(startfile,'(a)') "inicio.dat"
  call istart(seed, startfile, i)
  if (i /= 0) then
     write(STDERR, '(a)') "reading/setting seed faild"
     stop 2
  end if
  call random_seed(put = seed)

  call askInteger ( nChr, "number of chromosome genotype to be read (one file per chromosome)" )
  if(nChr <=0) then
     write(STDERR,'(a)') "Numbero of chromosomes was invalid"
     stop 2
  end if

  write (STDOUT, '(a)' ) " input genotype base file without numbers"
  read (STDIN,'(a20)') baseNameGenotype
  do ichr = 1, nchr
     write(filename, '(a,i3.3)') trim(baseNameGenotype), ichr
     inquire (file=filename, exist=l_exists)
     if (.not. l_exists) then
        write(STDERR,'(a4,1x,a,i3.3,1x,a)') "File", trim(baseNameGenotype),ichr,"does not exist"
        write(STDERR, '(a)') "exiting..."
        stop 2
     else
        write(STDOUT,'(a4,1x,a,i3.3,1x,a6)') "File", trim(baseNameGenotype),ichr, "exists."
        continue
     end if

  end do

  allocate(nloci(nChr))
  allocate(nblocks(nChr))
  call askInteger(nQTL, " Number of QTL")
  call askInteger(nSNP, " Number of SNP for SNP panel")
  nReq = nQTL + nSNP
  allocate(  iReq   (   nReq   ) )
  allocate( QTLlist (nChr, nQTL) )
  allocate( SNPlist (nChr, nSNP) )

  ! read files to extract nLoci
  do iChr = 1, nChr
     write(fileName, '(a,i3.3)') trim(baseNameGenotype), iChr
     open(newUnit = iun, file = fileName, status = 'old')
     read(iun, *, iostat = iostat) nanim, nLoci(iChr), nBlocks(iChr), iid
     if (nLoci(iChr) < nQTL) then
        write(STDERR,'(a)') " ERROR!"
        write(STDERR,'(a,i3)') " Number of QTL is less than nLoci for Chromosome", iChr
        write(STDERR,'(a,i3,i6)') " nQTL, nLoci", nQTL, nLoci(iChr)
        write(STDERR, '(a)') "exiting..."
        close(iun)
        stop 2
     end if
     close(iun)
  end do

  write(QTLfile, '(a20)') "QTLlist.txt"
  QTLfile = trim(QTLfile)
  write(SNPfile, '(a20)') "SNPlist.txt"
  SNPfile = trim(SNPfile)  

  randomSelection = 1
  call askInteger( randomSelection, " Selection of QTL is random with maf=0? (Y=1, N=0)")

  call askInteger(nComp, " Number of parameters")
  allocate( values  (nChr, nQTL,nComp), values_1D(nComp, nChr * nQTL), means(nComp), covMat(nComp,nComp))
  estatus = "o"
  call askFileName( corrStructFile, " Correlation structure file ", status, estatus)
  open(newUnit = iun3, file = corrStructFile)
  covMat(1:ncomp, 1:ncomp) = 0.0
  k = 1
  do i = 1, ncomp
     do j = 1, i
        read(iun3, *) covMat(i, j)
        covMat(j,i) = covMat(i,j)
     end do
  end do

  means(1:ncomp) = 0

  write(STDOUT, '(a)') "Correlation matrix was read"
  do i = 1, ncomp
     do j = 1, i
        write(6, '(f5.3,1x)',advance = 'no') covMat(i,j)
     end do
     write(6, *)
  end do
  call gnormal(means, covMat, ncomp, nQTL * nChr, values_1D)

  if (randomSelection .ne. 0) randomSelection = 1
  if (randomSelection .eq. 0) then

!     allocate(MAFArray(nChr))
!
!     write ( STDOUT, '(a)' ) ' input frequency base file without numbers'
!     read (STDIN,'(a20)') baseNameFreq
!     do ichr = 1, nchr
!        write(filename, '(a,i3.3)') trim(baseNameFreq), ichr
!        inquire (file=filename, exist=l_exists)
!        if (.not. l_exists) then
!           write(STDOUT,'(a4,1x,a,i3.3,1x,a)') "File", trim(baseNameFreq),ichr,"does not exist"
!           write(STDOUT, '(a)') "exiting..."
!           stop 2
!        else
!           write(6,'(a4,1x,a,i3.3,1x,a6)') "File", trim(baseNameFreq),ichr, "exists."
!           continue
!        end if
!     end do
!
!     write(STDOUT,'(a)') " Minimum frequency to filter for:"
!     read(STDIN,*) maf
!     if (maf > .5) then
!        write(STDERR, '(a)') " miniminm allele frequency cannot be more than 0.5"
!        write(STDERR, '(a)') " exiting..."
!        stop 2
!     end if

!     do iChr = 1, nChr
!        write(fileName,'(a,i3.3)') trim(baseNameFreq), ichr
!        open(newUnit = iun, file = fileName, status = 'old')
!        allocate(MAFArray(iChr)%ROW(nloci(iChr)))
!        do i = 1, nloci(iChr)
!           read(iun, *, iostat=iostat) iid, rand, iid, freq
!           MAFArray(iChr)%ROW(i) = min(freq, 1- freq)
!        end do
!        close(iun)
!
!        allocate(temp(nloci(iChr)))
!        call sortrx(nloci(ichr), MAFArray(ichr)%ROW,temp)
!
!        i = 0
!        do while(i < nLoci(iChr))
!           i = i + 1
!           if (MAFarray(ichr)%ROW(temp(i)) .gt. maf) then
!              i = i - 1
!              exit
!           end if
!        end do
!        deallocate(temp)
!
!        nAvail = nLoci(iChr) - i
!        if (nAvail < nReq) then
!681        format(" Number of QTL is less than nLoci for Chromosome ", i2," for MAF > ", f10.8)
!682        format(" iChr, NReq, Nloci, available, < maf", i2, 4x, i4, 3x, i4, 2x, i4, 2x, i4)
!           write(STDERR, 681) iChr, maf
!           write(STDERR, 682) iChr, NReq, nLoci(iChr), nAvail , i
!           write(STDERR, '(a)') "exiting..."
!           stop 2
!        end if
!
!        allocate(temp(nAvail))
!        do i = 1, nAvail
!           temp(i) = i
!        end do
!        call choice(temp, nAvail, nReq, iReq)
!        iQTL = iReq(1:nQTL)
!        iSNP = iReq(nQTL+1:nReq)
!        QTLlist(ichr,:) = iQTL
!        SNPlist(ichr,:) = iSNP
!        deallocate(temp)
!     end do
!
!     open(newUnit = iun , file = QTLfile, status = 'unknown')
!     open(newUnit = iun2, file = SNPfile, status = 'unknown')
!     do iChr = 1, nChr
!        do i = 1, nQTL
!           call normdev(rand)
!           values(ichr, i,1) = rand
!           call normdev(freq)
!           values(ichr, i,2) = freq
!           write(iun,'(i3,3x,i6,3x,f15.7,3x,f15.7)') iChr, QTLlist(ichr, i), rand, freq
!        end do
!        do i = 1, nSNP
!           WRITE(iun2, '(i3,3x,i6)') iChr, SNPlist(ichr, i)
!        end do
!     end do
!     close(iun )
!     close(iun2)

  else ! i.e., when the selection is random

     ! select nQTL loci randomly and assing values to them
     do iChr = 1, nChr
        allocate(temp(nLoci(iChr)))
        do i = 1, nLoci(iChr) ! making source for choice
           temp(i) = i
        end do
        call choice(temp, nLoci(iChr), nReq, iReq)
        iQTL = iReq(1:nQTL)
        iSNP = iReq(nQTL+1:nReq)
        QTLlist(iChr,:) = iQTL
        SNPlist(iChr,:) = iSNP
        deallocate(temp)
     end do

     open(newUnit = iun , file = QTLfile, status = 'unknown')
     open(newUnit = iun2, file = SNPfile, status = "unknown")
     do iChr = 1, nChr
        do i = 1, nQTL
           k = (iChr - 1) * nQTL + i 
           write(iun, '(i3,3x,i6,3x)', advance = 'no') iChr, QTLlist(iChr,i)
           do ii = 1, ncomp
              write(iun, '(f15.7)' , advance = 'no') values_1D(ii,k)
           end do
           write(iun, *)
        end do
        do i = 1, nSNP
           WRITE(iun2, '(i3,3x,i6)') iChr, SNPlist(ichr, i)
        end do
     end do
     close(iun )
     close(iun2)
  end if

  call ifinal(seed, startfile)

end program makeListExt
