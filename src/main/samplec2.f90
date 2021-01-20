program sample2
  use constants
  use rng_module
  use evolution_module
  
  implicit none
  
  integer, dimension(:), allocatable :: seed
  character(len = 20) :: startfile = "inicio.dat", prefixfilename, formato
  integer :: i, ifail, j
  integer :: nloci, nblock, maxloci, maxblock
  type(chromosome), DIMENSION(:), allocatable :: genome
  real, dimension(:,:), allocatable :: covMat
  integer, dimension(:,:), allocatable :: SNPlist
  type(QTL_Array) :: QTLlist

  double precision :: val1
  integer :: nchr, nqtl, nsnp, ncomp
  logical :: rmaf
  call istart(seed, startfile, i)
  if (i /= 0) then
     write(STDERR, '(a)') "reading/setting seed faild"
     stop 2
  end if
  call random_seed(put = seed)

  prefixfilename = "dumpedigree.ch"
  allocate(genome(1))
  call initialiseGenotypes(1, 900, 3, nloci, nblock, 1, genome,&
       maxloci, maxblock, ifail, prefixfilename)
  
 
  open(1, file = 'corr')
  allocate(covmat(5,5))
  covMat(1:5, 1:5) = 0.0
  do i = 1, 5
     do j = 1, i
        read(1, *) covMat(i, j)
        covMat(j,i) = covMat(i,j)
     end do
  end do
  close(1)

  nchr = 1
  nQTL = 500
  nSNP = 1000
  ncomp = 5
  rmaf = .false.
  val1 = 0.1d0
  prefixfilename = "freq.txt"
  call getQTLandSNP(nchr, nQTL, nSNP, nComp, rmaf, genome, QTLlist, &
       SNPlist, covMat, prefixfilename, val1)

  open(1,file = 'QTLlist.txt')
  open(2,file = 'SNPlist.txt')
  write(formato, '(a1,i2,a12)') "(", nComp, "(f15.7, x))"
  do i = 1, nQTL
     write(1, '(i3,3x,i6)', advance = 'no') 1, QTLlist%indices(1, i)
     write(1, formato) QTLlist%values(1, i, 1:nComp)
  end do
  do i = 1, nSNP
     write(2, '(i3,3x,i6)') 1, SNPlist(1, i)
  end do

  call ifinal(seed, startfile)
end program sample2

