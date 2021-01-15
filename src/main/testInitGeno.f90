program testInitGeno

  use constants
  use global_module
  use rng_module
  use evolution_module

  implicit none
  integer :: genestart, nanim, nchr, istore
  integer :: nloci, nblock, maxloci, maxblock, ifail
  character(len=256) :: prefixfilename, startFile
  type(chromosome), DIMENSION(:), allocatable :: genome
  integer, dimension(:), allocatable :: seed

  startfile = "inicio.dat"
  call istart(seed, startfile, ifail)
  if (ifail /= 0) then
     write(stderr, '(a)') "reading/setting seed faild"
     stop 2
  end if

  genestart = 3
  nanim = 100
  nchr = 1
  istore = 1
  nloci = 300
  prefixfilename = "dumpedigree.ch"
  allocate(genome(nchr))
  call initialiseGenotypes(nchr, nanim, genestart, nloci, nblock, istore, genome,&
       maxloci, maxblock, ifail, prefixfilename)

end program testInitGeno
