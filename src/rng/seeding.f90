subroutine istart(seed, startfile, returnVal)
  !subroutine to read seed from startfile, if it doesnot exist, 
  ! be created and written on.
  use constants
  implicit none
  integer, intent(inout), dimension(:), allocatable :: seed
  character(len = 11), intent(in) :: startfile
  integer, intent(out) :: returnval
  logical :: iex
  integer , allocatable, dimension(:) :: initial
  integer :: iun, n
  character(len=7) :: fmto
  call random_seed(size = n)
  write(fmto, "(a1,i2,a4)") "(", n, "i12)"
  call alloc1I(initial, n, 'seed', 'istart')
  call alloc1I(seed, n, 'seed', 'istart')
  inquire (file=startfile,exist=iex)
  if (iex) then
     open(newUnit = iun, file=startfile, form='formatted',status='old')
     do
        read(iun, fmt = fmto, err=52, end=53) initial(1:N)
        seed(1:N) = initial(1:N)
     end do
53   continue
     write(STDOUT, '(a)')"file found"
     write(STDOUT, '(a)', advance = 'no')"initial value = "
     write(STDOUT, fmt = fmto) seed(1:N)
     close(iun)
     call random_seed(put = seed)
     returnval = 0
     return
  else
     call random_seed(get = seed)
     print*,"file for initializing random generator does not exist"
     print*,"it will be initialized  with ",seed
     open(newUnit = iun, file = startfile, form = 'formatted', status = 'new')
     initial = seed
     write(iun, fmto) seed
     close(iun)
     call random_seed(put = seed) ! redundant
     returnval = 0
     return
  endif
52 write(STDERR,*) " Error when reading initializing file ",startfile
  write(STDERR,*) " File may be for different purpose"
  returnval = 2
  deallocate(initial)
  deallocate(seed)
end subroutine istart

!=======================================================================

subroutine ifinal(seed,startfile)
  !	storing new initializing value in file
  use constants
  implicit none
  integer, dimension(:), intent(inout) :: seed
  character(len = 11), intent(in) :: startfile
  integer :: i, n, iun, HUGEN
  real(KINDR), dimension(:), allocatable :: values
  logical :: iex
  character(len = 7) :: fmto
  call random_seed(size = n)
  HUGEN = huge(n)
  write(fmto, "(a1,i2,a4)") "(", n, "i12)"
  call alloc1D(values, n, "values", "ifinal")
  call random_number(values)
  forall (i = 1:n)
     seed(i) = int( (2*values(i) - 1) * HUGEN)
  end forall
  write(STDOUT, '(a)', advance = 'no') "next initial value = "
  write(STDOUT, fmt = fmto) seed

  inquire (file=startfile,exist=iex)
  if (iex) then
     open(newUnit = iun, file=startfile,form='formatted',status='old',position = 'append')
  else
     open(newunit = iun,file=startfile,form='formatted',status='new')
  endif
  write(iun, fmt = fmto) seed
  close(10)
  deallocate(values)
end subroutine ifinal

