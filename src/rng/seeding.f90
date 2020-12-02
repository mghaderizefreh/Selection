subroutine istart(seed, startfile, returnVal)
  !subroutine to read seed from startfile, if it doesnot exist, 
  ! be created and written on.
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
  allocate(initial(n), seed(n))
  inquire (file=startfile,exist=iex)
  if (iex) then
     open(newUnit = iun, file=startfile, form='formatted',status='old')
     do
        read(iun, fmt = fmto, err=52, end=53) initial(1:N)
        seed(1:N) = initial(1:N)
     end do
53   continue
     write(stdout, '(a)')"file found"
     write(stdout, '(a)', advance = 'no')"initial value = "
     write(stdout, fmt = fmto) seed(1:N)
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
52 write(stderr,*) " Error when reading initializing file ",startfile
  write(stderr,*) " File may be for different purpose"
  returnval = 2
  return
end subroutine istart

!=======================================================================

subroutine ifinal(seed,startfile)
  !	storing new initializing value in file
  implicit none
  integer, dimension(:), intent(inout) :: seed
  character(len = 11), intent(in) :: startfile
  integer :: i, n, iun, HUGEN
  real, dimension(:), allocatable :: values
  logical :: iex
  character(len = 7) :: fmto
  call random_seed(size = n)
  HUGEN = huge(n)
  write(fmto, "(a1,i2,a4)") "(", n, "i12)"
  allocate(values(n))
  call random_number(values)
  forall (i = 1:n)
     seed(i) = int( (2*values(i) - 1) * HUGEN)
  end forall
  write(stdout, '(a)', advance = 'no') "next initial value = "
  write(stdout, fmt = fmto) seed

  inquire (file=startfile,exist=iex)
  if (iex) then
     open(newUnit = iun, file=startfile,form='formatted',status='old',position = 'append')
  else
     open(newunit = iun,file=startfile,form='formatted',status='new')
  endif
  write(iun, fmt = fmto) seed
  close(10)

end subroutine ifinal

