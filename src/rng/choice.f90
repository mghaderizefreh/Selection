subroutine Choice(source, sourceSize, number, output, seed)
  !---------------------------------------------------------------------------------------------
  ! Samples "number" elements from "Source" with size "sourceSize" and stores in "output"

  integer , intent( in ) :: sourceSize, number
  integer, dimension(:), intent(in) :: source
  integer, dimension(:), intent(out) :: output
  integer, intent(in), dimension(:), optional :: seed

  integer, dimension(:), allocatable :: sourceCopy
  integer :: ind, temp, i
  real :: rand

  if (present(seed)) then
     call random_seed(put = seed)
  end if

  allocate(sourceCopy(sourceSize))

  ! First take a copy of source because arrays are passed by reference (are they?)
  forall (i = 1:sourceSize)
     sourceCopy(i) = source(i)
  end forall

  ! Then shuffle sourceCopy
  do i = 1, sourceSize
     call random_number(rand)
     do while(rand.eq.1) ! avoiding 1
        call random_number(rand)
     end do
     ind = int(rand * dble(sourceSize)) + 1 ! get a random index
     ! and swap index i with that random index
     temp = sourceCopy(i) 
     sourceCopy(i) = sourceCopy(ind)
     sourceCopy(ind) = temp
  end do

  ! finally get the first "number" and store in "output"
  forall (i= 1: number)
     output(i) = sourceCopy(i)
  end forall
end subroutine Choice

