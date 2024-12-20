subroutine Choice(source, sourceDim, sourceSize, number, output, outputDim, seed)
  !---------------------------------------------------------------------------------------------
  ! Samples "number" elements from "Source" with size "sourceSize" and stores in "output"
  ! sourceDim is the number of elements in source
  ! sourceSize is the first number of elements in source that will be used
  ! 'number' elements will be returned in ouput
  ! output has outputDim number of elements
  ! sourceDim >= sourceSize >= number and outputDim >= number

  use constants, only : KINDR, alloc1I
  implicit none

  integer, intent(in) :: sourceSize, number, sourceDim, outputDim
  integer, dimension(1:sourceDim), intent(in) :: source
  integer, dimension(1:outputDim), intent(out) :: output
  integer, intent(in), dimension(:), optional :: seed

  integer, dimension(:), allocatable :: sourceCopy
  integer :: ind, temp, i
  real(KINDR) :: rand

  if (present(seed)) then
     call random_seed(put = seed)
  end if
  call alloc1I(sourceCopy, sourceSize, "sourceCopy", "choice")
  ! First take a copy of source because source array is not changeable
  sourceCopy(1:sourceSize)= source(1:sourceSize)

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
  output(1: number) = sourceCopy(1: number)
end subroutine Choice

