!
! Ricardo Pong-Wong 28/08/2010
!
!
!----------------------------------------------------------------------------
! this subroutine counts the number of of lines in a <ascii> file
! after skiping the first <skip> lines. This number <skip> is passed to 
! the subroutine as a variable
!
! there is the option to count all lines and all non empty lines.
!
! This information may be needed to allocate the variable to the right size
! 
!
!--------------------------------------------------------------------------
!
! the way it is done is by reading the file.
!
! HENCE IT CAN BE VERY UNEFFICIENT 
! (but this is the only way to optain this information with fortran)
!
! AVOID A MUCH AS POSSIBLE TO USE THE SUBROUTINE 
!
! (it is more efficient if the user provides such information to the programme 
! as input)
!
! DUMMY VARIABLE 
!     filename (input) : name of file to be counted
!     skip     (input) :  numbre of lines to be skip before start counting
!     lines    (output): number of lines ounted in file
!     empties  (input/output) : as input. if > 0 then count the number of empty lines
!                                         if <=0 do not count number of empty lines
!                               as output: return number of empties lines if variable at input was set to be  > 0
!----------------------------------------------------------------------------

subroutine countNumberLines(filename, skip, lines, empties, ifail)
  use constants
  implicit none

  character(len=*) :: filename
  integer          :: skip
  integer          :: lines
  integer          :: empties 
  integer          :: ifail

  character(len=1) :: aline

  integer :: iun,i, nonempties

  open(newunit=iun, file=filename, err=100, status='old')

  ifail=0
  do i=1,skip
     read(iun,*,end=10)
  enddo

  lines=0
  do 
     read(iun,*,end=10)
     lines=lines+1
  enddo
10 continue
  !write(stdout,*)'lines= ',lines

  if(lines==0) then
     empties=0
     close(iun)
     return
  endif

  !===========================================
  if(empties > 0) then
     rewind(iun)
     do i=1,skip
        read(iun,*,end=11)
     enddo

     nonempties=0
     do 
        read(iun,*,end=11)aline
        nonempties=nonempties+1
     enddo
11   continue
     empties=lines-nonempties
     !write(stdout,*)'nonempties = ',nonempties

  endif
  close (iun)

  return


100 continue
  write ( stderr, * ) ' Error opening file'
  write ( stderr, * ) ' lines cannot be counted'

  ifail  = 1
  lines  =-1
  empties= 0
  return

end subroutine countNumberLines


