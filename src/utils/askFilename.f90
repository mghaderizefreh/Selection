subroutine askFilename(filename, question, theStatus, expectStatus)
  use constants, only : STDOUT, STDERR, STDIN
  implicit none

  character ( len = * ) :: filename  !need to be sure that it will not go out side
  character ( len = * ) :: question
  character ( len = * ) :: theStatus
  character ( len = * ) :: expectStatus
  logical :: L_EXISTS,doloop

  integer :: icount,maxtime

  icount=0
  maxtime=10  ! number of times to input a valid filename given expectStatus. if greter the 

  doloop=.true.
  do while ( doloop )
     theStatus = ''
     if(len_trim(question) < 1) then
        write(STDOUT, * ) ' filename'
     else
        write (STDOUT, * ) trim(question)
     endif

     read  (STDIN, * ) filename
     write (STDOUT, * ) trim(filename)
     inquire (file=filename, exist=l_exists)

     if ( l_exists ) then
        theStatus = 'o'
        if(expectStatus(1:1) .ne. 's' .and. expectStatus(1:1) .ne. 'o') &
            write(STDOUT,*)' file exists'
        if(expectStatus(1:1) .ne. 'n') doloop=.false.
     else
        theStatus = 'n'
        if(expectStatus(1:1) .ne. 's' .and. expectStatus(1:1) .ne. 'n') &
            write(STDERR,*)' file does NOT exists'
        if(expectStatus(1:1) .ne. 'o') doloop=.false.
     endif
     if(doloop) then
        icount=icount+1
        if(icount > maxtime) then  ! it will allow only maxtime mistake. If continue then it will go out with empty value
           write (STDERR, * )
           write (STDERR, * ) '============================================'
           write (STDERR, * )
           write (STDERR, * ) ' no valid file given after several attempts', icount
           write (STDERR, * ) ' check if input variable is too long '
           write (STDERR, * ) ' if it is the case then the read variable will be truncated'
           write (STDERR, * ) ' in this is the case the file need to be renamed with a shorter name'
           write (STDERR, * ) ' also filename MUST not have a space characters (as they will be triuncated)'  
           write (STDERR, * )
           write (STDERR, * ) ' last input was |', trim( filename ), '|'
           write (STDERR, * )
           write (STDERR, * ) '============================================'
           write (STDERR, * )
           filename = ''
           theStatus = 'x'
           doloop=.false.
        end if

     end if
  end do
  write(STDOUT,*)

end subroutine askFilename

!===============================================================================
subroutine askInteger ( intVal, question )
  use constants, only : STDOUT, STDIN
  implicit none

  integer :: intVal
  character ( len = * ) :: question
  integer :: ierr
  integer :: maxtime,ntimes

  maxtime=10
  ntimes=0

  do
     if(len_trim(question) < 1) then
        write ( STDOUT, * ) ' input integer'
     else
        write ( STDOUT, * ) trim(question)
     endif

     read  ( STDIN, *, iostat = ierr  ) intVal
     if(ierr == 0) exit
     write(STDOUT,*)' enter a valid integer', ierr
     ntimes=ntimes+1
     if(ntimes > maxtime) then
        write(STDOUT,*)' wrong input several times', ntimes
        write(STDOUT,*)' it will be given as default a -1'
        intVal=-1
        exit
     endif
  enddo

  write ( STDOUT, * ) intVal

end subroutine askInteger

!===============================================================================
subroutine askYesNoInteger ( intVal, question, default )
  use constants, only : STDOUT, STDIN
  implicit none

  integer :: intVal
  character ( len = * ) :: question
  integer :: default
  integer :: ierr,maxtime, ntimes

  maxtime=10
  ntimes=0

  do 
     if(len_trim(question) < 1) then
        write ( STDOUT, * ) ' yes or not? (Y=1, N=0)'
     else
        write ( STDOUT, * ) trim(question), ' (Y=1, N=0)'
     endif

     read  ( STDIN, *, iostat = ierr  ) intVal
     if(ierr /= 0) then
        ntimes=ntimes+1
        if( ntimes> maxtime) then
           write(STDOUT,*)' wrong input several times', ntimes
           write(STDOUT,*)' need it a 0/1 iput and non integer was given'
           write(STDOUT,*)' it will be set the default'
           intVal=default
        else
           write(STDOUT,*)' enter a valid integer'
           cycle
        endif
     endif


     if(default== 0)then
        if(intVal/=1) intVal = 0
        exit
     elseif(default==1) then
        if ( intVal /= 0 ) intVal = 1
        exit
     else
        if(intVal .ne.0 .and. intVal .ne. 1)then
           write(STDOUT,*)' input must be 0 (no) or 1 (yes)'
           write(STDOUT,*)
           ntimes=ntimes+1
           if(ntimes> maxtime) then
              write(STDOUT,*)' set to be the default'
              intVal=default
              exit
           endif
        else
           exit
        endif
     endif
  enddo

  write ( STDOUT, * ) intVal

end subroutine askYesNoInteger
