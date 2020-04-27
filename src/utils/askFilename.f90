subroutine askFilename(filename, question, theStatus, expectStatus)
  use constants
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
        write(stdout, * ) ' filename'
     else
        write (stdout, * ) trim(question)
     endif

     read  (stdin, * ) filename
     write (stdout, * ) trim(filename)
     INQUIRE (file=filename, EXIST=L_EXISTS)

     if ( L_EXISTS ) then
        theStatus = 'o'
        if(expectStatus(1:1) .ne. 's' .and. expectStatus(1:1) .ne. 'o') write(stdout,*)' file exists'
        if(expectStatus(1:1) .ne. 'n') doloop=.false.
     else
        theStatus = 'n'
        if(expectStatus(1:1) .ne. 's' .and. expectStatus(1:1) .ne. 'n') write(stdout,*)' file does NOT exists'
        if(expectStatus(1:1) .ne. 'o') doloop=.false.
     endif
     if(doloop) then
        icount=icount+1
        if(icount > maxtime) then  ! it will allow only maxtime mistake. If continue then it will go out with empty value
           write ( stdout, * )
           write ( stdout, * ) '============================================'
           write ( stdout, * )
           write ( stdout, * ) ' no valid file given after several attempts', icount
           write ( stdout, * ) ' check if input variable is too long '
           write ( stdout, * ) ' if it is the case then the read variable will be truncated'
           write ( stdout, * ) ' in this is the case the file need to be renamed with a shorter name'
           write ( stdout, * ) ' also filename MUST not have a space characters (as they will be triuncated)'  
           write ( stdout, * )
           write ( stdout, * ) ' last input was |', trim( filename ), '|'
           write ( stdout, * )
           write ( stdout, * ) '============================================'
           write ( stdout, * )
           filename = ''
           theStatus = 'x'
           doloop=.false.
        endif

     endif
  enddo
  write(stdout,*)

  return
end subroutine askFilename



!===============================================================================
subroutine askInteger ( intVal, question )
  use constants
  implicit none

  integer :: intVal
  character ( len = * ) :: question
  integer :: ierr
  integer :: maxtime,ntimes

  maxtime=10
  ntimes=0

  do
     if(len_trim(question) < 1) then
        write ( stdout, * ) ' input integer'
     else
        write ( stdout, * ) trim(question)
     endif

     read  ( stdin, *, iostat = ierr  ) intVal
     if(ierr == 0) exit
     write(stdout,*)' enter a valid integer', ierr
     ntimes=ntimes+1
     if(ntimes > maxtime) then
        write(stdout,*)' wrong input several times', ntimes
        write(stdout,*)' it will be given as default a -1'
        intVal=-1
        exit
     endif
  enddo

  write ( stdout, * ) intVal


  return
end subroutine askInteger

!===============================================================================
subroutine askYesNoInteger ( intVal, question, default )
  use constants
  implicit none

  integer :: intVal
  character ( len = * ) :: question
  integer :: default
  integer :: ierr,maxtime, ntimes

  maxtime=10
  ntimes=0

  do 
     if(len_trim(question) < 1) then
        write ( stdout, * ) ' yes or not? (Y=1, N=0)'
     else
        write ( stdout, * ) trim(question), ' (Y=1, N=0)'
     endif

     read  ( stdin, *, iostat = ierr  ) intVal
     if(ierr /= 0) then
        ntimes=ntimes+1
        if( ntimes> maxtime) then
           write(stdout,*)' wrong input several times', ntimes
           write(stdout,*)' need it a 0/1 iput and non integer was given'
           write(stdout,*)' it will be set the default'
           intVal=default
        else
           write(stdout,*)' enter a valid integer'
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
           write(stdout,*)' input must be 0 (no) or 1 (yes)'
           write(stdout,*)
           ntimes=ntimes+1
           if(ntimes> maxtime) then
              write(stdout,*)' st to be the default'
              intVal=default
              exit
           endif
        else
           exit
        endif
     endif
  enddo

  write ( stdout, * ) intVal

  return
end subroutine askYesNoInteger
