subroutine askFilename(filename, question, theStatus, expectStatus)
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
            write(6, * ) ' filename'
        else
            write (6, * ) trim(question)
        endif

        read  (5, * ) filename
        write (6, * ) trim(filename)
        INQUIRE (file=filename, EXIST=L_EXISTS)

        if ( L_EXISTS ) then
            theStatus = 'o'
            if(expectStatus(1:1) .ne. 's' .and. expectStatus(1:1) .ne. 'o') write(6,*)' file exists'
            if(expectStatus(1:1) .ne. 'n') doloop=.false.
        else
            theStatus = 'n'
            if(expectStatus(1:1) .ne. 's' .and. expectStatus(1:1) .ne. 'n') write(6,*)' file does NOT exists'
            if(expectStatus(1:1) .ne. 'o') doloop=.false.
        endif
        if(doloop) then
            icount=icount+1
            if(icount > maxtime) then  ! it will allow only maxtime mistake. If continue then it will go out with empty value
                write ( 6, * )
                write ( 6, * ) '============================================'
                write ( 6, * )
                write ( 6, * ) ' no valid file given after several attempts', icount
                write ( 6, * ) ' check if input variable is too long '
                write ( 6, * ) ' if it is the case then the read variable will be truncated'
                write ( 6, * ) ' in this is the case the file need to be renamed with a shorter name'
                write ( 6, * ) ' also filename MUST not have a space characters (as they will be triuncated)'  
                write ( 6, * )
                write ( 6, * ) ' last input was |', trim( filename ), '|'
                write ( 6, * )
                write ( 6, * ) '============================================'
                write ( 6, * )
                filename = ''
                theStatus = 'x'
                doloop=.false.
            endif

        endif
    enddo
    write(6,*)
 
    return
end subroutine



!===============================================================================
subroutine askInteger ( intVal, question )
    implicit none

    integer :: intVal
    character ( len = * ) :: question
    integer :: ierr
    integer :: maxtime,ntimes
    
    maxtime=10
    ntimes=0
    
    do
        if(len_trim(question) < 1) then
            write ( 6, * ) ' input integer'
        else
            write ( 6, * ) trim(question)
        endif

        read  ( 5, *, iostat = ierr  ) intVal
        if(ierr == 0) exit
        write(6,*)' enter a valid integer', ierr
        ntimes=ntimes+1
        if(ntimes > maxtime) then
           write(6,*)' wrong input several times', ntimes
           write(6,*)' it will be given as default a -1'
           intVal=-1
           exit
        endif
     enddo

        write ( 6, * ) intVal


    return
end subroutine

!===============================================================================
subroutine askYesNoInteger ( intVal, question, default )
    implicit none

    integer :: intVal
    character ( len = * ) :: question
    integer :: default
    integer :: ierr,maxtime, ntimes
    
    maxtime=10
    ntimes=0

    do 
        if(len_trim(question) < 1) then
            write ( 6, * ) ' yes or not? (Y=1, N=0)'
        else
            write ( 6, * ) trim(question), ' (Y=1, N=0)'
        endif

        read  ( 5, *, iostat = ierr  ) intVal
        if(ierr /= 0) then
           ntimes=ntimes+1
           if( ntimes> maxtime) then
              write(6,*)' wrong input several times', ntimes
              write(6,*)' need it a 0/1 iput and non integer was given'
              write(6,*)' it will be set the default'
              intVal=default
           else
              write(6,*)' enter a valid integer'
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
                write(6,*)' input must be 0 (no) or 1 (yes)'
                write(6,*)
                ntimes=ntimes+1
                if(ntimes> maxtime) then
                   write(6,*)' st to be the default'
                   intVal=default
                   exit
                endif
             else
                exit
             endif
        endif
    enddo

    write ( 6, * ) intVal

    return
end subroutine
