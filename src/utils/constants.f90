module constants
  use iso_fortran_env, only: KINDR => real64
  implicit none

  ! integer representation of alleles
  integer, parameter :: NBITS = 32

  !! derived types
  ! can make jagged arrays
  type JArr
     real(KINDR), dimension(:), allocatable :: array
  end type JArr

  type chromosome
     integer :: nloci, nblock
     real(KINDR):: chrL
     integer, dimension(:,:,:), pointer :: genotypes ! nanim x 2 x nblock
     ! if the loci are not equidist., then the position (no matter if Morgan or cM).
     real(KINDR), dimension(:), pointer :: positions ! nloci
  end type chromosome

  type variances
     real(KINDR), allocatable, dimension(:) :: A ! genetic part
     real(KINDR), allocatable, dimension(:) :: E ! environemntal part
     real(KINDR), allocatable, dimension(:) :: PE ! permanent environment
     real(KINDR), allocatable, dimension(:,:) :: corr ! (genetic) correlation
     real(KINDR), allocatable, dimension(:,:) :: cov ! (genetic) covariance
  end type variances

  type QTL_Array
     integer :: nQTL ! number of QTL (on one chromosome)
     integer :: nComp ! number of traits being affected
     integer, dimension(:,:), allocatable :: indices ! nchar x nQTL
     real(KINDR), dimension(:,:,:), allocatable :: values ! nChr x nQTL x nComp
  end type QTL_Array

  !   Handles
  integer, parameter :: STDIN  = 5
  integer, parameter :: STDOUT = 6
  integer, parameter :: STDERR = 0

  ! some other constants
  real(KINDR),    parameter :: ZERO   =  0_KINDR
  real(KINDR),    parameter :: ONE    =  1_KINDR
  real(KINDR),    parameter :: TWO    =  2_KINDR
  real(KINDR),    parameter :: FOURTH = .25_KINDR
  real(KINDR),    parameter :: HALF   =  .5_KINDR
  real(KINDR),    parameter :: THIRD  =  1._KINDR / 3._KINDR
  real(KINDR),    parameter :: PI     = 3.14159265358979323846_KINDR
  real(KINDR),    parameter :: SQRTPI = 1.7724538509055159_KINDR
  complex(KINDR), parameter :: I_C    = CMPLX(ZERO,    ONE,  KINDR)

  public :: alloc1L
  public :: alloc1I, alloc2I, alloc3Ip
  public :: alloc1D, alloc2D, alloc3D
  private :: HandleErr

contains
  subroutine alloc1L(array, len, name, from)
    implicit none
    logical, dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: len
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc1L
  subroutine alloc1I(array, len, name, from)
    implicit none
    integer, dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: len
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc1I
  subroutine alloc1D(array, len, name, from)
    implicit none
    real(KINDR), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: len
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc1D

  subroutine alloc2I(array, len1, len2, name, from)
    implicit none
    integer, dimension(:,:), allocatable, intent(inout) :: array
    integer, intent(in) :: len1, len2
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len1, len2), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc2I
  subroutine alloc2D(array, len1, len2, name, from)
    implicit none
    real(KINDR), dimension(:,:), allocatable, intent(inout) :: array
    integer, intent(in) :: len1, len2
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len1, len2), stat = i)
       if (i .ne. 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc2D

  subroutine alloc3Ip(array, len1, len2, len3, name, from)
    implicit none
    integer, dimension(:,:,:), pointer, intent(inout) :: array
    integer, intent(in) :: len1, len2, len3
    character(len=*), intent(in) :: name, from
    integer :: i
    allocate(array(len1, len2, len3), stat = i)
    if (i .ne. 0) call HandleErr(i, name, from)
  end subroutine alloc3Ip
  subroutine alloc3D(array, len1, len2, len3, name, from)
    implicit none
    real(KINDR), dimension(:,:,:), allocatable, intent(inout) :: array
    integer, intent(in) :: len1, len2, len3
    character(len=*), intent(in) :: name, from
    integer :: i
    if (allocated(array)) then
       i = -1
       call HandleErr(i, name, from)
    else
       allocate(array(len1, len2, len3), stat = i)
       if (i /= 0) call HandleErr(i, name, from)
    end if
  end subroutine alloc3D
  
  subroutine HandleErr(i, name, from)
    implicit none
    character(len=*), intent(in) :: name, from
    integer, intent(in) :: i
    if (i == -1) then
       write(STDERR, '(A)') "Error"
       write(STDERR,*) "array", trim(name), " in ", "is already allocated"
       stop 2
       if (i == 1) then
          write(STDERR,'(A)') "Error in system routine:"
          write(STDERR,*) "cannot allocate ", trim(name), " in ", trim(from)
          stop 2
       elseif (i == 2) then
          write(STDERR, '(A)') "Error (invalid data):"
          write(STDERR,*) "cannot allocate", trim(name), " in ", trim(from)
          stop 2
       elseif (i == 3) then
          write(STDERR, '(A)') "Error (invalid data and sys routine):"
          write(STDERR,*) "cannot allocate", trim(name), " in ", trim(from)
          stop 2
       else
          write(STDERR, '(A)') "Error:"
          write(STDERR, *) "cannot allocate", trim(name), " in ", trim(from),&
               "for unknown reason"
          stop 2
       end if
    end if
  end subroutine HandleErr
     
  
end module constants
   
