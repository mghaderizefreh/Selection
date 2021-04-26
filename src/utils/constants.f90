module constants
  implicit none

  !   Precision
  integer,        parameter :: KINDR = KIND(0d0)

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

end module constants
