module constants
  implicit none

  type doublePre_Array
     double precision, dimension(:), pointer :: level(:)
  end type doublePre_Array

  type JArrD
     double precision, dimension(:), allocatable :: array
  end type JArrD

  type JArrR
     real, dimension(:), allocatable :: array
  end type JArrR

  type chromosome
     integer :: nloci, nblock
     double precision  :: chrL
     integer, dimension(:,:,:), pointer :: genotypes ! nanim x 2 x nblock
     ! if the loci are not equidist., then the position (no matter if Morgan or cM).
     real , dimension(:), pointer :: positions ! nloci
  end type chromosome

  type variance
     double precision, allocatable, dimension(:) :: A ! genetic part
     double precision, allocatable, dimension(:) :: E ! environemntal part
     double precision, allocatable, dimension(:) :: PE ! permanent environment
     real, allocatable, dimension(:,:) :: corr ! (genetic) correlation
  end type variance

  type QTL_Array
     integer :: nQTL ! number of QTL (on one chromosome)
     integer :: nComp ! number of traits being affected
     integer, dimension(:,:), allocatable :: indices ! nchar x nQTL
     double precision, dimension(:,:,:), allocatable :: values ! nChr x nQTL x nComp
  end type QTL_Array

  !   Handles
  integer, parameter :: STDIN  = 5
  integer, parameter :: STDOUT = 6
  integer, parameter :: STDERR = 0

  !   Precision
  integer,        parameter :: KINDR = KIND(0d0)

  real(KINDR),    parameter :: ZERO   =  0.0e0_KINDR
  real(KINDR),    parameter :: ONE    =  1.0e0_KINDR
  real(KINDR),    parameter :: TWO    =  2.0e0_KINDR
  real(KINDR),    parameter :: FOURTH = 0.25e0_KINDR
  real(KINDR),    parameter :: HALF   =  0.5e0_KINDR
  real(KINDR),    parameter :: THIRD  =  1.0e0_KINDR / 3.0e0_KINDR
  real(KINDR),    parameter :: PI     = 3.14159265358979323846_KINDR
  real(KINDR),    parameter :: SQRTPI = 1.7724538509055159_KINDR
  complex(KINDR), parameter :: I_C    = CMPLX(ZERO,    ONE,  KINDR)

end module constants
