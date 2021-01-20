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
  end type variance

  type QTL_Array
     integer :: nQTL ! number of QTL (on one chromosome)
     integer :: nComp ! number of traits being affected
     integer, dimension(:,:), allocatable :: indices ! nchar x nQTL
     double precision, dimension(:,:,:), allocatable :: values ! nChr x nQTL x nComp
  end type QTL_Array

  !   Handles
  integer, parameter :: stdin  = 5
  integer, parameter :: stdout = 6
  integer, parameter :: stderr = 0

  !   Precision
  integer,        parameter :: KINDR = KIND(0d0)

  real(KINDR),    parameter :: ZERO        =  0.0e0_KINDR
  real(KINDR),    parameter :: ONE         =  1.0e0_KINDR
  real(KINDR),    parameter :: TWO         =  2.0e0_KINDR
  real(KINDR),    parameter :: THREE       =  3.0e0_KINDR
  real(KINDR),    parameter :: FOUR        =  4.0e0_KINDR
  real(KINDR),    parameter :: FIVE        =  5.0e0_KINDR
  real(KINDR),    parameter :: SIX         =  6.0e0_KINDR
  real(KINDR),    parameter :: SEVEN       =  7.0e0_KINDR
  real(KINDR),    parameter :: EIGHT       =  8.0e0_KINDR
  real(KINDR),    parameter :: NINE        =  9.0e0_KINDR
  real(KINDR),    parameter :: TEN         = 10.0e0_KINDR
  real(KINDR),    parameter :: ELEVEN      = 11.0e0_KINDR
  real(KINDR),    parameter :: TWELVE      = 12.0e0_KINDR
  real(KINDR),    parameter :: THIRTEEN    = 13.0e0_KINDR
  real(KINDR),    parameter :: FOURTEEN    = 14.0e0_KINDR
  real(KINDR),    parameter :: FIFTEEN     = 15.0e0_KINDR
  real(KINDR),    parameter :: SIXTEEN     = 16.0e0_KINDR
  real(KINDR),    parameter :: SEVENTEEN   = 17.0e0_KINDR
  real(KINDR),    parameter :: EIGHTEEN    = 18.0e0_KINDR
  real(KINDR),    parameter :: NINETEEN    = 19.0e0_KINDR
  real(KINDR),    parameter :: TWENTY      = 20.0e0_KINDR
  real(KINDR),    parameter :: TWENTYTWO   = 22.0e0_KINDR
  real(KINDR),    parameter :: TWENTYFOUR  = 24.0e0_KINDR
  real(KINDR),    parameter :: THIRTYTWO   = 32.0e0_KINDR
  real(KINDR),    parameter :: THIRTYSIX   = 36.0e0_KINDR
  real(KINDR),    parameter :: FORTYEIGHT  = 48.0e0_KINDR
  real(KINDR),    parameter :: SEVENTYTWO  = 72.0e0_KINDR
  real(KINDR),    parameter :: HUNDRED     =  1.0e2_KINDR
  real(KINDR),    parameter :: FOURTH      = 0.25e0_KINDR
  real(KINDR),    parameter :: HALF        =  0.5e0_KINDR
  real(KINDR),    parameter :: THREEFOURTH = 0.75e0_KINDR
  real(KINDR),    parameter :: THIRD       =  1.0e0_KINDR / 3.0e0_KINDR
  real(KINDR),    parameter :: TWOTHIRD    =  2.0e0_KINDR / 3.0e0_KINDR
  real(KINDR),    parameter :: SIXTH       =  1.0e0_KINDR / 6.0e0_KINDR
  real(KINDR),    parameter :: THREESECOND =  3.0e0_KINDR / 2.0e0_KINDR
  real(KINDR),    parameter :: FOURTHIRD   =  4.0e0_KINDR / 3.0e0_KINDR

  real(KINDR),    parameter :: PI            = 3.14159265358979323846_KINDR
  real(KINDR),    parameter :: SQRTPI        = 1.7724538509055159_KINDR
  real(KINDR),    parameter :: TWOINVPI      = TWO / PI
  real(KINDR),    parameter :: SQRTTWOINVPI  = 0.79788456080286541_KINDR
  real(KINDR),    parameter :: TWOINVSQRTPI  = TWO / SQRTPI

  complex(KINDR), parameter :: I_C       = CMPLX(ZERO,    ONE,  KINDR)
  complex(KINDR), parameter :: ZERO_C    = CMPLX(ZERO,    ZERO, KINDR)
  complex(KINDR), parameter :: ONE_C     = CMPLX(ONE,     ZERO, KINDR)
  complex(KINDR), parameter :: TWO_C     = CMPLX(TWO,     ZERO, KINDR)
  complex(KINDR), parameter :: THREE_C   = CMPLX(THREE,   ZERO, KINDR)
  complex(KINDR), parameter :: FOUR_C    = CMPLX(FOUR,    ZERO, KINDR)
  complex(KINDR), parameter :: FIVE_C    = CMPLX(FIVE,    ZERO, KINDR)
  complex(KINDR), parameter :: SIX_C     = CMPLX(SIX,     ZERO, KINDR)
  complex(KINDR), parameter :: HALF_C    = CMPLX(HALF,    ZERO, KINDR)
  complex(KINDR), parameter :: HUNDRED_C = CMPLX(HUNDRED, ZERO, KINDR)

end module constants
