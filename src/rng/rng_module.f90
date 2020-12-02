!     Created : MGZ  2 Dec 2020

!include 'mkl_vsl.f90'
module rng_module
  use constants
contains
  include "seeding.f90"
  include "choice.f90"
  include "gnormal.f90"
end module rng_module
!=======================================================================
