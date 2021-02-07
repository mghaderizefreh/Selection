!     Created : MGZ  2 Dec 2020

!include 'mkl_vsl.f90'
module rng_module

contains
  include "seeding.f90"
  include "choice.f90"
!  include "inormal.f90"
  include "gnormal.f90"
  include "poissonProb.f90"
end module rng_module
!=======================================================================
