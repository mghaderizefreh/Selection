module evolution_module
  use constants
!  use global_module
contains
  include "initialiseGenotypes.f90"
  include "getQTLandSNP.f90"
end module evolution_module
