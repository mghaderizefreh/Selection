module reml_module
  use global_module
contains
  include "getZGZMats.f90"
  include "calculateV.f90"
  include "calculateP.f90"
  include "calculateAImat.f90"
  include "calculateLogL.f90"
  include "calculaterhs.f90"
  include "updatetheta.f90"
  include "getEffects.f90"
end module reml_module
