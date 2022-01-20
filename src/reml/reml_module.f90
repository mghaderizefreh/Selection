module reml_module
contains
  include "calculateAImat.f90"
  include "calculateLogL.f90"
  include "calculaterhs.f90"
  include "updatetheta.f90"
  include "reml.f90"
end module reml_module
