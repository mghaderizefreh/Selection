module blup_module
  use constants
  use global_module
contains
  include "getZGZMats.f90"
  include "calculateV.f90"
  include "calculateP.f90"
  include "getEffects.f90"
end module blup_module