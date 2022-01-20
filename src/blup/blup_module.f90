module blup_module

contains
  include "getZGZMats.f90"
  include "calculateV.f90"
  include "calculateP.f90"
  include "getEffects.f90"
  include "blup.f90"
  include "BSRibsCalc.f90"
  include "leastSquare.f90"
  include "relationMatrix.f90"
  include "getVariance.f90"
end module blup_module
