module evolution_module
  private :: getPosition
contains
  include "initialiseGenotypes.f90"
  include "prepareMutRec.f90"
  include "getQTLandSNP.f90"
  include "simulateTBV.f90"
  include "simulatePhenotype.f90"
  include "covariate.f90"
  include "select.f90"
  include "SNP_map.f90"
  include "sampleGamete.f90"
  include "sampleMutation.f90"
end module evolution_module
