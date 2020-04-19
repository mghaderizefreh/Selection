module REML_module
  implicit none

  type doublePre_Array
     double precision, dimension(:), pointer :: level(:)
  end type doublePre_Array

contains
  include "getZGZMats.f90"
  include "dsptrf_Ldet.f90"
  include "countNumberLines.f90"
  include "askFilename.f90"
  include "RRiteration.f90"
  include "getEffects.f90"

end module REML_module

