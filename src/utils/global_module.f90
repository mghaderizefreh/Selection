module global_module
  use constants
  implicit none

contains
  include "askFilename.f90"
  include "countNumberLines.f90"
  include "detInv.f90"
  include "dunpack.f90"
  include "trace.f90"
  include "readData.f90"
end module global_module
