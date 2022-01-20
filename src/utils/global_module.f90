module global_module
  implicit none
  private
  public :: askFilename, detInv, traceA, traceAxB, traceAxBdiag
  public :: readInput, askInteger, askYesNoInteger
  public :: countNumberLines, dunpack  
contains
  include "askFilename.f90"
  include "countNumberLines.f90"
  include "detInv.f90"
  include "dunpack.f90"
  include "trace.f90"
  include "readData.f90"
end module global_module
