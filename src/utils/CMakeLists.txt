# utils library

SET(UTILS_src ${SRCUTILS}/constants.f90
	      ${SRCUTILS}/countNumberLines.f90
	      ${SRCUTILS}/askFilename.f90
	      ${SRCUTILS}/detInv.f90
	      ${SRCUTILS}/dunpack.f90
	      ${SRCUTILS}/trace.f90
	      ${SRCUTILS}/global_module.f90
)
ADD_LIBRARY(${UTILSLIB} STATIC ${UTILS_src})
