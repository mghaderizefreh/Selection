# bar library

SET(BAR_src ${SRCBAR}/countNumberLines.f90
	    ${SRCBAR}/askFilename.f90
	    ${SRCBAR}/detInv.f90
	    ${SRCBAR}/dunpack.f90
	    ${SRCBAR}/trace.f90
	    ${SRCBAR}/global_module.f90
)
ADD_LIBRARY(${BARLIB} STATIC ${BAR_src})
