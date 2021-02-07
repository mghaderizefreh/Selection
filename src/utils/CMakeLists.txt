# utils library

SET(UTIL_src ${SRCUTIL}/constants.f90
	     ${SRCUTIL}/countNumberLines.f90
	     ${SRCUTIL}/askFilename.f90
	     ${SRCUTIL}/detInv.f90
	     ${SRCUTIL}/dunpack.f90
	     ${SRCUTIL}/trace.f90
	     ${SRCUTIL}/readData.f90
	     ${SRCUTIL}/global_module.f90
)
ADD_LIBRARY(${UTILLIB} STATIC ${UTIL_src})

TARGET_LINK_LIBRARIES(${UTILLIB} ${BLAS_LIBRARIES}
				 ${LAPACK_LIBRARIES}
                                 ${CMAKE_THREAD_LIBS_INIT})
