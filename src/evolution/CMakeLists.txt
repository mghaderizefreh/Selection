# reml library

SET(EVOL_src ${SRCEVOL}/initialiseGenotypes.f90
	     ${SRCEVOL}/getQTLandSNP.f90
	     ${SRCEVOL}/evolution_module.f90
)
ADD_LIBRARY(${EVOLLIB} STATIC ${EVOL_src})

#TARGET_LINK_LIBRARIES(${EVOLLIB} ${BLAS_LIBRARIES}
#                                 ${LAPACK_LIBRARIES}
#                                 ${CMAKE_THREAD_LIBS_INIT})
