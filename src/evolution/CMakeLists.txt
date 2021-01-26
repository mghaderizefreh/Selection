# reml library

SET(EVOL_src ${SRCEVOL}/initialiseGenotypes.f90
	     ${SRCEVOL}/getQTLandSNP.f90
	     ${SRCEVOL}/covariate.f90
	     ${SRCEVOL}/simulateTBV.f90
	     ${SRCEVOL}/simulatePhenotype.f90
     	     ${SRCEVOL}/selection.f90
	     ${SRCEVOL}/evolution_module.f90
)
ADD_LIBRARY(${EVOLLIB} STATIC ${EVOL_src})

#TARGET_LINK_LIBRARIES(${EVOLLIB} ${BLAS_LIBRARIES}
#                                 ${LAPACK_LIBRARIES}
#                                 ${CMAKE_THREAD_LIBS_INIT})
