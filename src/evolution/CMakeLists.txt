# evolution library

SET(EVOL_src ${SRCEVOL}/initialiseGenotypes.f90
	     ${SRCEVOL}/prepareMutRec.f90
	     ${SRCEVOL}/getQTLandSNP.f90
	     ${SRCEVOL}/covariate.f90
	     ${SRCEVOL}/simulateTBV.f90
	     ${SRCEVOL}/simulatePhenotype.f90
	     ${SRCEVOL}/select.f90
	     ${SRCEVOL}/SNP_map.f90
	     ${SRCEVOL}/sampleGamete.f90
	     ${SRCEVOL}/sampleMutation.f90
	     ${SRCEVOL}/evolution_module.f90
)
ADD_LIBRARY(${EVOLLIB} STATIC ${EVOL_src})

TARGET_LINK_LIBRARIES(${EVOLLIB} ${RNGLIB} ${SORTLIB} ${CMAKE_THREAD_LIBS_INIT})
#${BLAS_LIBRARIES}
#                                 ${LAPACK_LIBRARIES}

