# reml library

SET(REML_src ${SRCREML}/calculateAImat.f90
	     ${SRCREML}/calculateLogL.f90
	     ${SRCREML}/calculaterhs.f90
	     ${SRCREML}/updatetheta.f90
	     ${SRCREML}/reml_module.f90
)
ADD_LIBRARY(${REMLLIB} STATIC ${REML_src})

TARGET_LINK_LIBRARIES(${REMLLIB} ${BLAS_LIBRARIES}
                                 ${LAPACK_LIBRARIES}
                                 ${CMAKE_THREAD_LIBS_INIT})
