# blup library

SET(BLUP_src ${SRCBLUP}/getZGZMats.f90
	     ${SRCBLUP}/calculateV.f90
	     ${SRCBLUP}/calculateP.f90
	     ${SRCBLUP}/getEffects.f90
    	     ${SRCBLUP}/BSRibsCalc.f90
    	     ${SRCBLUP}/leastSquare.f90
    	     ${SRCBLUP}/relationMatrix.f90
	     ${SRCBLUP}/getVariance.f90
	     ${SRCBLUP}/blup_module.f90
)
ADD_LIBRARY(${BLUPLIB} STATIC ${BLUP_src})

TARGET_LINK_LIBRARIES(${BLUPLIB} ${MATHLIB} ${SORTLIB}
				 ${BLAS_LIBRARIES}
                                 ${LAPACK_LIBRARIES}
                                 ${CMAKE_THREAD_LIBS_INIT})
