# rng library

SET(RNG_src ${SRCRNG}/seeding.f90
	     ${SRCRNG}/choice.f90
#	     ${SRCRNG}/inormal.f90
	     ${SRCRNG}/gnormal.f90
	     ${SRCRNG}/poissonProb.f90
	     ${SRCRNG}/rng_module.f90
)
ADD_LIBRARY(${RNGLIB} STATIC ${RNG_src})

TARGET_LINK_LIBRARIES(${RNGLIB} ${BLAS_LIBRARIES}
				 ${LAPACK_LIBRARIES}
                                 ${CMAKE_THREAD_LIBS_INIT})
