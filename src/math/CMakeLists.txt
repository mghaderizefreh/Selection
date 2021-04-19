# math library

SET(MATH_src ${SRCMATH}/stats.f90
	     ${SRCMATH}/math.f90
)
ADD_LIBRARY(${MATHLIB} STATIC ${MATH_src})

TARGET_LINK_LIBRARIES(${MATHLIB} ${UTILLIB} ${CMAKE_THREAD_LIBS_INIT})
#				 ${BLAS_LIBRARIES}
#				  ${LAPACK_LIBRARIES}
                                 
