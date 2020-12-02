########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(RR_src ${MAIN}/RRiteration.f90
	   ${MAIN}/RRReml.f90
)

SET(ST_src ${MAIN}/STReml.f90
)

SET(BL_src ${MAIN}/blup.f90
)

SET(SA_src ${MAIN}/sample.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${RR} ${RR_src})
ADD_EXECUTABLE(${ST} ${ST_src})
ADD_EXECUTABLE(${BL} ${BL_src})
ADD_EXECUTABLE(${SA} ${SA_src})
#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links reml to the utils library
TARGET_LINK_LIBRARIES(${RR} ${UTILLIB} ${BLUPLIB} ${REMLLIB})
TARGET_LINK_LIBRARIES(${ST} ${UTILLIB} ${SORTLIB} ${BLUPLIB})
TARGET_LINK_LIBRARIES(${BL} ${UTILLIB} ${BLUPLIB})
TARGET_LINK_LIBRARIES(${SA} ${UTILLIB} ${RNGLIB} ${SORTLIB})

# Uncomment if you need to link to BLAS and LAPACK
TARGET_LINK_LIBRARIES(${RR} ${BLAS_LIBRARIES}
			    ${LAPACK_LIBRARIES}
                            ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(${ST} ${BLAS_LIBRARIES}
                            ${LAPACK_LIBRARIES}
                            ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(${SA} ${BLAS_LIBRARIES}
                            ${LAPACK_LIBRARIES}
                            ${CMAKE_THREAD_LIBS_INIT})

# Uncomment if you have parallization
#IF(USE_OPENMP)
#    SET_TARGET_PROPERTIES(${RR} PROPERTIES
#                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
#                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
#ELSEIF(USE_MPI)
#    SET_TARGET_PROPERTIES(${RR} PROPERTIES
#                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
#                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
#    TARGET_LINK_LIBRARIES(${RR} ${MPI_Fortran_LIBRARIES})
#ENDIF(USE_OPENMP)

#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${RR} RUNTIME DESTINATION bin)
