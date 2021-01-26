########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(analysis_src ${MAIN}/analysis.f90)

SET(FE_src ${MAIN}/fixedEff.f90
)

SET(SA_src ${MAIN}/sample.f90
)
SET(SA2_src ${MAIN}/samplec2.f90
)

SET(SEL_src ${MAIN}/main.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${analysis} ${analysis_src})
ADD_EXECUTABLE(${SA} ${SA_src})
ADD_EXECUTABLE(${SA2} ${SA2_src})
ADD_EXECUTABLE(${FE} ${FE_src})
ADD_EXECUTABLE(${SEL} ${SEL_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links exes to the their libraries
TARGET_LINK_LIBRARIES(${analysis} ${UTILLIB} ${BLUPLIB} ${REMLLIB})
TARGET_LINK_LIBRARIES(${SA} ${UTILLIB} ${SORTLIB} ${RNGLIB})
TARGET_LINK_LIBRARIES(${SA2} ${UTILLIB} ${SORTLIB} ${RNGLIB} ${EVOLLIB})
TARGET_LINK_LIBRARIES(${FE} ${UTILLIB} ${SORTLIB} ${REMLLIB} ${BLUPLIB})
TARGET_LINK_LIBRARIES(${SEL} ${UTILLIB} ${SORTLIB} ${RNGLIB} ${BLUPLIB} ${REMLLIB} ${EVOLLIB})

# Uncomment if you need to link to BLAS and LAPACK
TARGET_LINK_LIBRARIES(${analysis} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(${SA} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(${SA2} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(${SEL} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
# Uncomment if you have parallization
#IF(USE_OPENMP)
#    SET_TARGET_PROPERTIES(${RRR} PROPERTIES
#                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
#                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
#ELSEIF(USE_MPI)
#    SET_TARGET_PROPERTIES(${RRR} PROPERTIES
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
INSTALL(TARGETS ${analysis} RUNTIME DESTINATION bin)
