########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(RRR_src ${MAIN}/doRRReml.f90
)

SET(RRB_src ${MAIN}/doRRBlup.f90
)

SET(STR_src ${MAIN}/doSTReml.f90
)

SET(STB_src ${MAIN}/doSTBlup.f90
)

SET(FE_src ${MAIN}/fixedEff.f90
)

SET(SA_src ${MAIN}/sample.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${RRR} ${RRR_src})
ADD_EXECUTABLE(${RRB} ${RRB_src})
ADD_EXECUTABLE(${STR} ${STR_src})
ADD_EXECUTABLE(${STB} ${STB_src})
ADD_EXECUTABLE(${SA} ${SA_src})
ADD_EXECUTABLE(${FE} ${FE_src})
#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links exes to the their libraries
TARGET_LINK_LIBRARIES(${RRR} ${UTILLIB} ${BLUPLIB} ${REMLLIB})
TARGET_LINK_LIBRARIES(${STR} ${UTILLIB} ${BLUPLIB} ${REMLLIB})
TARGET_LINK_LIBRARIES(${RRB} ${UTILLIB} ${BLUPLIB})
TARGET_LINK_LIBRARIES(${STB} ${UTILLIB} ${BLUPLIB})
TARGET_LINK_LIBRARIES(${SA} ${UTILLIB} ${SORTLIB} ${RNGLIB})
TARGET_LINK_LIBRARIES(${FE} ${UTILLIB} ${SORTLIB} ${REMLLIB} ${BLUPLIB})

# Uncomment if you need to link to BLAS and LAPACK
TARGET_LINK_LIBRARIES(${RRR} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(${RRB} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(${STR} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(${STB} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(${SA} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

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
INSTALL(TARGETS ${RRR} RUNTIME DESTINATION bin)
