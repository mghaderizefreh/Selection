########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(REML_src ${SRCREML}/getZGZMats.f90
	     ${SRCREML}/calculateV.f90
	     ${SRCREML}/calculateP.f90
	     ${SRCREML}/calculateLogL.f90
	     ${SRCREML}/calculateAImat.f90
	     ${SRCREML}/calculaterhs.f90
	     ${SRCREML}/updatetheta.f90
	     ${SRCREML}/getEffects.f90
	     ${SRCREML}/reml_module.f90
	     ${SRCREML}/RRiteration.f90
	     ${SRCREML}/RRReml.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${EXE} ${REML_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links reml to the utils library
TARGET_LINK_LIBRARIES(${EXE} ${UTILSLIB})

# Uncomment if you need to link to BLAS and LAPACK
TARGET_LINK_LIBRARIES(${EXE} ${BLAS_LIBRARIES}
                                ${LAPACK_LIBRARIES}
                                ${CMAKE_THREAD_LIBS_INIT})

# Uncomment if you have parallization
#IF(USE_OPENMP)
#    SET_TARGET_PROPERTIES(${EXE} PROPERTIES
#                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
#                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
#ELSEIF(USE_MPI)
#    SET_TARGET_PROPERTIES(${EXE} PROPERTIES
#                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
#                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
#    TARGET_LINK_LIBRARIES(${EXE} ${MPI_Fortran_LIBRARIES})
#ENDIF(USE_OPENMP)

#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${EXE} RUNTIME DESTINATION bin)
