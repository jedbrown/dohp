project (Dohp)

cmake_minimum_required (VERSION 2.6)

list (APPEND CMAKE_MODULE_PATH "${Dohp_SOURCE_DIR}/CMake")

find_package (PETSc REQUIRED)
find_package (iMesh REQUIRED)
# use the PETSc compiler as a hint when looking for an MPI compiler
set (MPI_COMPILER ${PETSC_COMPILER} CACHE FILEPATH "MPI compiler. Used only to detect MPI compilation flags.")
set (MPIEXEC ${PETSC_MPIEXEC} CACHE FILEPATH "Executable for running MPI programs.")
mark_as_advanced(MPI_COMPILER MPI_EXEC)

# Can we get things to stop here, MPI doesn't reconfigure well?
# PETSc knows how to link MPI so we don't need to find it separately.
#find_package (MPI REQUIRED)

# How do we identify whether all libs link against the same MPI?

#find_package (Doxygen)

include_directories (${PETSC_INCLUDE_PATH} ${IMESH_INCLUDE_PATH} ${MPI_INCLUDE_PATH})
add_definitions (${PETSC_DEFINITIONS} ${MPI_COMPILE_FLAGS})

option (Dohp_USE_DEBUG     "Compile with debugging support" ON)
option (Dohp_BUILD_TESTS   "Build tests"                    ON)
option (Dohp_BUILD_SANDBOX "Build programs in sandbox"      ON)
option (PEDANTIC_WARNINGS  "Compile with pedantic warnings" ON)
option (BUILD_SHARED_LIBS  "Build shared Dohp libraries"    ON)

add_definitions (-std=c99)
if (Dohp_USE_DEBUG)
  add_definitions (-DdUSE_DEBUG=1 -g3)
endif (Dohp_USE_DEBUG)

if (PEDANTIC_WARNINGS)
  set (DEFAULT_PEDANTIC_FLAGS "-pedantic -Wall -Wextra -Winline -Wshadow -Wconversion -Wlogical-op -Wmissing-prototypes -Wvla")
  set (DEFAULT_PEDANTIC_FLAGS "${DEFAULT_PEDANTIC_FLAGS} -Wno-sign-conversion -Wwrite-strings -Wstrict-aliasing -Wcast-align -fstrict-aliasing")
  set (DEFAULT_PEDANTIC_FLAGS "${DEFAULT_PEDANTIC_FLAGS} -Wdisabled-optimization -funit-at-a-time")
  #set (DEFAULT_PEDANTIC_FLAGS "${DEFAULT_PEDANTIC_FLAGS} -Wpadded")
  set (PEDANTIC_FLAGS ${DEFAULT_PEDANTIC_FLAGS} CACHE STRING "Compiler flags to enable pedantic warnings")
  add_definitions (${PEDANTIC_FLAGS})
endif (PEDANTIC_WARNINGS)

add_custom_target (etags
  etags `find include src sandbox -name *.[ch]`
  WORKING_DIRECTORY ${Dohp_SOURCE_DIR})

include_directories (${Dohp_SOURCE_DIR} ${Dohp_SOURCE_DIR}/include)

set (Dohp_LIBRARIES dohp ${IMESH_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_LIBRARIES})

add_subdirectory (include)
add_subdirectory (src)

if (Dohp_BUILD_SANDBOX)
  add_subdirectory (${Dohp_SOURCE_DIR}/sandbox)
endif (Dohp_BUILD_SANDBOX)
