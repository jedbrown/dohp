cmake_minimum_required (VERSION 2.6.3)
project (Dohp)

list (APPEND CMAKE_MODULE_PATH "${Dohp_SOURCE_DIR}/CMake")

# Normally PETSc is built with MPI, if not, use CC=mpicc, etc
find_package (PETSc REQUIRED)
find_package (ITAPS COMPONENTS MESH GEOM REL)

#find_package (Doxygen)

include_directories (${PETSC_INCLUDES} ${ITAPS_INCLUDES})
add_definitions (${PETSC_DEFINITIONS})

option (dUSE_DEBUG         "Compile Dohp with extra debugging" ON)
option (Dohp_BUILD_TESTS   "Build tests"                       ON)
option (Dohp_BUILD_SANDBOX "Build programs in sandbox"         ON)
option (PEDANTIC_WARNINGS  "Compile with pedantic warnings"    ON)
option (BUILD_SHARED_LIBS  "Build shared Dohp libraries"       ON)

list (APPEND CMAKE_C_FLAGS -std=c99)

if (Dohp_BUILD_TESTS)
  enable_testing ()
endif (Dohp_BUILD_TESTS)

if (PEDANTIC_WARNINGS)
  set (DEFAULT_PEDANTIC_FLAGS "-pedantic -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wbad-function-cast -Wcast-align -Wwrite-strings -Wconversion -Wlogical-op -Wsign-compare -Waggregate-return -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Winline -Wno-long-long -Wmissing-format-attribute -Wmissing-noreturn -Wpacked -Wdisabled-optimization -Wmultichar -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wendif-labels -Wdeclaration-after-statement -Wold-style-definition -Winvalid-pch -Wmissing-field-initializers -Wvariadic-macros -Wunsafe-loop-optimizations -Wvolatile-register-var -Wstrict-aliasing -funit-at-a-time -Wno-sign-conversion")
  #set (DEFAULT_PEDANTIC_FLAGS "-Wunreachable-code -Wfloat-equal -Wc++-compat")
  #set (DEFAULT_PEDANTIC_FLAGS "-pedantic -Wall -Wextra -Winline -Wshadow -Wconversion -Wlogical-op -Wmissing-prototypes -Wvla")
  #set (DEFAULT_PEDANTIC_FLAGS "${DEFAULT_PEDANTIC_FLAGS} -Wno-sign-conversion -Wwrite-strings -Wstrict-aliasing -Wcast-align -fstrict-aliasing")
  #set (DEFAULT_PEDANTIC_FLAGS "${DEFAULT_PEDANTIC_FLAGS} -Wdisabled-optimization -funit-at-a-time")
  #set (DEFAULT_PEDANTIC_FLAGS "${DEFAULT_PEDANTIC_FLAGS} -Wpadded")
  set (PEDANTIC_FLAGS ${DEFAULT_PEDANTIC_FLAGS} CACHE STRING "Compiler flags to enable pedantic warnings")
  add_definitions (${PEDANTIC_FLAGS})
endif (PEDANTIC_WARNINGS)

configure_file (${Dohp_SOURCE_DIR}/dohpconfig.h.in
  ${Dohp_BINARY_DIR}/include/dohpconfig.h @ONLY)

add_custom_target (etags
  etags `find include src sandbox -name *.[ch]`
  WORKING_DIRECTORY ${Dohp_SOURCE_DIR})

include_directories (${Dohp_SOURCE_DIR}/include ${Dohp_BINARY_DIR}/include)

set (Dohp_DEPENDENT_LIBRARIES "${ITAPS_REL_LIBRARIES}" "${ITAPS_MESH_LIBRARIES}" "${ITAPS_GEOM_LIBRARIES}" "${PETSC_LIBRARIES}")
set (Dohp_LIBRARIES dohp ${Dohp_DEPENDENT_LIBRARIES})

set (ARCHIVE_OUTPUT_DIRECTORY ${Dohp_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp archives")
set (LIBRARY_OUTPUT_DIRECTORY ${Dohp_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp libraries")
set (RUNTIME_OUTPUT_DIRECTORY ${Dohp_BINARY_DIR}/bin CACHE PATH "Output directory for Dohp executables")
mark_as_advanced (ARCHIVE_OUTPUT_DIRECTORY LIBRARY_OUTPUT_DIRECTORY RUNTIME_OUTPUT_DIRECTORY)

set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_subdirectory (include)
add_subdirectory (src)

if (Dohp_BUILD_SANDBOX)
  add_subdirectory (${Dohp_SOURCE_DIR}/sandbox)
endif (Dohp_BUILD_SANDBOX)
