# ADD_SUBDIRECTORY (interface)
#set (Dohp_SRCS) # interface/dohp.c)

# ADD_SUBDIRECTORY (quotient)
list (APPEND Dohp_SRCS
  quotient/interface/quotient.c
  quotient/impls/gauss/gauss.c)

# ADD_SUBDIRECTORY (jacobi)
list (APPEND Dohp_SRCS
  jacobi/interface/jacobi.c
  jacobi/impls/tensor/tensor.c
  jacobi/impls/tensor/ruletopo.c
  jacobi/impls/tensor/efstopo.c)

# ADD_SUBDIRECTORY (fs)
list (APPEND Dohp_SRCS
  fs/interface/fs.c
  fs/interface/fsreg.c
  fs/interface/fsrot.c
  fs/interface/fsgeom.c
  fs/interface/q1.c
  fs/impls/cont/cont.c)

# ADD_SUBDIRECTORY (fs/mesh)
list (APPEND Dohp_SRCS
  fs/mesh/interface/mesh.c
  fs/mesh/interface/meshreg.c
  fs/mesh/impls/pack/pack.c
  fs/mesh/impls/serial/serial.c)

# ADD_SUBDIRECTORY (fs/vec)
list (APPEND Dohp_SRCS
  vec/vecd.c)

list (APPEND Dohp_SRCS
  util/util.c)

add_library (dohp ${Dohp_SRCS})
target_link_libraries (dohp ${Dohp_DEPENDENT_LIBRARIES})
install (TARGETS dohp
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

if (Dohp_BUILD_TESTS)
  add_subdirectory (jacobi/tests)
  add_subdirectory (vec/tests)
  add_subdirectory (fs/tests)
endif (Dohp_BUILD_TESTS)

if (Dohp_BUILD_TOOLS)
  add_subdirectory (tools)
endif (Dohp_BUILD_TOOLS)

