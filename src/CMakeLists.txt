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
LIST (APPEND Dohp_SRCS
  fs/interface/fs.c
  fs/interface/fsreg.c
  fs/interface/q1.c
  fs/impls/cont/cont.c)

# ADD_SUBDIRECTORY (fs/mesh)
LIST (APPEND Dohp_SRCS
  fs/mesh/interface/mesh.c
  fs/mesh/interface/meshreg.c
  fs/mesh/impls/pack/pack.c
  fs/mesh/impls/serial/serial.c)

LIST (APPEND Dohp_SRCS
  util/util.c)

ADD_LIBRARY (dohp ${Dohp_SRCS})

if (Dohp_BUILD_TESTS)
  add_subdirectory (fs/tests)
  add_subdirectory (jacobi/tests)
endif (Dohp_BUILD_TESTS)
