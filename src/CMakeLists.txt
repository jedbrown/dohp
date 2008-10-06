# ADD_SUBDIRECTORY (interface)
#set (dohp_SRCS) # interface/dohp.c)

# ADD_SUBDIRECTORY (quotient)
list (APPEND dohp_SRCS
  quotient/interface/quotient.c
  quotient/impls/gauss/gauss.c)

# ADD_SUBDIRECTORY (jacobi)
list (APPEND dohp_SRCS
  jacobi/interface/jacobi.c
  jacobi/impls/tensor/tensor.c
  jacobi/impls/tensor/ruletopo.c
  jacobi/impls/tensor/efstopo.c)

# ADD_SUBDIRECTORY (fs)
LIST (APPEND dohp_SRCS
  fs/interface/fs.c
  fs/interface/fsreg.c
  fs/impls/cont/cont.c)

# ADD_SUBDIRECTORY (fs/mesh)
LIST (APPEND dohp_SRCS
  fs/mesh/interface/mesh.c
  fs/mesh/interface/meshreg.c
  fs/mesh/impls/pack/pack.c
  fs/mesh/impls/serial/serial.c)

ADD_LIBRARY (dohp ${dohp_SRCS})

if (dohp_BUILD_TESTS)
  add_subdirectory (fs/tests)
endif (dohp_BUILD_TESTS)
