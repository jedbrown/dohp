include ${HOME}/usr/lib/iMesh-Defs.inc
include ${PETSC_DIR}/conf/base

CFLAGS= $(iMesh_INCLUDES)

connect : connect.o
	@-echo 'Linking connect'
	@-$(CLINKER) -o $@ $< ${iMesh_LIBS} ${PETSC_LIB}
dohpenhance : dohpenhance.o
	@-echo 'Linking dohpenhance'
	@-$(CLINKER) -o $@ $< ${iMesh_LIBS} ${PETSC_LIB}
dohpblock : dohpblock.o dohpmesh.h
	@-echo 'Linking dohpblock'
	@-$(CLINKER) -o $@ $< ${iMesh_LIBS} ${PETSC_LIB}
