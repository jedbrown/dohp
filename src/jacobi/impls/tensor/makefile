ALL : lib

DIRS =

SOURCEC = tensor.c ruletopo.c efstopo.c
SOURCEH = tensor.h inlinepoly.h polylib.h
OBJSC = $(SOURCEC:.c=.o)
LIBBASE = libdohp
LOCDIR = src/jacobi/impls/tensor/
EXAMPLES_C = blisttest.c

include ${DOHP_DIR}/conf/base

blisttest: blisttest.o  chkopts
	-${CLINKER} -o blisttest blisttest.o  ${PETSC_LIB}
	${RM} blisttest.o
