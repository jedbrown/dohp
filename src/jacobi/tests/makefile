ALL : ex1

DIRS = 
LOCDIR = src/jacobi/tests/
EXAMPLES_C = ex1.c

include ${DOHP_DIR}/conf/base

ex1: ex1.o  chkopts
	-@${CLINKER} -o ex1 ex1.o  ${DOHP_LIB}
	${RM} ex1.o

runex1 :

include ${DOHP_DIR}/conf/test
