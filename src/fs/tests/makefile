ALL :

DIRS = 
LOCDIR = src/fs/tests/
EXAMPLES_C = ex1.c ex2.c

include ${DOHP_DIR}/conf/base

ex1: ex1.o  chkopts
	-@${CLINKER} -o ex1 ex1.o  ${DOHP_LIB}
	${RM} ex1.o
ex2: ex2.o  chkopts
	-@${CLINKER} -o ex2 ex2.o  ${DOHP_LIB}
	${RM} ex2.o

runex1 :

include ${DOHP_DIR}/conf/test
