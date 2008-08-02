ALL: all
LOCDIR = .
DIRS = src include doc
include ${DOHP_DIR}/conf/base

all: build shared

build:
	-@echo "BEGINNING TO COMPILE DOHP LIBRARIES IN ALL DIRECTORIES"
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH} ACTION=libfast  tree 
	${RANLIB} ${DOHP_LIB_DIR}/*.${AR_LIB_SUFFIX}
	-@echo "Completed building DOHP libraries"
	-@echo "========================================="
