#ifndef _DOHPVEC_H
#define _DOHPVEC_H

#include "petscvec.h"
#include "dohptype.h"

PETSC_EXTERN_CXX_BEGIN

#define VECDOHP "dohp"

EXTERN dErr VecDohpGetClosure(Vec,Vec*);
EXTERN dErr VecDohpRestoreClosure(Vec,Vec*);
EXTERN dErr VecCreateDohp(MPI_Comm,dInt,dInt,dInt,dInt,const dInt[],Vec*);

PETSC_EXTERN_CXX_BEGIN

#endif
