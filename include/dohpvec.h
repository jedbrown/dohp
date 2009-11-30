#ifndef _DOHPVEC_H
#define _DOHPVEC_H

#include <petscvec.h>
#include "dohptype.h"

dEXTERN_C_BEGIN

#define VECDOHP "dohp"

extern dErr VecDohpGetClosure(Vec,Vec*);
extern dErr VecDohpRestoreClosure(Vec,Vec*);
extern dErr VecCreateDohp(MPI_Comm,dInt,dInt,dInt,dInt,const dInt[],Vec*);

extern dErr VecDohpCreateDirichletCache(Vec gvec,Vec *dcache,VecScatter *dscat);

dEXTERN_C_BEGIN

#endif
