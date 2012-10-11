#ifndef _DOHPVEC_H
#define _DOHPVEC_H

#include "dohptype.h"
#include <petscvec.h>

dEXTERN_C_BEGIN

#define VECDOHP "dohp"

extern dErr VecDohpGetClosure(Vec,Vec*);
extern dErr VecDohpRestoreClosure(Vec,Vec*);
extern dErr VecDohpZeroEntries(Vec);
extern dErr VecCreateDohp(MPI_Comm,dInt,dInt,dInt,dInt,const dInt[],Vec*);

extern dErr VecDohpCreateDirichletCache(Vec gvec,Vec *dcache,VecScatter *dscat);

extern dErr VecDohpLoadIntoVector(PetscViewer,const char fieldname[],Vec);

extern dErr VecCreateRedimensioned(Vec X,dInt bs,Vec *Y);
extern dErr VecBlockView(Vec X,dViewer viewer);

dEXTERN_C_END
#endif
