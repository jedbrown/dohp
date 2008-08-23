#ifndef _DOHPQUOTIENT_H
#define _DOHPQUOTIENT_H

#include "dohpmesh.h"

PETSC_EXTERN_CXX_BEGIN

//typedef char *dQuotientType;
#define dQuotientType char*
typedef struct p_dQuotient *dQuotient;

extern PetscCookie dQUOTIENT_COOKIE,dMESH_COOKIE;
typedef dErr (*dQuotientSetDegreeFunc)(dQuotient,void*,dInt,dInt[]);

#define dQuotientGauss "gauss"

EXTERN dErr dQuotientUpdate(dQuotient);
EXTERN dErr dQuotientCreate(dMesh m,DohpESH loc,DohpTag qsizetag,dQuotient *inq);
EXTERN dErr dQuotientSetType(dQuotient q,const dQuotientType type);
EXTERN dErr dQuotientSetFromOptions(dQuotient q);
EXTERN dErr dQuotientSetUp(dQuotient q);
EXTERN dErr dQuotientDestroy(dQuotient q);
EXTERN dErr dQuotientRegister(const char sname[],const char path[],const char name[],dErr (*function)(dQuotient));
EXTERN dErr dQuotientRegisterAll(const char path[]);
EXTERN dErr dQuotientGetType(dQuotient q,const dQuotientType *type);
EXTERN dErr dQuotientView(dQuotient q,PetscViewer viewer);
EXTERN dErr dQuotientInitializePackage(const char path[]);
EXTERN dErr dQuotientSetSetDegree(dQuotient q,dQuotientSetDegreeFunc func,void *ctx);
EXTERN dErr dQuotientSetDegreeConst(dQuotient q,void *vval,dInt n,dInt *degree);

PETSC_EXTERN_CXX_END

#endif  /* _DOHPQUOTIENT_H */
