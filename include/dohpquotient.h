#ifndef _DOHPQUOTIENT_H
#define _DOHPQUOTIENT_H

#include "dohpmesh.h"

dEXTERN_C_BEGIN

#define dQuotientType char*
typedef struct p_dQuotient *dQuotient;

extern dCookie dQUOTIENT_COOKIE;
typedef dErr (*dQuotientSetDegreeFunc)(dQuotient,void*,dInt,dInt[]);

#define dQuotientGauss "gauss"

extern dErr dQuotientUpdate(dQuotient);
extern dErr dQuotientCreate(dMesh m,dMeshESH loc,dMeshTag qsizetag,dQuotient *inq);
extern dErr dQuotientSetType(dQuotient q,const dQuotientType type);
extern dErr dQuotientSetFromOptions(dQuotient q);
extern dErr dQuotientSetUp(dQuotient q);
extern dErr dQuotientGetMesh(dQuotient,dMesh*);

extern dErr dQuotientDestroy(dQuotient q);
extern dErr dQuotientRegister(const char sname[],const char path[],const char name[],dErr (*function)(dQuotient));
extern dErr dQuotientRegisterAll(const char path[]);
extern dErr dQuotientGetType(dQuotient q,const dQuotientType *type);
extern dErr dQuotientView(dQuotient q,dViewer viewer);
extern dErr dQuotientInitializePackage(const char path[]);
extern dErr dQuotientSetSetDegree(dQuotient q,dQuotientSetDegreeFunc func,void *ctx);
extern dErr dQuotientSetDegreeConst(dQuotient q,void *vval,dInt n,dInt *degree);
extern dErr dQuotientGetArrRule(dQuotient q,dInt,dMeshEH[],dRule**);

dEXTERN_C_END

#endif  /* _DOHPQUOTIENT_H */
