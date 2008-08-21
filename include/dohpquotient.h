#ifndef _DOHPQUOTIENT_H
#define _DOHPQUOTIENT_H

#include "dohpmesh.h"

PETSC_EXTERN_CXX_BEGIN

//typedef char *DohpQuotientType;
#define DohpQuotientType char*
typedef struct _p_DohpQuotient *DohpQuotient;

extern PetscCookie DOHP_QUOTIENT_COOKIE,DOHP_MESH_COOKIE;

#define DohpQuotientGauss "gauss"

EXTERN PetscErrorCode DohpQuotientUpdate(DohpQuotient);
EXTERN PetscErrorCode DohpQuotientCreate(DohpMesh m,DohpESH loc,DohpTag qsizetag,DohpQuotient *inq);
EXTERN PetscErrorCode DohpQuotientSetType(DohpQuotient q,const DohpQuotientType type);
EXTERN PetscErrorCode DohpQuotientSetFromOptions(DohpQuotient q);
EXTERN PetscErrorCode DohpQuotientSetUp(DohpQuotient q);
EXTERN PetscErrorCode DohpQuotientDestroy(DohpQuotient q);
EXTERN PetscErrorCode DohpQuotientRegister(const char sname[],const char path[],const char name[],PetscErrorCode (*function)(DohpQuotient));
EXTERN PetscErrorCode DohpQuotientRegisterAll(const char path[]);
EXTERN PetscErrorCode DohpQuotientGetType(DohpQuotient q,const DohpQuotientType *type);
EXTERN PetscErrorCode DohpQuotientView(DohpQuotient q,PetscViewer viewer);
EXTERN PetscErrorCode DohpQuotientInitializePackage(const char path[]);

PETSC_EXTERN_CXX_END

#endif  /* _DOHPQUOTIENT_H */
