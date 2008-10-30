#ifndef _DOHPFS_H
#define _DOHPFS_H
/**
* @file   dohpfs.h
* @author Jed Brown <jed@59A2.org>
* @date   Sun Sep  7 17:46:54 2008
*
* @brief  The function space object
*
*
*/

#include "dohpmesh.h"
#include "dohpjacobi.h"
#include "dohpquotient.h"
#include "petscpf.h"

PETSC_EXTERN_CXX_BEGIN

typedef struct _p_dFS *dFS;

typedef struct _p_dFSBoundary *dFSBoundary;

#define dFSType char *

#define dFSCONT "cont"

extern PetscCookie dFS_COOKIE;

EXTERN dErr dFSCreate(MPI_Comm,dFS*);
EXTERN dErr dFSSetMesh(dFS,dMesh,dMeshESH); /* mesh, active set */
EXTERN dErr dFSSetRuleTag(dFS,dJacobi,dMeshTag);
EXTERN dErr dFSSetDegree(dFS,dJacobi,dMeshTag);
EXTERN dErr dFSAddBdy(dFS,const char*,dMeshESH,dMeshTag,dBool,PF); /* name, facets, orientation tag, flip orientation?, normal -> constraints */
EXTERN dErr dFSSetFromOptions(dFS);
EXTERN dErr dFSSetType(dFS,const dFSType);
EXTERN dErr dFSCreateLocalVector(dFS,Vec*);
EXTERN dErr dFSCreateGlobalVector(dFS,Vec*);
EXTERN dErr dFSBuildSpace(dFS);
EXTERN dErr dFSGetBoundaryType(dFS,dInt,const dMeshEH[],dBdyType[]);
EXTERN dErr dFSCreateSubspace(dFS,dInt,dFSBoundary,dFS*);

EXTERN dErr dFSDestroy(dFS);
EXTERN dErr dFSView(dFS,PetscViewer);
#define dFSRegisterDynamic(a,b,c,d) dFSRegister(a,b,c,d)
EXTERN dErr dFSRegister(const char[],const char[],const char[],dErr(*)(dFS));
EXTERN dErr dFSRegisterAll(const char[]);
EXTERN dErr dFSInitializePackage(const char[]);

PETSC_EXTERN_CXX_END
#endif  /* _DOHPFS_H */
