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
EXTERN dErr dFSCreateExpandedVector(dFS,Vec*);
EXTERN dErr dFSCreateGlobalVector(dFS,Vec*);
EXTERN dErr dFSGlobalToExpandedBegin(dFS,Vec,InsertMode,Vec);
EXTERN dErr dFSGlobalToExpandedEnd(dFS,Vec,InsertMode,Vec);
EXTERN dErr dFSExpandedToGlobal(dFS,Vec,InsertMode,Vec);
EXTERN dErr dFSExpandedToGlobalBegin(dFS,Vec,InsertMode,Vec);
EXTERN dErr dFSExpandedToGlobalEnd(dFS,Vec,InsertMode,Vec);
EXTERN dErr dFSGetElements(dFS,dInt*,dInt**,s_dRule**,s_dEFS**,dInt**,dReal(**)[3]);
EXTERN dErr dFSRestoreElements(dFS,dInt*,dInt**,s_dRule**,s_dEFS**,dInt**,dReal(**)[3]);
EXTERN dErr dFSGetWorkspace(dFS,dReal(**)[3],dReal(**)[3][3],dReal**,dScalar**,dScalar**,dScalar**,dScalar**);
EXTERN dErr dFSRestoreWorkspace(dFS,dReal(**)[3],dReal(**)[3][3],dReal**,dScalar**,dScalar**,dScalar**,dScalar**);
EXTERN dErr dFSMatSetValuesExpanded(dFS,Mat,dInt,const dInt[],dInt,const dInt[],const dScalar[],InsertMode);
EXTERN dErr dFSGetMatrix(dFS,const MatType,Mat*);
EXTERN dErr dFSBuildSpace(dFS);
EXTERN dErr dFSGetBoundaryType(dFS,dInt,const dMeshEH[],dBdyType[]);
EXTERN dErr dFSCreateSubspace(dFS,dInt,dFSBoundary,dFS*);

EXTERN dErr dFSDestroy(dFS);
EXTERN dErr dFSView(dFS,PetscViewer);
#define dFSRegisterDynamic(a,b,c,d) dFSRegister(a,b,c,d)
EXTERN dErr dFSRegister(const char[],const char[],const char[],dErr(*)(dFS));
EXTERN dErr dFSRegisterAll(const char[]);
EXTERN dErr dFSInitializePackage(const char[]);

EXTERN dErr dQ1HexComputeQuadrature(const dReal x[8][3],dInt *n,const dReal (**qx)[3],const dReal **jw,const dReal **basis,const dReal **deriv);


PETSC_EXTERN_CXX_END
#endif  /* _DOHPFS_H */
