#if !defined(_DOHP_H)
#define _DOHP_H

#include "petscsnes.h"
#include "dohpquotient.h"
#include "dohpjacobi.h"

PETSC_EXTERN_CXX_BEGIN

// Needed by (automatic) user code to perform element operations

// Requested derivatives
enum DohpDeriv {
  DERIV_BASIS = 0x1,
  DERIV_X     = 0x2,
  DERIV_Y     = 0x4,
  DERIV_Z     = 0x8
};

/* Function to define the quadrature order on a given mesh. */
typedef dErr (*dQuotientFunction1)(const dReal*,dInt*);
/* Function to define an approximation space given a mesh and a quadrature order. */
typedef dErr (*DohpMFSFunction1)(const dReal*,const dInt*,dInt*);

/* Central user-visible objects. */
typedef struct p_dDM *DohpDM;
typedef struct p_dMFS *DohpMFS;
typedef struct p_dWF *DohpWF;
typedef struct p_dBlock *DohpBlock;


/* Element coordinate mappings, these should become private eventually. */

/* Example implementation, parametric maps normally need *much* more space since the Jacobian has
* different values at each quadrature point. */
typedef struct {
  dReal jac[9];
  dReal jinv[9];
  dReal jdet;
} DohpEMap_Affine3;

/* We can do a poor man's parametric map by using vertex coordinates and computing Jacobians at
* quadrature points on the fly.  For this, we need only store the vertex coordinates in interleaved
* ordering. */
typedef struct {
  dReal vtx[2*3];
} DohpEMap_Line;

typedef struct {
  dReal vtx[4*3];
} DohpEMap_Quad;

typedef struct {
  dReal vtx[8*3];
} DohpEMap_Hex;

/* Computing Jacobians on the fly is likely to be slow for matrix-free computations because extra
* derivatives and 3x3 inverses must be calculated at every quadrature point during every mat-vec
* product.  To overcome this, we can also store the values of the Jacobian, inverse, and determinant
* at all quadrature points.  The we need only map these values into the output vectors during
* element operations. */
typedef struct {
  dReal vtx[8*3];
  dReal *jac,*jinv,jdet;
} DohpEMap_HexStored;

/* Element basis functions.  Should become private. */


/* DohpDM: the distributed domain manager consists of
* - 1  DohpWF: continuum statement of the problem
* - 1  dMesh: domain and associated metadata
* - 1  dRule: quadrature rule associated with every element and every boundary facet
* - 1  DohpEMap: element coordinate mapping associated with every element with a quadrature rule
* - 1+ DohpMFS: scalar function space over a subdomain of the mesh
* - 1+ DohpField: may be scalar or vector valued, defined on a MFS
* */
EXTERN dErr DohpDMCreate(MPI_Comm,DohpDM);
EXTERN dErr DohpDMSetMesh(DohpDM,dMesh);
EXTERN dErr DohpDMGetMFS(DohpDM,const char[],DohpMFS*);
EXTERN dErr DohpDMCreateMFS(DohpDM,const char[],DohpMFS*);
EXTERN dErr DohpDMGetMesh(DohpDM,dMesh*);
EXTERN dErr DohpDMAddField(DohpDM,const char[],DohpMFS,dInt);
EXTERN dErr DohpDMSetUp(DohpDM);
EXTERN dErr DohpDMGetVec(DohpDM,Vec*); /* Get a global vector compatible with the parallel layout of the WF. */
EXTERN dErr DohpDMGetVecs(DohpDM,const char*[],Vec*[]);
EXTERN dErr DohpDMGetLocalVec(DohpDM,Vec*);
EXTERN dErr DohpDMGetLocalVecs(DohpDM,const char*[],Vec*[]);

EXTERN dErr DohpBlockGetMatrices(DohpBlock,const MatType,Mat*,Mat*);
EXTERN dErr DohpBlockMatMult(Mat,Vec,Vec);

/* EXTERN dErr dQuotientCreate(dQuotient); */
/* EXTERN dErr dQuotientSetMesh(dQuotient,dMesh); */
/* EXTERN dErr dQuotientSetFunction(dQuotient,dQuotientFunction1); */
/* EXTERN dErr dQuotientSetup(dQuotient); */
EXTERN dErr DohpMFSCreate(MPI_Comm,DohpMFS);
/* EXTERN dErr DohpMFSSetQuotient(DohpMFS,dQuotient); */
EXTERN dErr DohpMFSSetFunction(DohpMFS,DohpMFSFunction1);
EXTERN dErr DohpMFSSetUp(DohpMFS);
EXTERN dErr DohpMFSApplyMinimumRule(DohpMFS,const MeshListInt*);
EXTERN dErr DohpMFSSetUpElementBases(DohpMFS,const MeshListInt*);
EXTERN dErr DohpMFSSetUpElemFacetProjections(DohpMFS);
EXTERN dErr DohpMFSSetUpBoundaryTypes(DohpMFS mfs);

PETSC_EXTERN_CXX_END
#endif /* _DOHP_H */
