#if !defined(_DOHP_H)
#define _DOHP_H
#include "petscsnes.h"
#include "dohpquotient.h"
PETSC_EXTERN_CXX_BEGIN

// Needed by (automatic) user code to perform element operations

// Requested derivatives
enum DohpDeriv {
  DERIV_BASIS = 0x1,
  DERIV_X     = 0x2,
  DERIV_Y     = 0x4,
  DERIV_Z     = 0x8
};

#define DOHP_NAME_LEN 128

/* Function to define the quadrature order on a given mesh. */
typedef PetscErrorCode (*DohpQuotientFunction1)(const PetscReal*,PetscInt*);
/* Function to define an approximation space given a mesh and a quadrature order. */
typedef PetscErrorCode (*DohpMFSFunction1)(const PetscReal*,const PetscInt*,PetscInt*);

/* Central user-visible objects. */
typedef struct _p_DohpDM *DohpDM;
typedef struct _p_DohpMFS *DohpMFS;
typedef struct _p_DohpWF *DohpWF;
typedef struct _p_DohpBlock* DohpBlock;

typedef struct _p_DohpEFS* DohpEFS;

/* Quadrature rules.  Currently we support tensor-product type Gauss, Gauss-Lobatto, and
* Gauss-Radau.  Tensor product spaces and quadrature rules are central to the efficiency of the
* method, but the implementation should become private. */
typedef struct {
  PetscReal   *coord;
  PetscReal   *weight;
  PetscInt     size;
} DohpRule_Line;

typedef struct {
  DohpRule_Line l[2];
} DohpRule_Quad;

typedef struct {
  DohpRule_Line l[3];
} DohpRule_Hex;

/* Element coordinate mappings, these should become private eventually. */

/* Example implementation, parametric maps normally need *much* more space since the Jacobian has
* different values at each quadrature point. */
typedef struct {
  PetscReal jac[9];
  PetscReal jinv[9];
  PetscReal jdet;
} DohpEMap_Affine3;

/* We can do a poor man's parametric map by using vertex coordinates and computing Jacobians at
* quadrature points on the fly.  For this, we need only store the vertex coordinates in interlaced
* ordering. */
typedef struct {
  PetscReal vtx[2*3];
} DohpEMap_Line;

typedef struct {
  PetscReal vtx[4*3];
} DohpEMap_Quad;

typedef struct {
  PetscReal vtx[8*3];
} DohpEMap_Hex;

/* Computing Jacobians on the fly is likely to be slow for matrix-free computations because extra
* derivatives and 3x3 inverses must be calculated at every quadrature point during every mat-vec
* product.  To overcome this, we can also store the values of the Jacobian, inverse, and determinant
* at all quadrature points.  The we need only map these values into the output vectors during
* element operations. */
typedef struct {
  PetscReal vtx[8*3];
  PetscReal *jac,*jinv,jdet;
} DohpEMap_HexStored;

/* Element basis functions.  Should become private. */

/* We build these for a range of sizes and quadrature sizes.  Such a table can be compiled in or
* created at runtime depending on the size range needed. */
typedef struct {
  PetscScalar *basis;  // (size*qsize), basis[i*size+j] = phi_j(q_i)
  PetscScalar *deriv;  // (size*qsize), deriv[i*size+j] = phi_j'(q_i)
  PetscReal   *ncoord; // (size), nodes of Lagrange polynomial
  PetscInt     size;
} DohpBase;

/* At each element, we just need to store pointers to the parts of the tensor product. */
typedef struct {
  DohpBase *l;
} DohpElem_Line;

typedef struct {
  DohpBase *l[2];
} DohpElem_Quad;

typedef struct {
  DohpBase *l[3];
} DohpElem_Hex;


/* DohpDM: the distributed domain manager consists of
* - 1  DohpWF: continuum statement of the problem
* - 1  DohpMesh: domain and associated metadata
* - 1  DohpRule: quadrature rule associated with every element and every boundary facet
* - 1  DohpEMap: element coordinate mapping associated with every element with a quadrature rule
* - 1+ DohpMFS: scalar function space over a subdomain of the mesh
* - 1+ DohpField: may be scalar or vector valued, defined on a MFS
* */
EXTERN PetscErrorCode DohpDMCreate(MPI_Comm,DohpDM);
EXTERN PetscErrorCode DohpDMSetMesh(DohpDM,DohpMesh);
EXTERN PetscErrorCode DohpDMGetRule(DohpDM,DohpRule_Hex**);
EXTERN PetscErrorCode DohpDMGetMFS(DohpDM,const char[],DohpMFS*);
EXTERN PetscErrorCode DohpDMCreateMFS(DohpDM,const char[],DohpMFS*);
EXTERN PetscErrorCode DohpDMGetMesh(DohpDM,DohpMesh*);
EXTERN PetscErrorCode DohpDMAddField(DohpDM,const char[],DohpMFS,PetscInt);
EXTERN PetscErrorCode DohpDMSetUp(DohpDM);
EXTERN PetscErrorCode DohpDMGetVec(DohpDM,Vec*); /* Get a global vector compatible with the parallel layout of the WF. */
EXTERN PetscErrorCode DohpDMGetVecs(DohpDM,const char*[],Vec*[]);
EXTERN PetscErrorCode DohpDMGetLocalVec(DohpDM,Vec*);
EXTERN PetscErrorCode DohpDMGetLocalVecs(DohpDM,const char*[],Vec*[]);

EXTERN PetscErrorCode DohpBlockGetMatrices(DohpBlock,const MatType,Mat*,Mat*);
EXTERN PetscErrorCode DohpBlockMatMult(Mat,Vec,Vec);

/* EXTERN PetscErrorCode DohpQuotientCreate(DohpQuotient); */
/* EXTERN PetscErrorCode DohpQuotientSetMesh(DohpQuotient,DohpMesh); */
/* EXTERN PetscErrorCode DohpQuotientSetFunction(DohpQuotient,DohpQuotientFunction1); */
/* EXTERN PetscErrorCode DohpQuotientSetup(DohpQuotient); */
EXTERN PetscErrorCode DohpMFSCreate(MPI_Comm,DohpMFS);
/* EXTERN PetscErrorCode DohpMFSSetQuotient(DohpMFS,DohpQuotient); */
EXTERN PetscErrorCode DohpMFSSetFunction(DohpMFS,DohpMFSFunction1);
EXTERN PetscErrorCode DohpMFSSetUp(DohpMFS);
EXTERN PetscErrorCode DohpMFSApplyMinimumRule(DohpMFS,const MeshListInt*);
EXTERN PetscErrorCode DohpMFSSetUpElementBases(DohpMFS,const MeshListInt*);
EXTERN PetscErrorCode DohpMFSSetUpElemFacetProjections(DohpMFS);

/* DohpQuotient: this is not actually an object, just a combination of quadrature rule and element
* map.  The element map operations need to know about the quadrature rule so we keep them together. */
EXTERN PetscErrorCode DohpQuotientComputeElemJac_Hex(const DohpEMap_Hex*,const DohpRule_Hex*,PetscReal*,PetscReal*,PetscReal*);

PETSC_EXTERN_CXX_END
#endif /* _DOHP_H */
