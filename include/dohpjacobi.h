#ifndef _DOHPJACOBI_H
#define _DOHPJACOBI_H
/**
* @file   dohpjacobi.h
* @author Jed Brown <jed@59A2.org>
* @date   Fri Aug 22 19:40:10 2008
*
* @brief  Interface to dJacobi, dEFS, and dRule.
*
* The purpose of these objects is to abstract the generation of basis functions.  dJacobi serves as a cache and
* lookup table for quadrature rules and Basis operations (derivatives and basis evalution).  The implementation may also
* provide unrolled versions for small to modest basis sizes.  In particular, the tightest loops (dofs/node, basis in
* minor dimension, nodes in minor dimension) can be arbitrarily unrolled.  This should reduce overhead of the \a hp
* method for regions of low order.
*
* The idea is that when building a #dMFS we would like the maximum possible sharing and flexibility.  The #dJacobi
* assembles a #dRule and #dEFS for use by any part of the #dMFS.  All the real data (like interpolation and
* differentiation matrices) are shared.
*
*/

#include "dohptype.h"
#include "petscmat.h"

PETSC_EXTERN_CXX_BEGIN

extern PetscCookie dJACOBI_COOKIE;

/**
* Handle for manipulating EFS objects.  The EFS are stored directly in arrays so other components (like dFS) will have
* to be able to see the struct definition in order to know the size.  This is okay since the objects are so simple that
* they should never change.  Note that implementations must have the same size as the opaque prototype (this might
* require a union in some cases).  The tensor products are always okay.
*/
typedef struct s_dEFS *dEFS;

/**
* Handle for manipulating dRule objects (quadrature rules)
*/
typedef struct s_dRule *dRule;

/**
* The vtable is allocated by the dJacobi, one for each element topology and amount of unrolling.  We could move the
* vtable here if we were worried about inlining, but performance studies show that it doesn't actually matter.
*/
#define dEFSHEADER                              \
  struct _dEFSOps *ops;                         \
  dRule            rule

/** This struct should only be referenced in code that holds the array. */
typedef struct s_dEFS {
  dEFSHEADER;
  void *opaque[3];
} s_dEFS;


#define dRuleHEADER                             \
  struct _dRuleOps *ops

/** This struct should only be referenced in code that holds the array. */
typedef struct s_dRule {
  dRuleHEADER;
  void *opaque[3];
} s_dRule;

/**
* Indicates whether or not to apply the transpose of a interpolation/derivative matrix.
*
*/
typedef enum {
  dTRANSPOSE_NO=30001,
  dTRANSPOSE_YES
} dTransposeMode;

/**
* Indicates what type of basis operation to do, see dEFSApply().
*
*/
typedef enum {
  dAPPLY_INTERP=40001,
  dAPPLY_INTERP_TRANSPOSE,
  dAPPLY_GRAD,
  dAPPLY_GRAD_TRANSPOSE,
  dAPPLY_SYMGRAD,
  dAPPLY_SYMGRAD_TRANSPOSE
} dApplyMode;

/**
* Handle for setting up #dRule and #dEFS contexts.
*
*/
typedef struct p_dJacobi *dJacobi;

#define dMESHADJACENCY_HAS_CONNECTIVITY 1
typedef struct _p_dMeshAdjacency *dMeshAdjacency;
struct _p_dMeshAdjacency {
  dMeshESH set;
  dMeshTag indexTag;
  dInt nents;
  dInt toff[5];
  dInt *adjoff,*adjind,*adjperm;
#if defined(dMESHADJACENCY_HAS_CONNECTIVITY)
  dInt *connoff;
  dMeshEH *conn;
#endif
  dEntTopology *topo;
  dMeshEH *ents;
};

#define dJacobiType char *
#define dJACOBI_TENSOR "tensor"

EXTERN dErr dJacobiCreate(MPI_Comm,dJacobi*);
EXTERN dErr dJacobiSetType(dJacobi,dJacobiType);
EXTERN dErr dJacobiSetFromOptions(dJacobi);
EXTERN dErr dJacobiSetUp(dJacobi);
EXTERN dErr dJacobiDestroy(dJacobi);
EXTERN dErr dJacobiView(dJacobi,PetscViewer);
#define dJacobiRegisterDynamic(a,b,c,d) dJacobiRegister(a,b,c,d)
EXTERN dErr dJacobiRegister(const char[],const char[],const char[],dErr(*)(dJacobi));
EXTERN dErr dJacobiRegisterAll(const char[]);
EXTERN dErr dJacobiInitializePackage(const char[]);

EXTERN dErr dJacobiSetDegrees(dJacobi,dInt,dInt);
EXTERN dErr dJacobiGetRule(dJacobi,dInt,const dEntTopology[],const dInt[],dRule);
EXTERN dErr dJacobiGetEFS(dJacobi,dInt,const dEntTopology[],const dInt[],dRule,dEFS);

EXTERN dErr dRuleView(dRule rule,PetscViewer);
EXTERN dErr dRuleGetSize(dRule rule,dInt *dim,dInt *nnodes);
EXTERN dErr dRuleGetNodeWeight(dRule rule,dReal *coord,dReal *weight);
EXTERN dErr dRuleGetTensorNodeWeight(dRule rule,dInt *dim,dInt *nnodes,const dReal *coord[],const dReal *weight[]);
EXTERN dErr dRuleComputeGeometry(dRule rule,const dReal vtx[restrict][3],dReal[restrict][3],dReal jinv[restrict][3][3],dReal jdet[restrict]);

EXTERN dErr dEFSView(dEFS efs,PetscViewer viewer);
EXTERN dErr dEFSGetSizes(dEFS efs,dInt*,dInt *inodes,dInt *total);
EXTERN dErr dEFSGetTensorNodes(dEFS,dInt*,dInt*,dReal**,dReal**,const dReal**,const dReal**);
EXTERN dErr dEFSGetGlobalCoordinates(dEFS,const dReal vtx[restrict][3],dInt*,dInt[3],dReal(*)[3]);
EXTERN dErr dEFSGetRule(dEFS efs,dRule *rule);
EXTERN dErr dEFSApply(dEFS,const dReal[],dInt,const dScalar[],dScalar[restrict],dApplyMode,InsertMode);
EXTERN dErr dJacobiPropogateDown(dJacobi,dMeshAdjacency,dInt[]);
EXTERN dErr dJacobiGetNodeCount(dJacobi,dInt,const dEntTopology[],const dInt[],dInt[],dInt[]);

EXTERN dErr dJacobiGetConstraintCount(dJacobi,dInt,const dInt[],const dInt[],const dInt[],const dInt[],const dMeshAdjacency,dInt[],dInt[]);
EXTERN dErr dJacobiAddConstraints(dJacobi,dInt,const dInt[],const dInt[],const dInt[],const dInt[],const dMeshAdjacency,Mat,Mat);

PETSC_EXTERN_CXX_END

#endif /* _DOHPJACOBI_H */
