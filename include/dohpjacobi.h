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

dEXTERN_C_BEGIN

extern PetscCookie dJACOBI_COOKIE,dQUADRATURE_COOKIE;

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

typedef struct p_dQuadrature *dQuadrature;

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

#define dJacobiType char*
#define dJACOBI_TENSOR "tensor"
#define dJACOBI_MODAL  "modal"

#define dQuadratureType char*
#define dQUADRATURE_TENSOR "tensor"

extern dErr dJacobiCreate(MPI_Comm,dJacobi*);
extern dErr dJacobiSetType(dJacobi,dJacobiType);
extern dErr dJacobiSetFromOptions(dJacobi);
extern dErr dJacobiDestroy(dJacobi);
extern dErr dJacobiView(dJacobi,dViewer);
#define dJacobiRegisterDynamic(a,b,c,d) dJacobiRegister(a,b,c,d)
extern dErr dJacobiRegister(const char[],const char[],const char[],dErr(*)(dJacobi));
extern dErr dJacobiRegisterAll(const char[]);
extern dErr dJacobiInitializePackage(const char[]);

extern dErr dJacobiSetDegrees(dJacobi,dInt,dInt);
extern dErr dJacobiGetEFS(dJacobi,dInt,const dEntTopology[],const dInt[],dRule,dEFS);
extern dErr dJacobiGetQuadrature(dJacobi,dQuadrature*);

extern dErr dRuleView(dRule rule,dViewer);
extern dErr dRuleGetSize(dRule rule,dInt *dim,dInt *nnodes);
extern dErr dRuleGetNodeWeight(dRule rule,dReal *coord,dReal *weight);
extern dErr dRuleGetTensorNodeWeight(dRule rule,dInt *dim,dInt *nnodes,const dReal *coord[],const dReal *weight[]);
extern dErr dRuleComputeGeometry(dRule rule,const dReal vtx[restrict][3],dReal[restrict][3],dReal jinv[restrict][3][3],dReal jdet[restrict]);

extern dErr dEFSView(dEFS efs,dViewer viewer);
extern dErr dEFSGetSizes(dEFS efs,dInt*,dInt *inodes,dInt *total);
extern dErr dEFSGetTensorNodes(dEFS,dInt*,dInt*,dReal**,dReal**,const dReal**,const dReal**);
extern dErr dEFSGetGlobalCoordinates(dEFS,const dReal vtx[restrict][3],dInt*,dInt[3],dReal(*)[3]);
extern dErr dEFSGetRule(dEFS efs,dRule *rule);
extern dErr dEFSApply(dEFS,const dReal[],dInt,const dScalar[],dScalar[restrict],dApplyMode,InsertMode);
extern dErr dJacobiPropogateDown(dJacobi,dMeshAdjacency,dInt[]);
extern dErr dJacobiGetNodeCount(dJacobi,dInt,const dEntTopology[],const dInt[],dInt[],dInt[]);

extern dErr dJacobiGetConstraintCount(dJacobi,dInt,const dInt[],const dInt[],const dInt[],const dInt[],const dMeshAdjacency,dInt[],dInt[]);
extern dErr dJacobiAddConstraints(dJacobi,dInt,const dInt[],const dInt[],const dInt[],const dInt[],const dMeshAdjacency,Mat,Mat);

#define dQuadratureRegisterDynamic(a,b,c,d) dQuadratureRegister(a,b,c,d)
extern dErr dQuadratureRegister(const char[],const char[],const char[],dErr(*)(dQuadrature));
extern dErr dQuadratureRegisterAll(const char[]);

extern dErr dQuadratureCreate(MPI_Comm,dQuadrature*);
extern dErr dQuadratureView(dQuadrature,PetscViewer);
extern dErr dQuadratureSetType(dQuadrature,dQuadratureType);
extern dErr dQuadratureSetFromOptions(dQuadrature);
extern dErr dQuadratureDestroy(dQuadrature);
extern dErr dQuadratureGetRule(dQuadrature,dInt,const dEntTopology[],const dInt[],dRule);

typedef enum {
  dJACOBI_MODAL_P_CONFORMING,
  dJACOBI_MODAL_P_DISCONTINUOUS,
  dJACOBI_MODAL_Q_CONFORMING,
  dJACOBI_MODAL_Q_DISCONTINUOUS,
} dJacobiModalFamily;
extern const char *const dJacobiModalFamilies[];

extern dErr dJacobiModalSetFamily(dJacobi,dJacobiModalFamily);

dEXTERN_C_END

#endif /* _DOHPJACOBI_H */
