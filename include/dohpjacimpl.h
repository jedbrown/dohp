#ifndef _JACIMPL_H
#define _JACIMPL_H

#include "dohpjacobi.h"

dEXTERN_C_BEGIN

extern dLogEvent dLOG_EFSApply;

/**
* The Ops table is directly included in the dRule struct.  There should only be one ops table for each quadrature order
* (the dQuadrature object manages this), and it's common to have optimized versions of some of these operations.  Even
* if generic functions are used for all the ops, there will never be huge numbers of distinct rule orders, and we'd
* prefer to save the indirection.
**/
struct _dRuleOps {
  dErr (*view)(dRule,dViewer);
  dErr (*getSize)(dRule,dInt*,dInt*); /**< topological dimension of the space, total number of nodes */
  dErr (*getNodeWeight)(dRule,dReal[],dReal[]); /**< nodes and weights in interlaced ordering, arrays must be large enough */
  dErr (*getPatches)(dRule,dInt*,dInt*,const dInt **,const dReal**);
  dErr (*getTensorNodeWeight)(dRule,dInt*,dInt[],const dReal**,const dReal**); /**< topological dimension, number of
                                                                               * nodes in each direction, weights in
                                                                               * each direction.  Does not copy, may not
                                                                               * be implemented.  */
  dErr (*computeGeometry)(dRule,const dReal[restrict][3],dReal[restrict][3],dReal[restrict][3][3],dReal[restrict]);
};

#define dRuleHEADER                             \
  struct _dRuleOps ops

/** This struct is never instantiated. */
struct s_dRule {
  dRuleHEADER;
};

/**
* Operations required for an EFS.  Defined here so that these function calls can be inlined.
*
*/
struct _dEFSOps {
  dErr (*view)(dEFS,dViewer);
  dErr (*getSizes)(dEFS,dInt*,dInt*,dInt*); /**< topological dimension, number of interior nodes, total number of nodes */
  dErr (*getTensorNodes)(dEFS,dInt*,dInt*,dReal**,dReal**,const dReal**,const dReal**);
  dErr (*apply)(dEFS,const dReal[],dInt,const dScalar[],dScalar[],dApplyMode,InsertMode);
  dErr (*getGlobalCoordinates)(dEFS,const dReal(*)[3],dInt*,dInt[],dReal(*)[3]);
  dErr (*getExplicit)(dEFS,const dReal[],dInt*,dInt*,const dReal**,const dReal**);
  dErr (*getExplicitSparse)(dEFS,dInt,dInt,const dInt[],const dReal[],dInt,dInt*,dInt[],dReal[],dReal[]);
};

/**
* Like dRule, the ops table is inlined.
*/
#define dEFSHEADER                             \
  struct _dEFSOps ops;                         \
  dRule           rule

/** This struct is never instantiated */
struct s_dEFS {
  dEFSHEADER;
};

/**
* Generic operations provided by a #dJacobi implementation.
*
*/
struct _dJacobiOps {
  dErr (*SetUp)(dJacobi);
  dErr (*SetFromOptions)(dJacobi);
  dErr (*Destroy)(dJacobi);
  dErr (*View)(dJacobi,dViewer);
  dErr (*PropagateDown)(dJacobi,dMeshAdjacency,dPolynomialOrder[]);
  dErr (*GetEFS)(dJacobi,dInt,const dEntTopology[],const dPolynomialOrder[],const dRule[],dEFS[]);
  dErr (*GetNodeCount)(dJacobi,dInt,const dEntTopology[],const dPolynomialOrder[],dInt[],dInt[]);
  dErr (*GetConstraintCount)(dJacobi,dInt,const dInt[],const dInt[],const dInt[],const dPolynomialOrder[],dMeshAdjacency,dInt[],dInt[]);
  dErr (*AddConstraints)(dJacobi,dInt,const dInt[],const dInt[],const dInt[],const dPolynomialOrder[],dMeshAdjacency,Mat,Mat);
  dErr (*GetQuadrature)(dJacobi,dQuadratureMethod,dQuadrature*);
};

/**
* Private Jacobi table context.
*
*/
struct p_dJacobi {
  PETSCHEADER(struct _dJacobiOps);
  bool setupcalled;
  dQuadrature quad[dQUADRATURE_METHOD_INVALID];
  void *data;                   /**< private implementation context */
};

struct _dQuadratureOps {
  dErr (*View)(dQuadrature,PetscViewer);
  dErr (*GetRule)(dQuadrature,dInt,const dEntTopology[],const dPolynomialOrder[],dRule[]);
  dErr (*GetFacetRule)(dQuadrature,dInt,const dEntTopology[],const dInt[],const dPolynomialOrder[],dRule[]);
  dErr (*SetFromOptions)(dQuadrature);
  dErr (*SetMethod)(dQuadrature,dQuadratureMethod);
  dErr (*Destroy)(dQuadrature);
};

struct p_dQuadrature {
  PETSCHEADER(struct _dQuadratureOps);
  void *data;
};

/* These only need to be visible to dJacobiInitializePackage, but I don't like declaring extern functions in source files. */
extern dErr dJacobiCreate_Tensor(dJacobi);
extern dErr dJacobiCreate_Modal(dJacobi);

extern dErr dQuadratureCreate_Tensor(dQuadrature);

extern dErr dEFSGetExplicitSparse_Basic(dEFS efs,dInt npatches,dInt Q,const dInt qidx[],const dReal cjinv[],dInt eoffset,dInt *P,dInt eidx[],dReal interp[],dReal deriv[]);
dEXTERN_C_END

#endif
