#ifndef _JACIMPL_H
#define _JACIMPL_H

#include "dohpjacobi.h"

PETSC_EXTERN_CXX_BEGIN

extern PetscLogEvent dLOG_RuleComputeGeometry,dLOG_EFSApply;

/**
* There is exactly one #dRule on each element.  The ops table is normally shared across the domain.
*
*/
struct _dRuleOps {
  dErr (*view)(dRule,PetscViewer);
  dErr (*getSize)(dRule,dInt*,dInt*); /**< topological dimension of the space, total number of nodes */
  dErr (*getNodeWeight)(dRule,dReal[],dReal[]); /**< nodes and weights in interlaced ordering, arrays must be large enough */
  dErr (*getTensorNodeWeight)(dRule,dInt*,dInt[],const dReal**,const dReal**); /**< topological dimension, number of
                                                                               * nodes in each direction, weights in
                                                                               * each direction.  Does not copy, may not
                                                                               * be implemented.  */
  dErr (*computeGeometry)(dRule,const dReal[restrict][3],dReal[restrict][3],dReal[restrict][3][3],dReal[restrict]);
};

/**
* Operations required for an EFS.  Defined here so that these function calls can be inlined.
*
*/
struct _dEFSOps {
  dErr (*view)(dEFS,PetscViewer);
  dErr (*getSizes)(dEFS,dInt*,dInt*,dInt*); /**< topological dimension, number of interior nodes, total number of nodes */
  dErr (*getTensorNodes)(dEFS,dInt*,dInt*,dReal**,const dReal**,const dReal**);
  dErr (*apply)(dEFS,const dReal[],dInt,const dScalar[],dScalar[],dApplyMode,InsertMode);
  dErr (*getGlobalCoordinates)(dEFS,const dReal(*)[3],dInt*,dInt[],dReal(*)[3]);
  /**< dofs/node, work length, work, modal values, nodal values */
  dErr (*scatterInt)(dEFS,dInt,dInt,const dScalar[],dScalar[],InsertMode,ScatterMode); /**< dofs/node, offset of interior dofs, array, local array */
  /**
  * @bug It's not yet clear to me how to implement this.
  *
  */
  dErr (*scatterFacet)(dEFS,dEFS,dInt*,dScalar**restrict,const dScalar[],dScalar[],InsertMode,ScatterMode);
};

/**
* Generic operations provided by a #dJacobi implementation.
*
*/
struct _dJacobiOps {
  dErr (*SetUp)(dJacobi);
  dErr (*SetFromOptions)(dJacobi);
  dErr (*Destroy)(dJacobi);
  dErr (*View)(dJacobi,PetscViewer);
  dErr (*GetRuleSize)(dJacobi,dEntTopology,dInt*);
  dErr (*PropogateDown)(dJacobi,const struct dMeshAdjacency*,dInt[]);
  dErr (*GetRule)(dJacobi,dInt,const dEntTopology[],const dInt[],dRule);
  dErr (*GetEFS)(dJacobi,dInt,const dEntTopology[],const dInt[],dRule,dEFS);
  dErr (*GetNodeCount)(dJacobi,dInt,const dEntTopology[],const dInt[],dInt[],dInt[]);
  dErr (*GetConstraintCount)(dJacobi,dInt,const dInt[],const dInt[],const dInt[],const dInt[],const struct dMeshAdjacency*,dInt[],dInt[]);
  dErr (*AddConstraints)(dJacobi,dInt,const dInt[],const dInt[],const dInt[],const dInt[],const struct dMeshAdjacency*,Mat,Mat);
};

/**
* Private Jacobi table context.
*
*/
struct p_dJacobi {
  PETSCHEADER(struct _dJacobiOps);
  dInt basisdegree;         /**< the maximum degree basis functions to be supported */
  dInt ruleexcess;          /**< the amount of over-integration to be supported */
  dBool setupcalled;
  void *impl;                   /**< private implementation context */
};

EXTERN dErr dJacobiCreate_Tensor(dJacobi); /* should really only be visible to dJacobiInitializePackage */

PETSC_EXTERN_CXX_END

#endif
