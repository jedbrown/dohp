#ifndef _FSIMPL_H
#define _FSIMPL_H

#include "dohpjacobi.h"

PETSC_EXTERN_CXX_BEGIN

/**
* There is exactly one #dRule on each element.  The ops table is normally shared across the domain.
* 
*/
struct v_dRuleOps {
  dErr (*getSize)(dRule*,dInt*,dInt*); /**< topological dimension of the space, number of nodes */
  dErr (*getNodeWeight)(dRule*,dReal[],dReal[]); /**< nodes and weights in interlaced ordering, arrays must be large enough */
  dErr (*getTensorNodeWeight)(dRule*,dInt*,dInt[],const dReal**,const dReal**); /**< topological dimension, number of
                                                                               * nodes in each direction, weights in
                                                                               * each direction.  Does not copy, may not
                                                                               * be implemented.  */
};
struct p_dRule {
  struct v_dRuleOps *ops;
  void              *data;
};
#define dRuleGetSize(rule,dim,nnodes)                           (*rule->ops->getSize)(rule,dim,nnodes)
#define dRuleGetNodeWeight(rule,coord,weight)                   (*rule->ops->getNodeWeight)(rule,coord,weight)
#define dRuleGetTensorNodeWeight(rule,dim,nnodes,coord,weights) (*rule->ops->getTensorNodeWeight)(rule,dim,nnodes,coord,weights)

/**
* Operations required for an EFS.  Defined here so that these function calls can be inlined.
* 
*/
struct v_dEFSOps {
  dErr (*applyBasis)(dEFS,dInt,const dScalar *,dScalar *,dScalar *); /**< dofs/node, modal values, nodal values, work */
  dErr (*applyBasisT)(dEFS,dInt,const dScalar *,dScalar *,dScalar *);
  dErr (*applyDeriv)(dEFS,dInt,const dScalar *,dScalar *,dScalar *);
  dErr (*applyDerivT)(dEFS,dInt,const dScalar *,dScalar *,dScalar *);
  dErr (*getNodes)(dEFS,dInt*,dReal[]); /**< coordinates of the nodal basis, not always implemented */
  dErr (*scatterToInt)(dEFS,dInt,dInt,const dScalar[],dScalar[]); /**< dofs/node, offset of interior dofs, array, local array */
  dErr (*scatterFromInt)(dEFS,dInt,dInt,const dScalar[],dScalar[]);
  /**
  * @bug It's not yet clear to me how to implement these.
  * 
  */
  dErr (*facettoelem)(void *,void *,const dScalar *,dScalar *);
  dErr (*elemtofacet)(void *,void *,const dScalar *,dScalar *);
};

/**
* This is held once for every function space on every element.  This part of the implementation is not really private.
* 
*/
struct p_dEFS {
  struct v_dEFSOps *ops;        /**< The element operations, normally constructed by #dJacobi */
  dRule            *rule;       /**< The rule this EFS was built on, probably redundant */
  void             *data;       /**< Private storage to define the basis operations */
};

/**
* Generic operations provided by a #dJacobi implementation.
* 
*/
struct _dJacobiOps {
  dErr (*setup)(dJacobi);
  dErr (*setfromoptions)(dJacobi);
  dErr (*destroy)(dJacobi);
  dErr (*view)(dJacobi,PetscViewer);
  //dErr (*getrule)(dJacobi,dInt,const dInt[],dRule*,dInt*);            /**< put a dRule into the output buffer */
  //dErr (*getefs)(dJacobi,dInt,const dInt[],const dInt[],dEFS*,dInt*); /**< put a dEFS into the output buffer */
  dErr (*getrulesize)(dJacobi,dTopology,dInt*);
  dErr (*getrule)(dJacobi jac,dTopology top,const dInt rsize[],dRule *rule,void **base,dInt *index);
  dErr (*getefs)(dJacobi jac,dTopology top,const dInt bsize[],dRule *rule,dEFS *efs,void **base,dInt *index);
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
