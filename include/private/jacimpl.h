#ifndef _JACIMPL_H
#define _JACIMPL_H

#include "dohpjacobi.h"

PETSC_EXTERN_CXX_BEGIN

/**
* There is exactly one #dRule on each element.  The ops table is normally shared across the domain.
*
*/
struct v_dRuleOps {
  dErr (*view)(dRule,PetscViewer);
  dErr (*getSize)(dRule,dInt*,dInt*); /**< topological dimension of the space, total number of nodes */
  dErr (*getNodeWeight)(dRule,dReal[],dReal[]); /**< nodes and weights in interlaced ordering, arrays must be large enough */
  dErr (*getTensorNodeWeight)(dRule,dInt*,dInt[],const dReal**,const dReal**); /**< topological dimension, number of
                                                                               * nodes in each direction, weights in
                                                                               * each direction.  Does not copy, may not
                                                                               * be implemented.  */
};

#define dRuleHEADER                             \
  struct v_dRuleOps *ops

struct p_dRule {
  dRuleHEADER;
};

/**
* Operations required for an EFS.  Defined here so that these function calls can be inlined.
*
*/
struct v_dEFSOps {
  dErr (*view)(dEFS,PetscViewer);
  dErr (*getSizes)(dEFS,dInt*,dInt*,dInt*); /**< topological dimension, number of interior nodes, total number of nodes */
  dErr (*getTensorNodes)(dEFS,dInt*,dInt*,dReal**);
  dErr (*apply)(dEFS,dInt,dInt*,dScalar**restrict,const dScalar[],dScalar[],dApplyMode,InsertMode);
  /**< dofs/node, work length, work, modal values, nodal values */
  dErr (*propogatedown)(const dInt a[],const dMeshEH ev[],const dInt f[],const dEntTopology ftopo[],const dMeshEH fv[],dInt af[]);
  dErr (*scatterInt)(dEFS,dInt,dInt,const dScalar[],dScalar[],InsertMode,ScatterMode); /**< dofs/node, offset of interior dofs, array, local array */
  /**
  * @bug It's not yet clear to me how to implement this.
  *
  */
  dErr (*scatterFacet)(dEFS,dEFS,dInt*,dScalar**restrict,const dScalar[],dScalar[],InsertMode,ScatterMode);
};

/**
* This is held once for every function space on every element.  This part of the implementation is not really private.
*
*/
#define dEFSHEADER                              \
  struct v_dEFSOps *ops;                        \
  dRule             rule

struct p_dEFS {
  dEFSHEADER;
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
  dErr (*getrulesize)(dJacobi,dEntTopology,dInt*);
  dErr (*getrule)(dJacobi jac,dEntTopology top,const dInt rsize[],dRule *rule,void **base,dInt *index);
  dErr (*getefs)(dJacobi jac,dEntTopology top,const dInt bsize[],dRule rule,dEFS *efs,void **base,dInt *index);
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
