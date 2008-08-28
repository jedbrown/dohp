#ifndef _FSIMPL_H
#define _FSIMPL_H

#include "dohpjacobi.h"

PETSC_EXTERN_CXX_BEGIN

/**
* There is exactly one #dRule on each element.  The ops table is normally shared across the domain.
* 
*/
struct v_dRuleOps {
  dErr (*view)(dRule*,PetscViewer);
  dErr (*getSize)(dRule*,dInt*,dInt*); /**< topological dimension of the space, total number of nodes */
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
#define dRuleView(rule,view)                                    (*rule->ops->view)(rule,view)
#define dRuleGetSize(rule,dim,nnodes)                           (*rule->ops->getSize)(rule,dim,nnodes)
#define dRuleGetNodeWeight(rule,coord,weight)                   (*rule->ops->getNodeWeight)(rule,coord,weight)
#define dRuleGetTensorNodeWeight(rule,dim,nnodes,coord,weights) (*rule->ops->getTensorNodeWeight)(rule,dim,nnodes,coord,weights)

/**
* Operations required for an EFS.  Defined here so that these function calls can be inlined.
* 
*/
struct v_dEFSOps {
  dErr (*view)(dEFS*,PetscViewer);
  dErr (*getSizes)(dEFS*,dInt*,dInt*,dInt*); /**< topological dimension, number of interior nodes, total number of nodes */
  dErr (*apply)(dEFS*,dInt,dInt*,dScalar**restrict,const dScalar[],dScalar[],dApplyMode,InsertMode);
  /**< dofs/node, work length, work, modal values, nodal values */
  dErr (*scatterInt)(dEFS*,dInt,dInt,const dScalar[],dScalar[],InsertMode,ScatterMode); /**< dofs/node, offset of interior dofs, array, local array */
  /**
  * @bug It's not yet clear to me how to implement this.
  * 
  */
  dErr (*scatterFacet)(dEFS,dEFS,dInt*,dScalar**restrict,const dScalar[],dScalar[],InsertMode,ScatterMode);
};

#define dEFSView(efs,viewer) (*efs->ops->view)(efs,viewer)
#define dEFSGetSizes(efs,inodes,total) (*efs->ops->getSizes)(efs,inodes,total)
#define dEFSApply(efs,dofs,wlen,work,in,out,mtype,imode) (*efs->ops->apply)(ofs,dofs,wlen,work,in,out,mtype,imode)

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
