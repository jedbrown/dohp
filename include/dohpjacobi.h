#ifndef _DOHPJACOBI_H
#define _DOHPJACOBI_H
/**
* @file   dohpjacobi.h
* @author Jed Brown <jed@59A2.org>
* @date   Fri Aug 22 19:40:10 2008
* 
* @brief  Interface to dJacobi, dBasis, and dRule.
*
* The purpose of these objects is to abstract the generation of basis functions.  dJacobi serves as a cache and
* lookup table for quadrature rules and Basis operations (derivatives and basis evalution).  The implementation may also
* provide unrolled versions for small to modest basis sizes.  In particular, the tightest loops (dofs/node, basis in
* minor dimension, nodes in minor dimension) can be arbitrarily unrolled.  This should reduce overhead of the \a hp
* method for regions of low order.
* 
*/

#include "dohptype.h"
#include "private/petscimpl.h"

PETSC_EXTERN_CXX_BEGIN

/**
* This isn't actually an opaque structure.  Note that it is necessarily \c const since it will normally point to shared
* memory space.
* 
*/
typedef const struct m_dRule * dRule;

/**
* Minimal struct to hold a 1-D quadrature rule.
*
* Although this does not support any functions, we still make it opaque.  This is so we can work with a dRule by
* passing pointers and it will look like any other Petsc/Dohp object.
* 
*/
struct m_dRule {
  const dReal *node;              /**< nodal coordinates on the reference element [-1,1] */
  const dReal *weight;            /**< quadrature weights at the nodes  */
  dInt size;                /**< total number of points */
};


/**
* The handle for basis evaluation.
*
* The basis must always be compatible with a quadrature rule.  Eventually we'll have wrapper functions, but for now,
* just use the member functions.  Note that it is \c const since it will normally point to shared space.
* 
*/
typedef const struct m_dBasis *dBasis;

/**
* Holds interpolation and differentiation context.
*
* This is implemented as pointers to the matrices and (possibly) unrolled versions of functions to apply those vectors
* across a 3D block of values.  This is the mutable handle, it should only be used during setup.
* 
*/
struct m_dBasis {
  dErr (*mult)(dInt,dInt,dInt,dInt,dInt,dInt,dInt,const dReal[],const dReal[],const dReal[],dReal[]);
  dErr (*multtrans)(dInt,dInt,dInt,dInt,dInt,dInt,dInt,const dReal[],const dReal[],const dReal[],dReal[]);
  const dReal *bmat;        /**< pointer to the basis matrix */
  const dReal *dmat;        /**< pointer to the derivative matrix */
  dRule rule;                /**< the rule on which the basis has been constructed, redundant */
  dInt bsize;               /**< basis size, should be redundant */
  dInt qsize;               /**< quadrature size, should be redundant */
};

/**
* Handle for setting up #dRule and #dBasis.
* 
*/
typedef struct p_dJacobi *dJacobi;


/**
* Generic operations provided by a #dJacobi implementation.
* 
*/
struct _dJacobiOps {
  dErr (*setup)(dJacobi);
  dErr (*setfromoptions)(dJacobi);
  dErr (*destroy)(dJacobi);
  dErr (*view)(dJacobi,PetscViewer);
  dErr (*getrule)(dJacobi,dInt,dRule*);
  dErr (*getbasis)(dJacobi,dInt,dInt,dBasis*);
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

#define dJacobiType char *
#define dJACOBI_LGL "lgl"

EXTERN dErr dJacobiCreate_LGL(dJacobi);

EXTERN dErr dJacobiCreate(MPI_Comm,dJacobi*);
EXTERN dErr dJacobiSetType(dJacobi,dJacobiType);
EXTERN dErr dJacobiSetFromOptions(dJacobi);
EXTERN dErr dJacobiSetUp(dJacobi);
EXTERN dErr dJacobiSetFromOptions(dJacobi);
EXTERN dErr dJacobiDestroy(dJacobi);
EXTERN dErr dJacobiView(dJacobi,PetscViewer);
#define dJacobiRegisterDynamic(a,b,c,d) dJacobiRegister(a,b,c,d)
EXTERN dErr dJacobiRegister(const char[],const char[],const char[],dErr(*)(dJacobi));
EXTERN dErr dJacobiRegisterAll(const char[]);
EXTERN dErr dJacobiInitializePackage(const char[]);

EXTERN dErr dJacobiSetDegrees(dJacobi,dInt,dInt);
EXTERN dErr dJacobiGetRule(dJacobi,dInt,dRule*);
EXTERN dErr dJacobiGetBasis(dJacobi,dInt,dInt,dBasis*);


PETSC_EXTERN_CXX_END

#endif /* _DOHPJACOBI_H */
