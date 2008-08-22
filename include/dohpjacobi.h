#ifndef _DOHPJACOBI_H
#define _DOHPJACOBI_H
/**
* @file   dohpjacobi.h
* @author Jed Brown <jed@59A2.org>
* @date   Fri Aug 22 19:40:10 2008
* 
* @brief  Interface to DohpJacobi, DohpBasis, and DohpRule.
*
* The purpose of these objects is to abstract the generation of basis functions.  DohpJacobi serves as a cache and
* lookup table for quadrature rules and Basis operations (derivatives and basis evalution).  The implementation may also
* provide unrolled versions for small to modest basis sizes.  In particular, the tightest loops (dofs/node, basis in
* minor dimension, nodes in minor dimension) can be arbitrarily unrolled.  This should reduce overhead of the \a hp
* method for regions of low order.
* 
*/

#include "private/petscimpl.h"

PETSC_EXTERN_CXX_BEGIN

/**
* This isn't actually an opaque structure.  Note that it is necessarily \c const since it will normally point to shared
* memory space.
* 
*/
typedef const struct _p_DohpRule * DohpRule;

/**
* Minimal struct to hold a 1-D quadrature rule.
*
* Although this does not support any functions, we still make it opaque.  This is so we can work with a DohpRule by
* passing pointers and it will look like any other Petsc/Dohp object.
* 
*/
struct _p_DohpRule {
  const PetscReal *node;              /**< nodal coordinates on the reference element [-1,1] */
  const PetscReal *weight;            /**< quadrature weights at the nodes  */
  PetscInt size;                /**< total number of points */
};


/**
* The handle for basis evaluation.
*
* The basis must always be compatible with a quadrature rule.  Eventually we'll have wrapper functions, but for now,
* just use the member functions.  Note that it is \c const since it will normally point to shared space.
* 
*/
typedef const struct _p_DohpBasis *DohpBasis;

/**
* Holds interpolation and differentiation context.
*
* This is implemented as pointers to the matrices and (possibly) unrolled versions of functions to apply those vectors
* across a 3D block of values.
* 
*/
struct _p_DohpBasis {
  PetscErrorCode (*basis)(DohpBasis,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscReal*,PetscReal*);
  PetscErrorCode (*basist)(DohpBasis,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscReal*,PetscReal*);
  PetscErrorCode (*deriv)(DohpBasis,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscReal*,PetscReal*);
  PetscErrorCode (*derivt)(DohpBasis,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscReal*,PetscReal*);
  const PetscReal *bmat;        /**< pointer to the basis matrix */
  const PetscReal *dmat;        /**< pointer to the derivative matrix */
  DohpRule rule;                /**< the rule on which the basis has been constructed */
  PetscInt bsize;               /**< basis size, should be redundant */
  PetscInt qsize;               /**< quadrature size, should be redundant */
};

/**
* Handle for setting up #DohpRule and #DohpBasis.
* 
*/
typedef struct _p_DohpJacobi *DohpJacobi;


/**
* Generic operations provided by a #DohpJacobi implementation.
* 
*/
struct _DohpJacobiOps {
  PetscErrorCode (*setup)(DohpJacobi);
  PetscErrorCode (*setfromoptions)(DohpJacobi);
  PetscErrorCode (*destroy)(DohpJacobi);
  PetscErrorCode (*view)(DohpJacobi,PetscViewer);
  PetscErrorCode (*getrule)(DohpJacobi,PetscInt,DohpRule*);
  PetscErrorCode (*getbasis)(DohpJacobi,PetscInt,PetscInt,DohpBasis*);
};

/**
* Private Jacobi table context.
* 
*/
struct _p_DohpJacobi {
  PETSCHEADER(struct _DohpJacobiOps);
  PetscInt basisdegree;         /**< the maximum degree basis functions to be supported */
  PetscInt ruleexcess;          /**< the amount of over-integration to be supported */
  PetscTruth setupcalled;
  void *impl;                   /**< private implementation context */
};

#define DohpJacobiType char *
#define DOHP_JACOBI_LGL "lgl"

EXTERN PetscErrorCode DohpJacobiCreate_LGL(DohpJacobi);

EXTERN PetscErrorCode DohpJacobiCreate(MPI_Comm,DohpJacobi*);
EXTERN PetscErrorCode DohpJacobiSetType(DohpJacobi,DohpJacobiType);
EXTERN PetscErrorCode DohpJacobiSetFromOptions(DohpJacobi);
EXTERN PetscErrorCode DohpJacobiSetUp(DohpJacobi);
EXTERN PetscErrorCode DohpJacobiSetFromOptions(DohpJacobi);
EXTERN PetscErrorCode DohpJacobiDestroy(DohpJacobi);
EXTERN PetscErrorCode DohpJacobiView(DohpJacobi,PetscViewer);
#define DohpJacobiRegisterDynamic(a,b,c,d) DohpJacobiRegister(a,b,c,d)
EXTERN PetscErrorCode DohpJacobiRegister(const char[],const char[],const char[],PetscErrorCode(*)(DohpJacobi));
EXTERN PetscErrorCode DohpJacobiRegisterAll(const char[]);
EXTERN PetscErrorCode DohpJacobiInitializePackage(const char[]);

EXTERN PetscErrorCode DohpJacobiSetDegrees(DohpJacobi,PetscInt,PetscInt);
EXTERN PetscErrorCode DohpJacobiGetRule(DohpJacobi,PetscInt,DohpRule*);
EXTERN PetscErrorCode DohpJacobiGetBasis(DohpJacobi,PetscInt,PetscInt,DohpBasis*);

PETSC_EXTERN_CXX_END

#endif /* _DOHPJACOBI_H */
