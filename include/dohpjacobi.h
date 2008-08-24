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
* The idea is that when building a #dMFS we would like the maximum possible sharing and flexibility.  The #dJacobi
* assembles a #dRule and #dEFS for use by any part of the #dMFS.  All the real data (like interpolation and
* differentiation matrices) are shared.
* 
*/

#include "dohptype.h"
// #include "private/petscimpl.h"

PETSC_EXTERN_CXX_BEGIN

/**
* Handle for manipulating EFS objects.  Since these are actually stored in arrays, the handle is the actual object
* rather than just a pointer.
* 
*/
typedef struct m_dEFS dEFS;

/**
* As above, the handle is the actual object.
* 
*/
typedef struct m_dRule dRule;


/**
* Handle for setting up #dRule and #dEFS contexts.
* 
*/
typedef struct p_dJacobi *dJacobi;

#define dJacobiType char *
#define dJACOBI_TENSOR "tensor"

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
EXTERN dErr dJacobiGetRule(dJacobi jac,dTopology top,const dInt rsize[],dInt left,dRule *rule,dInt *bytes);
EXTERN dErr dJacobiGetEFS(dJacobi jac,dTopology top,const dInt bsize[],const dRule *rule,dInt left,dEFS *efs,dInt *bytes);

  // EXTERN dErr dJacobiGetRule(dJacobi jac,dTopology top,const dInt rsize[],dInt left,dRule *rule,dInt *bytes);
// EXTERN dErr dJacobiGetEFS(dJacobi jac,dTopology top,const dInt bsize[],const dRule *rule,dInt left,dEFS *efs,dInt *bytes);

PETSC_EXTERN_CXX_END

#endif /* _DOHPJACOBI_H */
