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
#include "petscsys.h"
// #include "private/petscimpl.h"

PETSC_EXTERN_CXX_BEGIN

/**
* Handle for manipulating EFS objects.  Since these are actually stored in arrays, the handle is the actual object
* rather than just a pointer.  This means that a certain amount of library code (dFS) will have to be able to see the
* struct definition (in order to know the size, to make an array).  This is okay since the objects are so simple
* (everything is implementation-dependent, the details of which are still hidden).
* 
*/
typedef struct p_dEFS dEFS;

/**
* As above, the handle is the actual object.
* 
*/
typedef struct p_dRule dRule;

/**
* Indicates whether or not to apply the transpose of a interpolation/derivative matrix.
* 
*/
typedef enum { dTRANSPOSE_NO=113634,dTRANSPOSE_YES=853467 } dTransposeMode;

/**
* Indicates what type of basis operation to do, see dEFSApply().
* 
*/
typedef enum { dAPPLY_INTERP=40001,dAPPLY_INTERP_TRANSPOSE=40002,dAPPLY_GRAD=40003,dAPPLY_GRAD_TRANSPOSE=40004 } dApplyMode;

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
EXTERN dErr dJacobiGetRule(dJacobi jac,dTopology top,const dInt rsize[],dRule *rule,void **base,dInt *index);
EXTERN dErr dJacobiGetEFS(dJacobi jac,dTopology top,const dInt bsize[],dRule *rule,dEFS *efs,void **base,dInt *index);

PETSC_EXTERN_CXX_END

#endif /* _DOHPJACOBI_H */
