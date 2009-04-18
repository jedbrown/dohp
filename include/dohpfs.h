#ifndef _DOHPFS_H
#define _DOHPFS_H
/**
* @file   dohpfs.h
* @author Jed Brown <jed@59A2.org>
* @date   Sun Sep  7 17:46:54 2008
*
* @brief  The function space object
*
*
*/

#include "dohpmesh.h"
#include "dohpjacobi.h"
#include "dohpquotient.h"

PETSC_EXTERN_CXX_BEGIN

typedef struct _p_dFS *dFS;

/** User-provided constraint function.
* @param ctx User context
* @param x coordinates of node to generate constraints for (vector of length 3)
* @param b basis at node (3 vectors of length 3), the first vector \c b[0] is the normal vector, the others (\c b[1] and \c b[2]) are tangent vectors.
* @param T constraint matrix, the first \p g vectors correspond to global basis functions, the remaining vectors represent Dirichlet functions
* @param g Number of values to be represented in global system
*
* @note This function \b must be pure (no side effects, only output is T, g, and possible error code).
*
* @note The number of values in the global system \p g should \b not depend on the location or the normal (i.e. it is
*   constant for each boundary type).  Rationale: setting the number of global dofs separately makes it possible for the
*   constraint function to become out of sync with the number of global dofs.
**/
typedef dErr (*dFSConstraintFunction)(void *ctx,const dReal x[],const dReal b[],dReal T[],dInt *g);

enum { dFSBSTATUS_MASK = 0xffff }; /* Reserved for number of strongly enforced dofs */
typedef enum {
  dFSBSTATUS_INTERIOR  = 0,
  dFSBSTATUS_DIRICHLET = (dFSBSTATUS_MASK+1) << 0, /* Do not represent in global space */
  dFSBSTATUS_WEAK      = (dFSBSTATUS_MASK+1) << 1, /* A weak form must be evaluated on this element */
} dFSBStatus;

typedef enum {dFS_HOMOGENEOUS, dFS_INHOMOGENEOUS} dFSHomogeneousMode;

#define dFSType char *

#define dFSCONT "cont"

EXTERN dErr dFSCreate(MPI_Comm,dFS*);
EXTERN dErr dFSSetMesh(dFS,dMesh,dMeshESH); /* mesh, active set */
EXTERN dErr dFSSetRuleTag(dFS,dJacobi,dMeshTag);
EXTERN dErr dFSSetDegree(dFS,dJacobi,dMeshTag);
EXTERN dErr dFSSetBlockSize(dFS,dInt);
EXTERN dErr dFSRegisterBoundary(dFS,dInt,dFSBStatus,dFSConstraintFunction,void*);
EXTERN dErr dFSSetFromOptions(dFS);
EXTERN dErr dFSSetType(dFS,const dFSType);
EXTERN dErr dFSCreateExpandedVector(dFS,Vec*);
EXTERN dErr dFSCreateGlobalVector(dFS,Vec*);
EXTERN dErr dFSCreateDirichletVector(dFS,Vec*);
EXTERN dErr dFSGlobalToExpandedBegin(dFS,Vec,dFSHomogeneousMode,Vec);
EXTERN dErr dFSGlobalToExpandedEnd(dFS,Vec,dFSHomogeneousMode,Vec);
EXTERN dErr dFSExpandedToGlobal(dFS,Vec,InsertMode,Vec);
EXTERN dErr dFSExpandedToGlobalBegin(dFS,Vec,InsertMode,Vec);
EXTERN dErr dFSExpandedToGlobalEnd(dFS,Vec,InsertMode,Vec);
EXTERN dErr dFSExpandedToDirichlet(dFS,Vec,InsertMode,Vec);
EXTERN dErr dFSRotateLocalVector(dFS,Vec xloc,Vec save_strong,Vec correct_strong);
EXTERN dErr dFSUnRotateLocalVector(dFS,Vec,Vec,Vec);
EXTERN dErr dFSGetElements(dFS,dInt*,dInt*restrict*,s_dRule*restrict*,s_dEFS*restrict*,dInt*restrict*,dReal(*restrict*)[3]);
EXTERN dErr dFSRestoreElements(dFS,dInt*,dInt*restrict*,s_dRule*restrict*,s_dEFS*restrict*,dInt*restrict*,dReal(*restrict*)[3]);
EXTERN dErr dFSGetWorkspace(dFS,const char[],dReal(*restrict*)[3],dReal(*restrict*)[3][3],dReal*restrict*,dScalar*restrict*,dScalar*restrict*,dScalar*restrict*,dScalar*restrict*);
EXTERN dErr dFSRestoreWorkspace(dFS,const char[],dReal(*restrict*)[3],dReal(*restrict*)[3][3],dReal*restrict*,dScalar*restrict*,dScalar*restrict*,dScalar*restrict*,dScalar*restrict*);
EXTERN dErr dFSMatSetValuesBlockedExpanded(dFS,Mat,dInt,const dInt[],dInt,const dInt[],const dScalar[],InsertMode);
EXTERN dErr dFSGetMatrix(dFS,const MatType,Mat*);
EXTERN dErr dFSBuildSpace(dFS);

EXTERN dErr dFSDestroy(dFS);
EXTERN dErr dFSView(dFS,PetscViewer);
#define dFSRegisterDynamic(a,b,c,d) dFSRegister(a,b,c,d)
EXTERN dErr dFSRegister(const char[],const char[],const char[],dErr(*)(dFS));
EXTERN dErr dFSRegisterAll(const char[]);
EXTERN dErr dFSInitializePackage(const char[]);



/** The Q1 stuff doesn't really belong here, but it is used at the same level of abstraction and I'm too lazy to
* separate it out yet.  This macro is a rather dirty way to avoid a bunch of code duplication without much runtime cost.
* When there's time, I should make this more elegant.
*
* @param c An const array is declared with this name, it gives the element-wise index for each corner (basis function)
* on the Q1 sub-element.  The indexing is determined using \a p (3-array of number of basis functions in each Cartesian
* direction) and \a i, \a j, \a k, the indices in the tensor element.
*
* @param rowcol A const array with this name is declared, it is just \a c shifted by \a off and is normally only used to
* provide indices for dFSMatSetValuesExpanded.
*
* @param corners A const array \c dReal[8][3] with this name is declared which holds the global coordinates of each corner.
* In many cases, it will simply be passed to dQ1HexComputeQuadrature.
*
* @param off Offset of this element in expanded vector.
* @param nx Coordinates of all the nodes on the element, these will be used to create the \a corners array
* @param p Array \c dInt[3] with the number of Lagrange nodes in each Cartesian direction
* @param i Index in element basis x-direction
* @param j Index in element basis y-direction
* @param k Index in element basis z-direction
**/
#define dQ1CORNER_CONST_DECLARE(c,rowcol,corners,off,nx,p,i,j,k)        \
          const dInt (c)[8] = {((((i)+0)*(p)[1]+(j)+0)*(p)[2]+(k)),     \
                               ((((i)+1)*(p)[1]+(j)+0)*(p)[2]+(k)),     \
                               ((((i)+1)*(p)[1]+(j)+1)*(p)[2]+(k)),     \
                               ((((i)+0)*(p)[1]+(j)+1)*(p)[2]+(k)),     \
                               ((((i)+0)*(p)[1]+(j)+0)*(p)[2]+(k)+1),   \
                               ((((i)+1)*(p)[1]+(j)+0)*(p)[2]+(k)+1),   \
                               ((((i)+1)*(p)[1]+(j)+1)*(p)[2]+(k)+1),   \
                               ((((i)+0)*(p)[1]+(j)+1)*(p)[2]+(k)+1)};  \
          const dInt (rowcol)[8] = {(off)+(c)[0],(off)+(c)[1],(off)+(c)[2],(off)+(c)[3], \
                                    (off)+(c)[4],(off)+(c)[5],(off)+(c)[6],(off)+(c)[7]}; \
          const dReal (corners)[8][3] = {{(nx)[(c)[0]][0],(nx)[(c)[0]][1],(nx)[(c)[0]][2]}, \
                                         {(nx)[(c)[1]][0],(nx)[(c)[1]][1],(nx)[(c)[1]][2]}, \
                                         {(nx)[(c)[2]][0],(nx)[(c)[2]][1],(nx)[(c)[2]][2]}, \
                                         {(nx)[(c)[3]][0],(nx)[(c)[3]][1],(nx)[(c)[3]][2]}, \
                                         {(nx)[(c)[4]][0],(nx)[(c)[4]][1],(nx)[(c)[4]][2]}, \
                                         {(nx)[(c)[5]][0],(nx)[(c)[5]][1],(nx)[(c)[5]][2]}, \
                                         {(nx)[(c)[6]][0],(nx)[(c)[6]][1],(nx)[(c)[6]][2]}, \
                                         {(nx)[(c)[7]][0],(nx)[(c)[7]][1],(nx)[(c)[7]][2]}}

EXTERN dErr dQ1HexComputeQuadrature(const dReal x[8][3],dInt *n,const dReal (**qx)[3],const dReal **jw,const dReal **basis,const dReal **deriv);


PETSC_EXTERN_CXX_END
#endif  /* _DOHPFS_H */
