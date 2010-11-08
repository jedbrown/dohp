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

dEXTERN_C_BEGIN

extern dClassId dFSROT_CLASSID;

typedef struct _p_dFS *dFS;
typedef struct _p_dFSRotation *dFSRotation;
typedef struct _n_dRuleset *dRuleset;

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
typedef enum {dFS_CLOSURE, dFS_INTERIOR} dFSClosureMode;
typedef enum {dFS_ROTATE_FORWARD, dFS_ROTATE_REVERSE} dFSRotateMode;

extern const char *const dFSHomogeneousModes[];
extern const char *const dFSClosureModes[];
extern const char *const dFSRotateModes[];

#define dFSType char *

#define dFSCONT "cont"

extern dErr dFSCreate(MPI_Comm,dFS*);
extern dErr dFSSetMesh(dFS,dMesh,dMeshESH); /* mesh, active set */
extern dErr dFSGetMesh(dFS,dMesh*);
extern dErr dFSGetJacobi(dFS,dJacobi*);
extern dErr dFSSetRuleTag(dFS,dJacobi,dMeshTag);
extern dErr dFSSetDegree(dFS,dJacobi,dMeshTag);
extern dErr dFSSetBlockSize(dFS,dInt);
extern dErr dFSGetBlockSize(dFS,dInt*);
extern dErr dFSSetFieldName(dFS,dInt,const char*);
extern dErr dFSRegisterBoundary(dFS,dInt,dFSBStatus,dFSConstraintFunction,void*);
extern dErr dFSRegisterBoundarySet(dFS,dMeshESH,dFSBStatus,dFSConstraintFunction,void*);
extern dErr dFSSetFromOptions(dFS);
extern dErr dFSSetType(dFS,const dFSType);
extern dErr dFSCreateExpandedVector(dFS,Vec*);
extern dErr dFSCreateGlobalVector(dFS,Vec*);
extern dErr dFSExpandedToLocal(dFS,Vec,Vec,InsertMode);
extern dErr dFSLocalToExpanded(dFS,Vec,Vec,InsertMode);
extern dErr dFSInhomogeneousDirichletCommit(dFS fs,Vec gc);
extern dErr dFSRotateGlobal(dFS,Vec,dFSRotateMode,dFSHomogeneousMode);
#if 0
extern dErr dFSGetElements(dFS,dInt*,dInt*restrict*,dRule*restrict*,dEFS*restrict*,dInt*restrict*,dReal(*restrict*)[3]);
extern dErr dFSRestoreElements(dFS,dInt*,dInt*restrict*,dRule*restrict*,dEFS*restrict*,dInt*restrict*,dReal(*restrict*)[3]);
#endif
extern dErr dFSMatSetValuesBlockedExpanded(dFS,Mat,dInt,const dInt[],dInt,const dInt[],const dScalar[],InsertMode);
extern dErr dFSGetMatrix(dFS,const MatType,Mat*);
extern dErr dFSSetOptionsPrefix(dFS,const char[]);
extern dErr dFSBuildSpace(dFS);
extern dErr dFSGetSubElementMeshSize(dFS,dInt*,dInt*,dInt *);
extern dErr dFSGetSubElementMesh(dFS,dInt nelem,dInt nvtx,dEntTopology topo[],dInt off[],dInt ind[]);
extern dErr dFSGetDomain(dFS,dMeshESH*);

extern dErr dFSGetPreferredQuadratureRuleSet(dFS,dMeshESH,dEntType,dEntTopology,dQuadratureMethod,dRuleset*);
extern dErr dFSGetEFS(dFS,dRuleset,dInt*,const dEFS**);
extern dErr dFSRestoreEFS(dFS,dRuleset,dInt*,const dEFS**);
extern dErr dRulesetGetWorkspace(dRuleset rset,dScalar **q,dScalar **cjac,dScalar **cjinv,dScalar **jw,dInt dof,...);
extern dErr dRulesetRestoreWorkspace(dRuleset rset,dScalar **q,dScalar **cjac,dScalar **cjinv,dScalar **jw,dInt dof,...);
extern dErr dRulesetDestroy(dRuleset);

extern dErr dFSDestroy(dFS);
extern dErr dFSView(dFS,dViewer);
extern dErr dFSLoadIntoFS(dViewer,const char[],dFS);
#define dFSRegisterDynamic(a,b,c,d) dFSRegister(a,b,c,d)
extern dErr dFSRegister(const char[],const char[],const char[],dErr(*)(dFS));
extern dErr dFSRegisterAll(const char[]);
extern dErr dFSInitializePackage(const char[]);

/* These are purely for convenience */
extern dErr dFSGlobalToExpanded(dFS,Vec,Vec,dFSHomogeneousMode,InsertMode);
extern dErr dFSExpandedToGlobal(dFS,Vec,Vec,dFSHomogeneousMode,InsertMode);

extern dErr dFSRotationCreate(dFS,IS,dReal[],dInt[],Vec,dFSRotation*);
extern dErr dFSRotationDestroy(dFSRotation);
extern dErr dFSRotationView(dFSRotation,dViewer);
extern dErr dFSRotationApply(dFSRotation,Vec,dFSRotateMode,dFSHomogeneousMode);
extern dErr dFSRotationApplyLocal(dFSRotation,Vec,dFSRotateMode,dFSHomogeneousMode);

extern dErr dFSSetRotation(dFS,dFSRotation);
extern dErr dFSGetRotation(dFS,dFSRotation*);

extern dErr dFSGetCoordinateFS(dFS,dFS*);
extern dErr dFSGetGeometryVectorExpanded(dFS,Vec*);
extern dErr dFSGetGeometryVectorGlobal(dFS,Vec*);

extern dErr dFSGetNodalCoordinateFS(dFS,dFS*);
extern dErr dFSGetNodalCoordinatesExpanded(dFS,Vec*);
extern dErr dFSGetNodalCoordinatesGlobal(dFS,Vec*);

extern dErr dFSRedimension(dFS,dInt,dFSClosureMode,dFS*);

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

#define dQ1SCALE_DECLARE(tscale,scale,i,j,k)                            \
  const dReal dUNUSED (scale)[8] = {(tscale)[0][(i)+0]*(tscale)[1][(j)+0]*(tscale)[2][(k)+0], \
                                    (tscale)[0][(i)+1]*(tscale)[1][(j)+0]*(tscale)[2][(k)+0], \
                                    (tscale)[0][(i)+1]*(tscale)[1][(j)+1]*(tscale)[2][(k)+0], \
                                    (tscale)[0][(i)+0]*(tscale)[1][(j)+1]*(tscale)[2][(k)+0], \
                                    (tscale)[0][(i)+0]*(tscale)[1][(j)+0]*(tscale)[2][(k)+1], \
                                    (tscale)[0][(i)+1]*(tscale)[1][(j)+0]*(tscale)[2][(k)+1], \
                                    (tscale)[0][(i)+1]*(tscale)[1][(j)+1]*(tscale)[2][(k)+1], \
                                    (tscale)[0][(i)+0]*(tscale)[1][(j)+1]*(tscale)[2][(k)+1]}

extern dErr dQ1HexComputeQuadrature(const dReal x[8][3],dInt *n,const dReal (**qx)[3],const dReal **jw,const dReal **basis,const dReal **deriv);


dEXTERN_C_END
#endif  /* _DOHPFS_H */
