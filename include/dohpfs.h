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
#include "dohpunits.h"

dEXTERN_C_BEGIN

extern dClassId dFSROT_CLASSID;
extern PetscBool dFSRegisterAllCalled;

typedef struct _p_dFS *dFS;
typedef struct _p_dFSRotation *dFSRotation;
typedef struct _n_dRuleset *dRuleset;
typedef struct _n_dRulesetIterator *dRulesetIterator;

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
extern dErr dFSSetDegree(dFS,dJacobi,dMeshTag);
extern dErr dFSSetBlockSize(dFS,dInt);
extern dErr dFSGetBlockSize(dFS,dInt*);
extern dErr dFSSetFieldName(dFS,dInt,const char*);
extern dErr dFSGetFieldName(dFS,dInt,const char**);
extern dErr dFSSetFieldUnit(dFS fs,dInt fn,dUnit unit);
extern dErr dFSGetFieldUnit(dFS fs,dInt fn,dUnit *unit);
extern dErr dFSRegisterBoundary(dFS,dInt,dFSBStatus,dFSConstraintFunction,void*);
extern dErr dFSRegisterBoundarySet(dFS,dMeshESH,dFSBStatus,dFSConstraintFunction,void*);
extern dErr dFSSetFromOptions(dFS);
extern dErr dFSSetType(dFS,const dFSType);
extern dErr dFSSetOrderingType(dFS,const MatOrderingType);
extern dErr dFSCreateExpandedVector(dFS,Vec*);
extern dErr dFSCreateGlobalVector(dFS,Vec*);
extern dErr dFSExpandedToLocal(dFS,Vec,Vec,InsertMode);
extern dErr dFSLocalToExpanded(dFS,Vec,Vec,InsertMode);
extern dErr dFSInhomogeneousDirichletCommit(dFS fs,Vec gc);
extern dErr dFSDirichletProject(dFS fs,Vec X,dFSHomogeneousMode hmode);
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
extern dErr dFSGetSubElementMesh(dFS,dInt nelem,dInt nconn,dEntTopology topo[],dInt off[],dInt ind[]);
extern dErr dFSSubElementMeshView(dFS,dViewer);
extern dErr dFSGetDomain(dFS,dMeshESH*);

extern dErr dFSGetPreferredQuadratureRuleSet(dFS,dMeshESH,dEntType,dEntTopology,dQuadratureMethod,dRuleset*);
extern dErr dFSGetEFS(dFS,dRuleset,dInt*,const dEFS**);
extern dErr dFSRestoreEFS(dFS,dRuleset,dInt*,const dEFS**);
extern dErr dRulesetGetMaxQ(dRuleset rset,dInt *maxQ,dInt *maxnpatches,dInt *maxQelem);
extern dErr dRulesetGetSize(dRuleset rset,dInt *size);
extern dErr dRulesetGetWorkspace(dRuleset rset,dScalar **q,dScalar **cjac,dScalar **cjinv,dScalar **jw,dInt dof,...);
extern dErr dRulesetRestoreWorkspace(dRuleset rset,dScalar **q,dScalar **cjac,dScalar **cjinv,dScalar **jw,dInt dof,...);
extern dErr dRulesetDestroy(dRuleset*);

extern dErr dRulesetCreateIterator(dRuleset rset,dFS cfs,dRulesetIterator *iter);
extern dErr dRulesetIteratorDestroy(dRulesetIterator*);
extern dErr dRulesetIteratorSetMode(dRulesetIterator it,InsertMode imode);
extern dErr dRulesetIteratorAddFS(dRulesetIterator it,dFS fs);
extern dErr dRulesetIteratorStart(dRulesetIterator it,Vec X,...);
extern dErr dRulesetIteratorNextPatch(dRulesetIterator it);
extern bool dRulesetIteratorHasPatch(dRulesetIterator it);
extern dErr dRulesetIteratorSetupElement(dRulesetIterator it);
extern dErr dRulesetIteratorGetPatch(dRulesetIterator it,dRule *rule,dEFS *efs,dScalar **ex,dScalar **ey,...);
extern dErr dRulesetIteratorGetPatchSpace(dRulesetIterator it,dScalar **cjinv,dReal **jw,dScalar **u,dScalar **du,dScalar **v,dScalar **dv,...);
extern dErr dRulesetIteratorRestorePatchSpace(dRulesetIterator it,dScalar **cjinv,dReal **jw,dScalar **u,dScalar **du,dScalar **v,dScalar **dv,...);
extern dErr dRulesetIteratorCommitPatch(dRulesetIterator it,dScalar *v,...);
extern dErr dRulesetIteratorGetPatchApplied(dRulesetIterator it,dInt *Q,const dReal **jw,dScalar **u,dScalar **du,dScalar **v,dScalar **dv,...);
extern dErr dRulesetIteratorCommitPatchApplied(dRulesetIterator it,InsertMode imode,const dScalar *v,const dScalar *dv,...);
extern dErr dRulesetIteratorGetPatchExplicit(dRulesetIterator it,dInt *P,const dReal **interp,const dReal **deriv,...);
extern dErr dRulesetIteratorGetPatchAssembly(dRulesetIterator it,dInt *P,const dInt **rowcol,const dReal **interp,const dReal **deriv,...);
extern dErr dRulesetIteratorRestorePatchAssembly(dRulesetIterator it,dInt *P,const dInt **rowcol,const dReal **interp,const dReal **deriv,...);
extern dErr dRulesetIteratorGetElement(dRulesetIterator it,dRule *rule,dEFS *efs,dScalar **ex,dScalar **ey,...);
extern dErr dRulesetIteratorGetMatrixSpaceSplit(dRulesetIterator it,dScalar **K,...);
extern dErr dRulesetIteratorGetMatrixSpaceSizes(dRulesetIterator it,dInt *nrows,dInt *ncols,const dInt **sizes);
extern dErr dRulesetIteratorFinish(dRulesetIterator);
extern dErr dRulesetIteratorAddStash(dRulesetIterator it,dInt patchbytes,dInt nodebytes);
extern dErr dRulesetIteratorGetStash(dRulesetIterator,void *patchstash,void *nodestash);

extern dErr dFSDestroy(dFS*);
extern dErr dFSView(dFS,dViewer);
extern dErr dFSLoadIntoFS(dViewer,const char[],dFS);

#if defined PETSC_USE_DYNAMIC_LIBRARIES
#  define dFSRegisterDynamic(a,b,c,d) dFSRegister(a,b,c,0)
#else
#  define dFSRegisterDynamic(a,b,c,d) dFSRegister(a,b,c,d)
#endif

extern dErr dFSRegister(const char[],const char[],const char[],dErr(*)(dFS));
extern dErr dFSRegisterAll(const char[]);
extern dErr dFSInitializePackage(const char[]);
extern dErr dFSFinalizePackage(void);

/* These are purely for convenience */
extern dErr dFSGlobalToExpanded(dFS,Vec,Vec,dFSHomogeneousMode,InsertMode);
extern dErr dFSExpandedToGlobal(dFS,Vec,Vec,dFSHomogeneousMode,InsertMode);

extern dErr dFSRotationCreate(dFS,IS,dReal[],dInt[],Vec,dFSRotation*);
extern dErr dFSRotationDestroy(dFSRotation*);
extern dErr dFSRotationView(dFSRotation,dViewer);
extern dErr dFSRotationApply(dFSRotation,Vec,dFSRotateMode,dFSHomogeneousMode);
extern dErr dFSRotationApplyLocal(dFSRotation,Vec,dFSRotateMode,dFSHomogeneousMode);

extern dErr dFSSetRotation(dFS,dFSRotation);
extern dErr dFSGetRotation(dFS,dFSRotation*);

extern dErr dFSGetCoordinateFS(dFS,dFS*);
extern dErr dFSGetGeometryVectorExpanded(dFS,Vec*);
extern dErr dFSGetBoundingBox(dFS fs,dReal bbox[3][2]);

extern dErr dFSGetNodalCoordinateFS(dFS,dFS*);
extern dErr dFSGetNodalCoordinatesExpanded(dFS,Vec*);
extern dErr dFSGetNodalCoordinatesGlobal(dFS,Vec*);

extern dErr dFSRedimension(dFS,dInt,dFSClosureMode,dFS*);
extern dErr VecDohpGetFS(Vec,dFS*);

/* The normal MatAXPY() does not work for arguments of mixed type, so we redeclare this "private" function here. */
extern PetscErrorCode MatAXPY_Basic(Mat Y,PetscScalar a,Mat X,MatStructure str);

dEXTERN_C_END
#endif  /* _DOHPFS_H */
