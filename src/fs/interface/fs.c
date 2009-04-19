#include "private/fsimpl.h"

dErr dFSSetMesh(dFS fs,dMesh mesh,dMeshESH active)
{
  dMesh qmesh;
  dErr  err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(mesh,dMESH_COOKIE,2);
  if (fs->quotient) {
    err = dQuotientGetMesh(fs->quotient,&qmesh);dCHK(err);
    if (mesh != qmesh) fs->quotient = 0; /* The Quotient is stale */
  }
  fs->mesh = mesh;
  fs->activeSet = active;
  err = dMeshGetTag(mesh,fs->bdyTagName,&fs->bdyTag);dCHK(err);
  err = dMeshTagCreateTemp(mesh,"boundary_status",1,dDATA_INT,&fs->bstatusTag);dCHK(err);
  err = dMeshTagCreateTemp(mesh,"boundary_constraint",sizeof(struct dFSConstraintCtx),dDATA_BYTE,&fs->bdyConstraintTag);dCHK(err);
  err = dMeshTagCreateTemp(mesh,"global_offset",1,dDATA_INT,&fs->goffsetTag);dCHK(err);
  err = dMeshTagCreateTemp(mesh,"global_dirichlet_offset",1,dDATA_INT,&fs->gdoffsetTag);dCHK(err);
  err = dMeshTagCreateTemp(mesh,"local_offset",1,dDATA_INT,&fs->loffsetTag);dCHK(err);
  err = dMeshTagCreate(mesh,"global_closure_offset",1,dDATA_INT,&fs->gcoffsetTag);dCHK(err);
  err = dMeshSetCreate(mesh,&fs->ownedExplicitSet);dCHK(err);
  err = dMeshSetCreate(mesh,&fs->ghostExplicitSet);dCHK(err);
  err = dMeshSetCreate(mesh,&fs->ownedDirichletSet);dCHK(err);
  err = dMeshSetCreate(mesh,&fs->ghostDirichletSet);dCHK(err);
  err = dMeshSetCreate(mesh,&fs->weakFaceSet);dCHK(err);
  err = dMeshSetCreate(mesh,&fs->boundaries);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSSetRuleTag(dFS fs,dJacobi jac,dMeshTag rtag)
{

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  fs->ruletag = rtag;
  if (jac && fs->jacobi && fs->jacobi != jac) dERROR(1,"cannot change dJacobi");
  if (jac) fs->jacobi = jac;
  dFunctionReturn(0);
}

dErr dFSSetDegree(dFS fs,dJacobi jac,dMeshTag deg)
{

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(jac,dJACOBI_COOKIE,2);
  fs->degreetag = deg;
  if (jac && fs->jacobi && fs->jacobi != jac) dERROR(1,"cannot change dJacobi");
  if (jac) fs->jacobi = jac;
  dFunctionReturn(0);
}

/** Set the block size for a function space.
*
* @param fs function space
* @param bs block size (number of dofs per node)
**/
dErr dFSSetBlockSize(dFS fs,dInt bs)
{

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  fs->bs = bs;
  dFunctionReturn(0);
}

/** Register a boundary condition with the function space.
* After all boundary conditions are registered, dFSBuildSpace (called by dFSSetFromOptions) can be used.
*
* @param fs function space object
* @param mid Boundary ID, usually the value of the NEUMANN_SET tag
* @param bstat Boundary status, determines how boundary should be represented (weak, global, dirichlet)
* @param cfunc constraint function, only makes sense if !strong, but some degrees of freedom must be removed anyway
* @param user context for constraint function
*
* @note Collective on \p fs
*
* @note The constraint function \b must be a pure function (no side-effects, only writes to it's output matrix) with the
* same definition on every process.  The constraint matrix \b must be invertible and currently must be orthogonal.
* Support for general constraint matrices is easy, but of doubtful usefulness.  The number of dofs declared global and
* local should be the same at every point (this is not actually essential, but it's convenient).  The reason it is not
* declared statically (outside of the function definition) is merely to avoid duplicating information that must be kept
* consistent.
**/
dErr dFSRegisterBoundary(dFS fs,dInt mid,dFSBStatus bstat,dFSConstraintFunction cfunc,void *user)
{
  dMeshESH         bset;
  iMesh_Instance   mi;
  dErr             err;
  dIInt            ierr;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  err = dFSBStatusValid(bstat);dCHK(err);
  if (dFSBStatusStrongCount(bstat) > fs->bs) dERROR(1,"Cannot impose strong conditions on more dofs than the block size");
  err = dMeshGetTaggedSet(fs->mesh,fs->bdyTag,&mid,&bset);dCHK(err);
  err = dMeshTagSSetData(fs->mesh,fs->bstatusTag,&bset,1,&bstat,1,dDATA_INT);dCHK(err);
  if (cfunc) {
    struct dFSConstraintCtx ctx;
    ctx.cfunc = cfunc;
    ctx.user = user;
    err = dMeshTagSSetData(fs->mesh,fs->bdyConstraintTag,&bset,1,&ctx,sizeof(ctx),dDATA_BYTE);dCHK(err);
  }
  err = dMeshGetInstance(fs->mesh,&mi);dCHK(err);
  iMesh_addEntSet(mi,bset,fs->boundaries,&ierr);dICHK(mi,ierr);
  dFunctionReturn(0);
}

dErr dFSView(dFS fs,dViewer viewer)
{
  dBool iascii;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  if (!viewer) {
    err = PetscViewerASCIIGetStdout(((dObject)fs)->comm,&viewer);dCHK(err);
  }
  dValidHeader(viewer,PETSC_VIEWER_COOKIE,2);
  PetscCheckSameComm(fs,1,viewer,2);

  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);dCHK(err);
  if (iascii) {
    err = PetscViewerASCIIPrintf(viewer,"dFS object:(%s)\n",
                                  ((dObject)fs)->prefix ? ((dObject)fs)->prefix : "no prefix");dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"type: %s\n",
                                  ((dObject)fs)->type_name ? ((dObject)fs)->type_name : "type not set");dCHK(err);
    if (!fs->spacebuilt) {
      err = PetscViewerASCIIPrintf(viewer,"Function Space has not been built.\n");dCHK(err);
    }
    {
      dInt nents[4];
      err = PetscViewerASCIIPrintf(viewer,"General information about the mesh topology.\n");dCHK(err);
      for (dEntType type=dTYPE_VERTEX; type<dTYPE_ALL; type++) {
        err = dMeshGetNumEnts(fs->mesh,fs->activeSet,type,dTOPO_ALL,&nents[type]);dCHK(err);
      }
      err = PetscViewerASCIIPrintf(viewer,"number of vertices=%d edges=%d faces=%d regions=%d\n",nents[0],nents[1],nents[2],nents[3]);dCHK(err);
    }
    {                           /* print aggregate sizes */
      PetscMPIInt gm[2],lm[2];
      err = MatGetSize(fs->E,&lm[0],&lm[1]);dCHK(err);
      err = MPI_Reduce(lm,gm,2,MPI_INT,MPI_SUM,0,((dObject)fs)->comm);dCHK(err);
      err = PetscViewerASCIIPrintf(viewer,"on rank 0, %d/%d element dofs constrained against %d/%d local dofs\n",
                                   lm[0],gm[0],lm[1],gm[1]);dCHK(err);
    }
    if (fs->ops->view) {
      err = (*fs->ops->view)(fs,viewer);dCHK(err);
    } else {
      err = PetscViewerASCIIPrintf(viewer,"Internal info not available.\n");dCHK(err);
    }
    err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  } else if (fs->ops->view) {
    err = (*fs->ops->view)(fs,viewer);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dFSDestroy(dFS fs)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  if (fs->ops->impldestroy) {
    err = (*fs->ops->impldestroy)(fs);dCHK(err);
  }
  for (dInt i=0; i<dFS_MAX_WORKSPACES; i++) {
    s_dFSWorkspace *w = &fs->workspace[i];
    switch (w->status) {
      case 0: continue;         /* Nothing here, nothing allocated */
      case 1:                   /* Allocated, but is not currently checked out */
        err = dFree7(w->q,w->jinv,w->jw,w->u,w->v,w->du,w->dv);dCHK(err);
        break;
      case 2:                   /* Still checked out */
        dERROR(1,"Workspace checked out under name `%s' not restored",w->name);
      default: dERROR(1,"Invalid status %d",w->status);
    }
  }
  if (fs->slice) {err = SlicedDestroy(fs->slice);dCHK(err);}
  if (fs->dslice) {err = SlicedDestroy(fs->dslice);dCHK(err);}
  err = VecGhostRestoreLocalForm(fs->d,&fs->dl);dCHK(err);
  err = VecDestroy(fs->d);dCHK(err);
  err = VecDestroy(fs->gc);dCHK(err);
  err = VecScatterDestroy(fs->ctod);dCHK(err);
  err = VecScatterDestroy(fs->ctog);dCHK(err);
  err = MatDestroy(fs->E);dCHK(err);
  err = MatDestroy(fs->Ep);dCHK(err);
  err = MatDestroy(fs->Ed);dCHK(err);
  err = dFree3(fs->rule,fs->efs,fs->off);dCHK(err);
  err = dMeshRestoreVertexCoords(fs->mesh,fs->nelem,NULL,&fs->vtxoff,&fs->vtx);dCHK(err);
  err = PetscHeaderDestroy(fs);dCHK(err);
  dFunctionReturn(0);
}

/**
* Builds a function space.  Enforcement of constraints is implementation dependent.
*
* @param fs the function space
*
* @return err
*/
dErr dFSBuildSpace(dFS fs)
{
  Vec x,g;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  if (fs->spacebuilt) dERROR(1,"The space is already built, rebuilding is not implemented");
  if (fs->ops->buildspace) {
    err = (*fs->ops->buildspace)(fs);dCHK(err);
  }

  /* Determine the number of elements in which each dof appears */
  err = dFSCreateExpandedVector(fs,&x);dCHK(err);
  err = dFSCreateGlobalVector(fs,&g);dCHK(err);
  err = VecSet(x,1);dCHK(err);
  err = VecZeroEntries(g);dCHK(err);
  err = dFSExpandedToGlobal(fs,x,ADD_VALUES,g);dCHK(err);
  err = VecGhostUpdateBegin(g,ADD_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecGhostUpdateEnd(g,ADD_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecDestroy(x);dCHK(err);

  /* \todo Use g to set sparsity pattern */
  err = VecDestroy(g);dCHK(err);

  fs->spacebuilt = dTRUE;
  dFunctionReturn(0);
}

dErr dFSCreateExpandedVector(dFS fs,Vec *x)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidPointer(x,2);
  err = MatGetVecs(fs->E,NULL,x);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSCreateGlobalVector(dFS fs,Vec *g)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidPointer(g,2);
  err = SlicedCreateGlobalVector(fs->slice,g);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSCreateDirichletVector(dFS fs,Vec *d)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidPointer(d,2);
  err = SlicedCreateGlobalVector(fs->dslice,d);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSGlobalToExpandedBegin(dFS dUNUSED fs,Vec g,dFSHomogeneousMode dUNUSED hmode,Vec dUNUSED x)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(g,VEC_COOKIE,2);
  dValidHeader(x,VEC_COOKIE,4);
  err = VecGhostUpdateBegin(g,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSGlobalToExpandedEnd(dFS fs,Vec g,dFSHomogeneousMode hmode,Vec x)
{
  Vec lform;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(g,VEC_COOKIE,2);
  dValidHeader(x,VEC_COOKIE,4);
  switch (hmode) {
    case dFS_HOMOGENEOUS:
      err = VecZeroEntries(x);dCHK(err);
      break;
    case dFS_INHOMOGENEOUS:
      err = MatMult(fs->Ed,fs->dl,x);dCHK(err);
      break;
    default:
      dERROR(1,"dFSHomogeneousMode %d not supported",hmode);
  }
  err = VecGhostUpdateEnd(g,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecGhostGetLocalForm(g,&lform);dCHK(err);
  err = MatMultAdd(fs->E,lform,x,x);dCHK(err);
  err = VecGhostRestoreLocalForm(g,&lform);dCHK(err);
  dFunctionReturn(0);
}

/** Assemble expanded vector into global form.
*
* @param fs Function space
* @param x  Expanded vector
* @param imode sum into local form or overwrite it, see note below
* @param g  Global vector (must be a VecGhost, as obtained with dFSCreateGlobalVector()), see note 2 below
*
* @note \a imode is treated differently here than in most of PETSc, it never runs a scatter with INSERT_VALUES since
* there are \e always multiple inputs mapped to the same output (at least for continuous spaces).  Instead \a imode
* determines whether to clear the local form before assembly (INSERT_VALUES) or whether to just sum into the existing
* vector (ADD_VALUES).
*
* @note Not collective!  While the name indicates that it only updates the global vector, it also updates the local
* form.  For the collective operation which sums contributions from multiple processes, see dFSExpandedToGlobalBegin().
*
**/
dErr dFSExpandedToGlobal(dFS fs,Vec x,InsertMode imode,Vec g)
{
  Vec lform;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(x,VEC_COOKIE,2);
  dValidHeader(g,VEC_COOKIE,4);
  err = VecGhostGetLocalForm(g,&lform);dCHK(err);
  switch (imode) {
    case INSERT_VALUES:
      err = MatMultTranspose(fs->E,x,lform);dCHK(err);
      break;
    case ADD_VALUES:
      err = MatMultTransposeAdd(fs->E,x,lform,lform);dCHK(err);
      break;
    default:
      dERROR(1,"InsertMode %d not supported",imode);
  }
  err = VecGhostRestoreLocalForm(g,&lform);dCHK(err);
  dFunctionReturn(0);
}

/** Updates the global values, does \b not broadcast the global values back to the ghosts.
*
* Call VecGhostUpdate{Begin,End}(g,INSERT_VALUES,SCATTER_FORWARD) after this to update ghost values.
**/
dErr dFSExpandedToGlobalBegin(dFS fs,Vec x,InsertMode imode,Vec g)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(x,VEC_COOKIE,2);
  dValidHeader(g,VEC_COOKIE,4);
  err = dFSExpandedToGlobal(fs,x,imode,g);dCHK(err);
  err = VecGhostUpdateBegin(g,ADD_VALUES,SCATTER_REVERSE);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSExpandedToGlobalEnd(dFS dUNUSED fs,Vec dUNUSED x,InsertMode dUNUSED imode,Vec g)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(x,VEC_COOKIE,2);
  dValidHeader(g,VEC_COOKIE,4);
  err = VecGhostUpdateEnd(g,ADD_VALUES,SCATTER_REVERSE);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSExpandedToDirichlet(dFS fs,Vec x,InsertMode imode,Vec d)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(d,VEC_COOKIE,2);
  dValidHeader(x,VEC_COOKIE,4);
  switch (imode) {
    case INSERT_VALUES:
      err = MatMultTranspose(fs->Ed,x,d);dCHK(err);
      break;
    case ADD_VALUES:
      err = MatMultTransposeAdd(fs->Ed,x,d,d);dCHK(err);
      break;
    default:
      dERROR(1,"InsertMode %d not supported",imode);
  }
  dFunctionReturn(0);
}

dErr dFSClosureToGlobal(dFS fs,Vec c,Vec g,InsertMode imode,dFSHomogeneousMode hmode)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(c,VEC_COOKIE,2);
  dValidHeader(g,VEC_COOKIE,3);
  err = VecScatterBegin(fs->ctog,c,g,imode,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd(fs->ctog,c,g,imode,SCATTER_FORWARD);dCHK(err);
  err = dFSRotationApply(fs->rot,g,hmode,dFS_ROTATE_FORWARD);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSGlobalToClosure(dFS fs,Vec g,Vec c,InsertMode imode,dFSHomogeneousMode dUNUSED hmode)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(g,VEC_COOKIE,2);
  dValidHeader(c,VEC_COOKIE,3);
  err = dFSRotationApply(fs->rot,g,hmode,dFS_ROTATE_REVERSE);dCHK(err);
  err = VecScatterBegin(fs->ctog,g,c,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd(fs->ctog,g,c,imode,SCATTER_REVERSE);dCHK(err);
  err = dFSRotationApply(fs->rot,g,dFS_INHOMOGENEOUS,dFS_ROTATE_FORWARD);dCHK(err);
  /* \bug If the input vector did not satisfy the strongly imposed conditions before this, it will after */
  dFunctionReturn(0);
}

dErr dFSGetElements(dFS fs,dInt *n,dInt *restrict*off,s_dRule *restrict*rule,s_dEFS *restrict*efs,dInt *restrict*geomoff,dReal (*restrict*geom)[3])
{

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  if (n)       *n       = fs->nelem;
  if (off)     *off     = fs->off;
  if (rule)    *rule    = fs->rule;
  if (efs)     *efs     = fs->efs;
  if (geomoff) *geomoff = fs->vtxoff;
  if (geom)    *geom    = fs->vtx;
  dFunctionReturn(0);
}

dErr dFSRestoreElements(dFS dUNUSED fs,dInt *n,dInt *restrict*off,s_dRule *restrict*rule,s_dEFS *restrict*efs,dInt *restrict*geomoff,dReal (*restrict*geom)[3])
{

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  if (n)       *n       = 0;
  if (off)     *off     = NULL;
  if (rule)    *rule    = NULL;
  if (efs)     *efs     = NULL;
  if (geomoff) *geomoff = NULL;
  if (geom)    *geom    = NULL;
  dFunctionReturn(0);
}

/** Get arrays sufficiently large to hold the necessary quantities on the highest order element present in the mesh.
*
* @param fs
* @param q Pointer which will hold array of quadrature points, in physical space (not just the local tensor product)
* @param jinv Inverse of element Jacobian evaluated at quadrature points, normally just passed on to dEFSApply
* @param jw Array to store jacobian determinant times quadrature weights at quadrature points
* @param u first array to hold \c bs values per quadrature point
* @param v second array to hold \c bs values per quadrature point
* @param du first array to hold \c 3*bs values per quadrature point
* @param dv second array to hold \c 3*bs values per quadrature point
*/
dErr dFSGetWorkspace(dFS fs,const char name[],dReal (*restrict*q)[3],dReal (*restrict*jinv)[3][3],dReal *restrict*jw,dScalar *restrict*u,dScalar *restrict*v,dScalar *restrict*du,dScalar *restrict*dv)
{
  s_dFSWorkspace *w;
  dInt Q;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  if (!fs->maxQ) {
    Q = 0;
    for (dInt i=0; i<fs->nelem; i++) {
      dInt nnodes;
      err = dRuleGetSize(&fs->rule[i],NULL,&nnodes);dCHK(err);
      if (nnodes > Q) Q = nnodes;
    }
    fs->maxQ = Q;
  }
  Q = fs->maxQ;
  for (dInt i=0; i<dFS_MAX_WORKSPACES; i++) {
    const dInt bs = fs->bs;
    w = &fs->workspace[i];
    switch (w->status) {
      case 0:                   /* Not allocated */
        err = dMallocA7(Q,&w->q,Q,&w->jinv,Q,&w->jw,Q*bs,&w->u,Q*bs,&w->v,Q*bs*3,&w->du,Q*bs*3,&w->dv);dCHK(err);
        w->status = 1;
      case 1:                   /* Available */
        goto found;
      case 2:                   /* Checked out */
        continue;
      default: dERROR(1,"Invalid status %d\n",w->status);
    }
  }
  dERROR(1,"No workspaces available");

  found:
  err = dStrcpyS(w->name,sizeof(w->name),name);dCHK(err);
  w->status = 2;                /* Checked out */
  if (q)    *q    = w->q;
  if (jinv) *jinv = w->jinv;
  if (jw)   *jw   = w->jw;
  if (u)    *u    = w->u;
  if (v)    *v    = w->v;
  if (du)   *du   = w->du;
  if (dv)   *dv   = w->dv;
  dFunctionReturn(0);
}

/* These arrays are persistent for the life of the dFS so we just nullify the pointers */
dErr dFSRestoreWorkspace(dFS fs,const char name[],dReal (*restrict*q)[3],dReal (*restrict*jinv)[3][3],dReal *restrict*jw,dScalar *restrict*u,dScalar *restrict*v,dScalar *restrict*du,dScalar *restrict*dv)
{
  s_dFSWorkspace *w;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
#define dCHECK_NULLIFY(var) do {                                        \
    if (var) {                                                          \
      if (*var == fs->workspace->var) *var = NULL;                      \
      else dERROR(1,"Unable to restore array " #var " because it was not gotten or has been modified"); \
    }                                                                   \
  } while (false)
  for (dInt i=0; i<dFS_MAX_WORKSPACES; i++) {
    w = &fs->workspace[i];
    if (w->status == 2) {
      if (!strcmp(name,w->name)) {
        dCHECK_NULLIFY(q);
        dCHECK_NULLIFY(jinv);
        dCHECK_NULLIFY(jw);
        dCHECK_NULLIFY(u);
        dCHECK_NULLIFY(v);
        dCHECK_NULLIFY(du);
        dCHECK_NULLIFY(dv);
        err = dMemzero(w->name,sizeof(w->name));dCHK(err);
        w->status = 1;
        dFunctionReturn(0);
      }
    }
  }
  dERROR(1,"Did not find a workspace matching `%s' checked out",name);
#undef dCHECK_NULLIFY
  dFunctionReturn(0);
}

dErr dFSGetMatrix(dFS fs,const MatType mtype,Mat *J)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidPointer(J,3);
  err = SlicedGetMatrix(fs->slice,mtype,J);dCHK(err);
  dFunctionReturn(0);dCHK(err);
}

/* We call these directly because otherwise MatGetArray spends huge amounts of time in PetscMallocValidate (unless error
* checking is disabled) */
extern PetscErrorCode MatGetArray_SeqAIJ(Mat A,PetscScalar *array[]);
extern PetscErrorCode MatRestoreArray_SeqAIJ(Mat A,PetscScalar *array[]);

dErr dFSMatSetValuesBlockedExpanded(dFS fs,Mat A,dInt m,const dInt idxm[],dInt n,const dInt idxn[],const dScalar v[],InsertMode imode)
{
  dInt lidxms[128],lidxns[128];
  dScalar lvs[1024],lvts[1024];
  Mat E;
  dInt lm,ln,*lidxm = lidxms,*lidxn = lidxns;
  dInt bs,i,j,li,lj,row,col,cn,*ci,*cj;
  dScalar *lv = lvs,*lvt = lvts,*ca;
  dTruth done;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidHeader(A,MAT_COOKIE,2);
  dValidPointer(idxm,4);
  dValidPointer(idxn,6);
  dValidPointer(v,7);
  bs = fs->bs;
  err = PetscLogEventBegin(dLOG_FSMatSetValuesExpanded,fs,A,0,0);dCHK(err);
  err = MatMAIJGetAIJ(fs->assemblefull?fs->E:fs->Ep,&E);dCHK(err); /* Does not reference so do not destroy or return E */
  err = MatGetRowIJ(E,0,dFALSE,dFALSE,&cn,&ci,&cj,&done);dCHK(err);
  if (!done) dERROR(1,"Could not get indices");
  err = MatGetArray_SeqAIJ(E,&ca);dCHK(err);
  for (i=0,lm=0; i<m; i++) {
    /* Count the number of columns in constraint matrix for each row of input matrix, this will be the total number of
    * rows in result matrix */
    lm += ci[idxm[i]+1] - ci[idxm[i]];
  }
  for (j=0,ln=0; j<n; j++) {
    /* Count the number of columns in constraint matrix for each column of input matrix, this will be the total number of
    * columns in result matrix */
    ln += ci[idxn[j]+1] - ci[idxn[j]];
  }
  if (lm > 128) {err = dMallocA(lm,&lidxm);dCHK(err);}
  if (ln > 128) {err = dMallocA(ln,&lidxn);dCHK(err);}
  if (lm*ln > 1024) {err = dMallocA(lm*ln,&lv);dCHK(err);}
  if (m*ln > 1024) {err = dMallocA(lm*n,&lvt);dCHK(err);}

  /* Expand columns into temporary matrix \a lvt */
  for (j=0,lj=0; j<n; j++) {         /* columns in input matrix */
    col = idxn[j];
    for (dInt k=ci[col]; k<ci[col+1]; k++) { /* become columns in temporary matrix */
      for (i=0; i<m*bs; i++) {       /* every scalar row */
        for (dInt kk=0; kk<bs; kk++) {
          lvt[(i*ln+lj)*bs+kk] = ca[k] * v[(i*n+j)*bs+kk];
        }
      }
      lidxn[lj++] = cj[k];
    }
    err = PetscLogFlops((ci[col+1]-ci[col])*m);dCHK(err);
  }
  /* Expand rows of temporary matrix \a lvt into \a lv */
  for (i=0,li=0; i<m; i++) {         /* rows of temporary matrix */
    row = idxm[i];
    for (dInt k=ci[row]; k<ci[row+1]; k++) { /* become rows of new matrix */
      for (dInt ii=0; ii<bs; ii++) {
        for (j=0; j<ln*bs; j++) { /* each scalar column */
          lv[(li*bs+ii)*ln*bs+j] = ca[k] * lvt[(i*bs+ii)*ln*bs+j];
        }
      }
      lidxm[li++] = cj[k];
    }
    err = PetscLogFlops((ci[row+1]-ci[row])*ln);dCHK(err);
  }

  err = MatRestoreArray_SeqAIJ(E,&ca);dCHK(err);
  err = MatRestoreRowIJ(E,0,dFALSE,dFALSE,&cn,&ci,&cj,&done);dCHK(err);
  if (!done) dERROR(1,"Failed to return indices");
  err = MatSetValuesBlockedLocal(A,lm,lidxm,ln,lidxn,lv,imode);dCHK(err);

  if (lidxm != lidxms) {err = dFree(lidxm);dCHK(err);}
  if (lidxn != lidxns) {err = dFree(lidxn);dCHK(err);}
  if (lv != lvs)       {err = dFree(lv);dCHK(err);}
  if (lvt != lvts)     {err = dFree(lvt);dCHK(err);}
  err = PetscLogEventEnd(dLOG_FSMatSetValuesExpanded,fs,A,0,0);dCHK(err);
  dFunctionReturn(0);
}
