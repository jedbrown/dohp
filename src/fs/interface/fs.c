#include "private/fsimpl.h"

dErr dFSSetMesh(dFS fs,dMesh mesh,dMeshESH active)
{
  dMesh qmesh;
  dErr  err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(mesh,dMESH_COOKIE,2);
  if (fs->quotient) {
    err = dQuotientGetMesh(fs->quotient,&qmesh);dCHK(err);
    if (mesh != qmesh) fs->quotient = 0; /* The Quotient is stale */
  }
  fs->mesh = mesh;
  fs->active = active;
  dFunctionReturn(0);
}

dErr dFSSetRuleTag(dFS fs,dJacobi jac,dMeshTag rtag)
{

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  fs->ruletag = rtag;
  if (jac && fs->jacobi && fs->jacobi != jac) dERROR(1,"cannot change dJacobi");
  if (jac) fs->jacobi = jac;
  dFunctionReturn(0);
}

dErr dFSSetDegree(dFS fs,dJacobi jac,dMeshTag deg)
{

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(jac,dJACOBI_COOKIE,2);
  fs->degreetag = deg;
  if (jac && fs->jacobi && fs->jacobi != jac) dERROR(1,"cannot change dJacobi");
  if (jac) fs->jacobi = jac;
  dFunctionReturn(0);
}

dErr dFSSetBlockSize(dFS fs,dInt D)
{

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  if (fs->D != 1) dERROR(1,"Block size has already been set");
  fs->newD = D;
  dFunctionReturn(0);
}

/** Register a boundary condition with the function space.  After all boundary conditions are registered and the block
* size set, dFSUpdate can be used.
*
* @param fs function space object
* @param man manifold on which to apply the boundary condition
* @param flip FALSE to use the manifold orientation, TRUE to use the opposite orientation
* @param cfunc constraint function
* @param user context for constraint function
*
* @note The constraint function \b must be a pure function (no side-effects, only writes to it's output matrix) with the
* same definition on every process.  The constraint matrix \b must be invertible and should probably be orthogonal.  The
* number of dofs declared global and local should be the same at every point (this is not actually essential, but it's
* convenient).  The reason it is not declared statically (outside of the function definition) is merely to avoid
* duplicating information that must be kept consistent.
**/
dErr dFSRegisterBoundary(dFS fs,dMeshManifold man,dTruth flip,dFSBoundaryConstraintFunction cfunc,void *user)
{
  dFSBoundary bdy;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  err = dNew(struct _p_dFSBoundary,&bdy);dCHK(err);
  if (!fs->newD) dERROR(1,"Cannot register boundary before setting block size");
  {
    dReal x[3] = {0,0,0},b[3][3] = {{1,0,0},{0,1,0},{0,0,1}},T[dSqr(fs->newD)];
    dInt g;
    /* Run the user function on dummy data to determine the number of global dofs per node (maybe not needed?) */
    err = cfunc(user,x,b,T,&g);dCHK(err);
    bdy->nGlobal = g;
    bdy->nDirichlet = fs->newD - g;
  }
  bdy->manifold = man;
  bdy->flip = flip;
  bdy->cfunc = cfunc;
  bdy->user = user;
  bdy->next = fs->bdylist;      /* Cons with list */
  fs->bdylist = bdy;
  dFunctionReturn(0);
}

/** Commit changes requested by dFSSetBlockSize and dFSRegisterBoundary.
*
* Any Mat/Vec previously obtained from this function space will no longer be compatible with the new dFS.
*
* Internally, new constraint matrices are created that support a vector-valued field (block size) with given boundary
* conditions.  By default, the nonzero pattern of the matrix will be set so that each block is fully coupled using the
* same nonzero pattern as the scalar case.
**/
dErr dFSUpdate(dFS fs)
{
  struct { dFSBoundaryConstraintFunction cfunc; void *user } *bc;
  struct dMeshAdjacency ma;
  dFSBoundary b;
  dMeshTag bTag;
  dReal *tmat;
  dInt i,nb,D,D2;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  nb = 0; b = fs->boundaryList;
  while (b) {
    nb++;
    b = b->next;
  }
  err = dMallocA(nb,&bc);dCHK(err);
  err = dMeshTagCreateTemp(fs->mesh,"boundary_label",1,dDATA_INT,&bTag);dCHK(err);
  for (i=0,b=fs->boundaryList; i<nb; i++,b=b->next) {
    dInt toff[5];
    const dMeshEH *ents;
    const char *orient;
    const dInt *label;
    err = dMeshManifoldGetElements(b->manifold,toff,&ents,&orient);dCHK(err);
    err = dMallocA(toff[dTYPE_ALL],&label);dCHK(err);
    for (dInt j=toff[0]; j<toff[4]; j++) label[j] = i;
    err = dMeshTagSetData(mesh,bTag,ents,toff[4],label,toff[4],dDATA_INT);dCHK(err);
    err = dFree(label);dCHK(err);
    err = dMeshManifoldRestoreElements(b->manifold,toff,&ents,&orient);dCHK(err);
    bc[i].cfunc = b->cfunc;
    bc[i].user = b->user;
  }
  D = fs->newD; D2 = D*D;
  err = dMallocA2(fs->nlocal*D2,&tmat);dCHK(err);
  /* We need the old and new  */
  err = dMeshGetAdjacency(fs->mesh,fs->active,&ma);dCHK(err);
  err = 
  /* \todo Get coordinates and normals at each node */
  for (i=0; i<fs->nlocal; i++) {
    dReal *T = tmat + i*D2;
    
    err = 
  }
  err = dFree(bc);dCHK(err);
  dFunctionReturn(0);
}

static inline dErr dBdyIntersect(dBdyType *a,dBdyType b)
{
  dFunctionBegin;
  /* FIXME: make this handle complicated cases like dBDYTYPE_NORMAL \cap dBDYTYPE_SLIP (= dBDYTYPE_DIRICHLET)*/
  *a = dMax(*a,b);
  dFunctionReturn(0);
}

/** Get an array of boundary types associated with a list of entities
*
* @param fs
* @param[in] nents Number of entities
* @param[in] ents handles
* @param[out] btype array of boundary types, should be at least \a nents in length
*/
dErr dFSGetBoundaryType(dFS fs,dInt nents,const dMeshEH ents[],dBdyType btype[])
{
  dFSBoundary bdy;
  dMeshEH *eh;
  dInt i,n,n_alloc;
  dBdyType *bt;
  dMeshTag tag;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidPointer(ents,3);
  dValidPointer(btype,4);
  bdy = fs->bdylist;
  err = dMeshTagCreateTemp(fs->mesh,"boundary_type",1,dDATA_INT,&tag);dCHK(err);

  n_alloc = nents;
  err = dMallocA2(nents,&eh,nents,&bt);dCHK(err); /* Usually plenty of space */
  for (i=0; i<nents; i++) { bt[i] = dBDYTYPE_NO; }
  err = dMeshTagSetData(fs->mesh,tag,ents,nents,bt,nents,dDATA_INT);dCHK(err); /* Clears all the values */
  while (bdy) {
    err = dMeshGetNumEnts(fs->mesh,bdy->entset,dTYPE_ALL,dTOPO_ALL,&n);dCHK(err); /* Should be superfluous range checking */
    if (n > n_alloc) {
      n_alloc = n;
      err = dFree2(eh,bt);dCHK(err);
      err = dMallocA2(n_alloc,&eh,n_alloc,&bt);dCHK(err);
    }
    /* Overwrite boundary type for all ents in set */
    err = dMeshGetEnts(fs->mesh,bdy->entset,dTYPE_ALL,dTOPO_ALL,eh,n_alloc,&n);dCHK(err);
    err = dMeshTagGetData(fs->mesh,tag,eh,n,bt,n,dDATA_INT);dCHK(err);
    for (i=0; i<n; i++) {
      err = dBdyIntersect(&bt[i],bdy->btype);dCHK(err);
    }
    err = dMeshTagSetData(fs->mesh,tag,eh,n,bt,n,dDATA_INT);dCHK(err);
    bdy = bdy->next;
  }
  err = dFree2(eh,bt);dCHK(err);

  err = dMeshTagGetData(fs->mesh,tag,ents,nents,btype,n,dDATA_INT);dCHK(err);
  err = dMeshTagDestroy(fs->mesh,tag);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSView(dFS fs,dViewer viewer)
{
  dBool iascii;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
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
        err = dMeshGetNumEnts(fs->mesh,fs->active,type,dTOPO_ALL,&nents[type]);dCHK(err);
      }
      err = PetscViewerASCIIPrintf(viewer,"number of vertices=%d edges=%d faces=%d regions=%d\n",nents[0],nents[1],nents[2],nents[3]);dCHK(err);
    }
    {                           /* print aggregate sizes */
      PetscMPIInt gm[2],lm[2];
      lm[0] = fs->m; lm[1] = fs->nlocal; /* set local `element' size and `local' size */
      err = MPI_Reduce(lm,gm,2,MPI_INT,MPI_SUM,0,((dObject)fs)->comm);dCHK(err);
      err = PetscViewerASCIIPrintf(viewer,"on rank 0, %d/%d element dofs constrained against %d/%d local dofs\n",
                                   lm[0],gm[0],lm[1],gm[1]);dCHK(err);
      err = PetscViewerASCIIPrintf(viewer,"on rank 0, %d local (%d owned,%d ghosted) out of %d global",
                                   fs->nlocal,fs->n,fs->nlocal-fs->n,fs->N);dCHK(err);
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
  dValidHeader(fs,dFS_COOKIE,1);
  if (fs->ops->impldestroy) {
    err = (*fs->ops->impldestroy)(fs);dCHK(err);
  }
  if (fs->workspace) {
    s_dFSWorkspace *w = fs->workspace;
    err = dFree7(w->q,w->jinv,w->jw,w->u,w->v,w->du,w->dv);dCHK(err);
    err = dFree(fs->workspace);dCHK(err);
  }
  if (fs->sliced) {
    err = SlicedDestroy(fs->sliced);dCHK(err);
  }
  err = MatDestroy(fs->C);dCHK(err);
  err = MatDestroy(fs->Cp);dCHK(err);
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
  dValidHeader(fs,dFS_COOKIE,1);
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
  dValidHeader(fs,dFS_COOKIE,1);
  dValidPointer(x,2);
  err = MatGetVecs(fs->C,NULL,x);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSCreateGlobalVector(dFS fs,Vec *g)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidPointer(g,2);
  err = SlicedCreateGlobalVector(fs->sliced,g);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSGlobalToExpandedBegin(dFS dUNUSED fs,Vec g,InsertMode imode,Vec dUNUSED x)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(g,VEC_COOKIE,2);
  dValidHeader(x,VEC_COOKIE,4);
  err = VecGhostUpdateBegin(g,imode,SCATTER_FORWARD);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSGlobalToExpandedEnd(dFS fs,Vec g,InsertMode imode,Vec x)
{
  Vec lform;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(g,VEC_COOKIE,2);
  dValidHeader(x,VEC_COOKIE,4);
  err = VecGhostUpdateEnd(g,imode,SCATTER_FORWARD);dCHK(err);
  err = VecGhostGetLocalForm(g,&lform);dCHK(err);
  switch (imode) {
    case INSERT_VALUES:
      err = MatMult(fs->C,lform,x);dCHK(err);
      break;
    case ADD_VALUES:
      err = MatMultAdd(fs->C,lform,x,x);dCHK(err);
      break;
    default:
      dERROR(1,"InsertMode %d not supported",imode);
  }
  err = VecGhostRestoreLocalForm(g,&lform);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSExpandedToGlobal(dFS fs,Vec x,InsertMode imode,Vec g)
{
  Vec lform;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(x,VEC_COOKIE,2);
  dValidHeader(g,VEC_COOKIE,4);
  err = VecGhostGetLocalForm(g,&lform);dCHK(err);
  switch (imode) {
    case INSERT_VALUES:
      err = MatMultTranspose(fs->C,x,lform);dCHK(err);
      break;
    case ADD_VALUES:
      err = MatMultTransposeAdd(fs->C,x,lform,lform);dCHK(err);
      break;
    default:
      dERROR(1,"InsertMode %d not supported",imode);
  }
  err = VecGhostRestoreLocalForm(g,&lform);dCHK(err);
  dFunctionReturn(0);
}

/** Updates the global values, does \b not broadcast the global values back to the ghosts */
dErr dFSExpandedToGlobalBegin(dFS fs,Vec x,InsertMode imode,Vec g)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
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
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(x,VEC_COOKIE,2);
  dValidHeader(g,VEC_COOKIE,4);
  err = VecGhostUpdateEnd(g,ADD_VALUES,SCATTER_REVERSE);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSGetElements(dFS fs,dInt *n,dInt *restrict*off,s_dRule *restrict*rule,s_dEFS *restrict*efs,dInt *restrict*geomoff,dReal (*restrict*geom)[3])
{

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
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
  dValidHeader(fs,dFS_COOKIE,1);
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
* @param u first array to hold D values per quadrature point
* @param v second array to hold D values per quadrature point
* @param du first array to hold 3*D values per quadrature point
* @param dv second array to hold 3*D values per quadrature point
*/
dErr dFSGetWorkspace(dFS fs,dReal (*restrict*q)[3],dReal (*restrict*jinv)[3][3],dReal *restrict*jw,dScalar *restrict*u,dScalar *restrict*v,dScalar *restrict*du,dScalar *restrict*dv)
{
  s_dFSWorkspace *w;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  if (!fs->workspace) {
    dInt Q = 0, D = fs->D;
    err = dNew(s_dFSWorkspace,&w);dCHK(err);
    for (dInt i=0; i<fs->nelem; i++) {
      dInt nnodes;
      err = dRuleGetSize(&fs->rule[i],NULL,&nnodes);dCHK(err);
      if (nnodes > Q) Q = nnodes;
    }
    err = dMallocA7(Q,&w->q,Q,&w->jinv,Q,&w->jw,Q*D,&w->u,Q*D,&w->v,Q*D*3,&w->du,Q*D*3,&w->dv);dCHK(err);
    fs->workspace = w;
  }
  w = fs->workspace;
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
dErr dFSRestoreWorkspace(dFS fs,dReal (*restrict*q)[3],dReal (*restrict*jinv)[3][3],dReal *restrict*jw,dScalar *restrict*u,dScalar *restrict*v,dScalar *restrict*du,dScalar *restrict*dv)
{

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
#define dCHECK_NULLIFY(var) do {                                        \
    if (var) {                                                          \
      if (*var == fs->workspace->var) *var = NULL;                      \
      else dERROR(1,"Unable to restore array " #var " because it was not gotten or has been modified"); \
    }                                                                   \
  } while (false)
  dCHECK_NULLIFY(q);
  dCHECK_NULLIFY(jinv);
  dCHECK_NULLIFY(jw);
  dCHECK_NULLIFY(u);
  dCHECK_NULLIFY(v);
  dCHECK_NULLIFY(du);
  dCHECK_NULLIFY(dv);
#undef dCHECK_NULLIFY
  dFunctionReturn(0);
}

dErr dFSGetMatrix(dFS fs,const MatType mtype,Mat *J)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidPointer(J,3);
  err = SlicedGetMatrix(fs->sliced,mtype,J);dCHK(err);
  dFunctionReturn(0);dCHK(err);
}

/* We call these directly because otherwise MatGetArray spends huge amounts of time in PetscMallocValidate (unless error
* checking is disabled) */
extern PetscErrorCode MatGetArray_SeqAIJ(Mat A,PetscScalar *array[]);
extern PetscErrorCode MatRestoreArray_SeqAIJ(Mat A,PetscScalar *array[]);

dErr dFSMatSetValuesExpanded(dFS fs,Mat A,dInt m,const dInt idxm[],dInt n,const dInt idxn[],const dScalar v[],InsertMode imode)
{
  dInt lidxms[128],lidxns[128];
  dScalar lvs[1024],lvts[1024];
  Mat C;
  dInt lm,ln,*lidxm = lidxms,*lidxn = lidxns;
  dInt i,j,k,li,lj,row,col,cn,*ci,*cj;
  dScalar *lv = lvs,*lvt = lvts,*ca;
  dTruth done;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(A,MAT_COOKIE,2);
  dValidPointer(idxm,4);
  dValidPointer(idxn,6);
  dValidPointer(v,7);
  err = PetscLogEventBegin(dLOG_FSMatSetValuesExpanded,fs,A,0,0);dCHK(err);
  C = (fs->assemblefull) ? fs->C : fs->Cp;
  err = MatGetRowIJ(C,0,dFALSE,dFALSE,&cn,&ci,&cj,&done);dCHK(err);
  if (!done) dERROR(1,"Could not get indices");
  err = MatGetArray_SeqAIJ(C,&ca);dCHK(err);
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
  if (m*ln > 1024) {err = dMallocA(lm*n,&lv);dCHK(err);}

  /* Expand columns into temporary matrix \a lvt */
  for (j=0,lj=0; j<n; j++) {         /* columns in input matrix */
    col = idxn[j];
    for (k=ci[col]; k<ci[col+1]; k++) { /* become columns in temporary matrix */
      for (i=0; i<m; i++) {       /* row */
        lvt[i*ln+lj] = ca[k] * v[i*n+j];
      }
      lidxn[lj++] = cj[k];
    }
    err = PetscLogFlops((ci[col+1]-ci[col])*m);dCHK(err);
  }
  /* Expand rows of temporary matrix \a lvt into \a lv */
  for (i=0,li=0; i<m; i++) {         /* rows of temporary matrix */
    row = idxm[i];
    for (k=ci[row]; k<ci[row+1]; k++) {
      for (j=0; j<ln; j++) {
        lv[li*ln+j] = ca[k] * lvt[i*ln+j];
      }
      lidxm[li++] = cj[k];
    }
    err = PetscLogFlops((ci[row+1]-ci[row])*ln);dCHK(err);
  }

  err = MatRestoreArray_SeqAIJ(C,&ca);dCHK(err);
  err = MatRestoreRowIJ(C,0,dFALSE,dFALSE,&cn,&ci,&cj,&done);dCHK(err);
  if (!done) dERROR(1,"Failed to return indices");
  err = MatSetValuesLocal(A,lm,lidxm,ln,lidxn,lv,imode);dCHK(err);

  if (lidxm != lidxms) {err = dFree(lidxm);dCHK(err);}
  if (lidxn != lidxns) {err = dFree(lidxn);dCHK(err);}
  if (lv != lvs)       {err = dFree(lv);dCHK(err);}
  if (lvt != lvts)     {err = dFree(lvt);dCHK(err);}
  err = PetscLogEventEnd(dLOG_FSMatSetValuesExpanded,fs,A,0,0);dCHK(err);
  dFunctionReturn(0);
}

/**
* Create a function space based on a simpler space
*
* @param fs Function space on which to base the subspace
* @param D Number of dofs per node in the new function space
* @param bdy Boundary conditions to build into the new function space
* @param infs The new function space
*/
dErr dFSCreateSubspace(dFS fs,dInt D,dFSBoundary bdy,dFS *infs)
{
  dFS nfs;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidPointer(infs,3);
  *infs = 0;
  if (fs->bdylist) dERROR(1,"Can only subspace a function space with no boundaries");
  if (D == fs->D && !bdy) dERROR(1,"Attempting to subspace with no modification");
  if (D < fs->D) dERROR(1,"Attempting to subspace by reducing the number of dofs per node");

  dERROR(1,"not implemented");

#if 0
  err = dFSGetBoundaryType(fs,ma.nents,ma.ents,bdytype);dCHK(err);

  err = dMallocA2(ma.nents,&idofs,ma.nents,&bdofs);dCHK(err); /* interior and boundary dofs */
  /**
  * Count the number of dofs that are owned, local, boundary.
  */
  fs->n = fs->nlocal = fs->nbdofs = 0;
  for (type=dTYPE_VERTEX; type<dTYPE_ALL; type++) {
    for (i=ma.toff[type]; i<ma.toff[type+1]; i++) {
      switch (bdytype[i]) {     /* determine if nodes are on a boundary, remove dofs if appropriate */
        case dBDYTYPE_NO:       /* Not on a boundary */
        case dBDYTYPE_WEAK:     /* Neumann boundary, all dofs are in the function space */
          idofs[i] = inodes[i] * cont->D;
          bdofs[i] = 0;
          break;
        case dBDYTYPE_NORMAL:   /* The tangent components are strong (usually 0), normal component is weak */
          idofs[i] = inodes[i];
          bdofs[i] = inodes[i] * (cont->D - 1);
          break;
        case dBDYTYPE_SLIP:     /* The normal component is strong, tangent components are weak */
          idofs[i] = inodes[i] * (cont->D - 1);
          bdofs[i] = inodes[i];
          break;
        case dBDYTYPE_STRONG:   /* Dirichlet conditions */
          idofs[i] = 0;
          bdofs[i] = inodes[i] * cont->D;
          break;
      }
    }
  }
#endif
  *infs = nfs;
  dFunctionReturn(0);
}
