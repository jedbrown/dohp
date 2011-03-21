#include <dohpfsimpl.h>
#include <dohpvec.h>
#include <dohpstring.h>

extern dErr VecView_Dohp_FSCont(Vec,PetscViewer);

const char *const dFSHomogeneousModes[] = {"HOMOGENOUS","INHOMOGENEOUS","dFSHomogeneousMode","dFS_",0};
const char *const dFSClosureModes[] = {"CLOSURE","INTERIOR","dFSClosureMode","dFS_",0};
const char *const dFSRotateModes[] = {"FORWARD","REVERSE","dFSRotateMode","dFS_ROTATE_",0};

dErr dFSSetMesh(dFS fs,dMesh mesh,dMeshESH active)
{
  dErr  err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(mesh,dMESH_CLASSID,2);
  err = PetscObjectReference((PetscObject)mesh);dCHK(err);
  if (fs->mesh) {err = dMeshDestroy(fs->mesh);dCHK(err);}
  fs->mesh = mesh;
  fs->set.active= active;
  err = dMeshGetTag(mesh,fs->bdyTagName,&fs->tag.boundary);dCHK(err);
  err = dMeshTagCreateTemp(mesh,"boundary_status",1,dDATA_INT,&fs->tag.bstatus);dCHK(err);
  err = dMeshTagCreateTemp(mesh,"boundary_constraint",sizeof(struct dFSConstraintCtx),dDATA_BYTE,&fs->tag.bdyConstraint);dCHK(err);
  err = dMeshTagCreateTemp(mesh,"global_offset",1,dDATA_INT,&fs->tag.goffset);dCHK(err);
  err = dMeshTagCreateTemp(mesh,"local_offset",1,dDATA_INT,&fs->tag.loffset);dCHK(err);
  err = dMeshTagCreate(mesh,"global_closure_offset",1,dDATA_INT,&fs->tag.gcoffset);dCHK(err);
  err = dMeshSetCreate(mesh,dMESHSET_ORDERED,&fs->set.ordered);dCHK(err);
  err = dMeshSetCreate(mesh,dMESHSET_UNORDERED,&fs->set.explicit);dCHK(err);
  err = dMeshSetCreate(mesh,dMESHSET_UNORDERED,&fs->set.dirichlet);dCHK(err);
  err = dMeshSetCreate(mesh,dMESHSET_UNORDERED,&fs->set.ghost);dCHK(err);
  err = dMeshSetCreate(mesh,dMESHSET_UNORDERED,&fs->set.boundaries);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSGetMesh(dFS fs,dMesh *mesh)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(mesh,2);
  if (!fs->mesh) {
    err = dMeshCreate(((dObject)fs)->comm,&fs->mesh);dCHK(err);
    err = PetscObjectIncrementTabLevel((PetscObject)fs->mesh,(PetscObject)fs,1);dCHK(err);
    err = PetscLogObjectParent(fs,fs->mesh);dCHK(err);
  }
  *mesh = fs->mesh;
  dFunctionReturn(0);
}

dErr dFSGetJacobi(dFS fs,dJacobi *jacobi)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(jacobi,2);
  if (!fs->jacobi) {
    err = dJacobiCreate(((dObject)fs)->comm,&fs->jacobi);dCHK(err);
    err = PetscObjectIncrementTabLevel((PetscObject)fs->jacobi,(PetscObject)fs,1);dCHK(err);
    err = PetscLogObjectParent(fs,fs->jacobi);dCHK(err);
  }
  *jacobi = fs->jacobi;
  dFunctionReturn(0);
}

dErr dFSSetRuleTag(dFS fs,dJacobi jac,dMeshTag rtag)
{
  dErr err;
  char *name;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  fs->tag.rule = rtag;
  if (fs->jacobi && fs->jacobi != jac) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"cannot change dJacobi");
  if (!fs->mesh) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"You must call dFSSetMesh() before setting rule tags");
  err = dMeshGetTagName(fs->mesh,rtag,&name);dCHK(err);
  if (!name || (name[0] == '_' && name[1] == '_'))
    dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"The quadrature Rule tag must be persistent, it cannot start with '__'");
  err = dFree(name);dCHK(err);
  if (!fs->jacobi) {
    err = PetscObjectReference((PetscObject)jac);dCHK(err);
    fs->jacobi = jac;
  }
  dFunctionReturn(0);
}

dErr dFSSetDegree(dFS fs,dJacobi jac,dMeshTag deg)
{
  dErr err;
  char *name;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(jac,dJACOBI_CLASSID,2);
  fs->tag.degree = deg;
  if (fs->jacobi && fs->jacobi != jac) dERROR(PETSC_COMM_SELF,1,"cannot change dJacobi");
  if (!fs->mesh) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"You must call dFSSetMesh() before setting rule tags");
  err = dMeshGetTagName(fs->mesh,deg,&name);dCHK(err);
  if (!name || (name[0] == '_' && name[1] == '_'))
    dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"The element Degree tag must be persistent, it cannot start with '__'");
  err = dFree(name);dCHK(err);
  if (!fs->jacobi) {
    err = PetscObjectReference((PetscObject)jac);dCHK(err);
    fs->jacobi = jac;
  }
  dFunctionReturn(0);
}

/** Set the block size for a function space.
*
* @param fs function space
* @param bs block size (number of dofs per node)
**/
dErr dFSSetBlockSize(dFS fs,dInt bs)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  for (dInt i=0; i<fs->bs; i++) {err = dFree(fs->fieldname[i]);dCHK(err);}
  err = dFree(fs->fieldname);dCHK(err);
  err = dCallocA(bs,&fs->fieldname);dCHK(err);
  fs->bs = bs;
  dFunctionReturn(0);
}

dErr dFSGetBlockSize(dFS fs,dInt *bs)
{
  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidIntPointer(bs,2);
  *bs = fs->bs;
  dFunctionReturn(0);
}

/** Set the name of a field managed by the function space.
*
* @param fs function space
* @param fn field number
* @param fname field name
*
* @note You must call dFSsetBlockSize() before this if you have multiple fields.
*/
dErr dFSSetFieldName(dFS fs,dInt fn,const char *fname)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  if (fn < 0 || fs->bs <= fn) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Field number %d out of range",fn);
  if (fs->fieldname[fn]) {err = dFree(fs->fieldname[fn]);dCHK(err);}
  err = PetscStrallocpy(fname,&fs->fieldname[fn]);dCHK(err);
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
  dErr             err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  err = dMeshGetTaggedSet(fs->mesh,fs->tag.boundary,&mid,&bset);dCHK(err);
  err = dFSRegisterBoundarySet(fs,bset,bstat,cfunc,user);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSRegisterBoundarySet(dFS fs,dMeshESH bset,dFSBStatus bstat,dFSConstraintFunction cfunc,void *user)
{
  iMesh_Instance   mi;
  dErr             err;
  dIInt            ierr;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  if (!dFSBStatusValid(bstat)) dERROR(PETSC_COMM_SELF,1,"Boundary status %x invalid",bstat);
  if (dFSBStatusStrongCount(bstat) > fs->bs) dERROR(PETSC_COMM_SELF,1,"Cannot impose strong conditions on more dofs than the block size");
  err = dMeshTagSSetData(fs->mesh,fs->tag.bstatus,&bset,1,&bstat,sizeof(bstat),dDATA_BYTE);dCHK(err);
  if (cfunc) {
    struct dFSConstraintCtx ctx;
    ctx.cfunc = cfunc;
    ctx.user = user;
    err = dMeshTagSSetData(fs->mesh,fs->tag.bdyConstraint,&bset,1,&ctx,sizeof(ctx),dDATA_BYTE);dCHK(err);
  }
  err = dMeshGetInstance(fs->mesh,&mi);dCHK(err);
  iMesh_addEntSet(mi,bset,fs->set.boundaries,&ierr);dICHK(mi,ierr);
  dFunctionReturn(0);
}

dErr dFSView(dFS fs,dViewer viewer)
{
  dBool iascii;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  if (!viewer) {
    err = PetscViewerASCIIGetStdout(((dObject)fs)->comm,&viewer);dCHK(err);
  }
  dValidHeader(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(fs,1,viewer,2);

  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);dCHK(err);
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
        err = dMeshGetNumEnts(fs->mesh,fs->set.active,type,dTOPO_ALL,&nents[type]);dCHK(err);
      }
      err = PetscViewerASCIIPrintf(viewer,"number of vertices=%d edges=%d faces=%d regions=%d\n",nents[0],nents[1],nents[2],nents[3]);dCHK(err);
    }
    {                           /* print aggregate sizes */
      dInt lm[4],gm[4];
      err = MatGetSize(fs->E,&lm[0],&lm[1]);dCHK(err);
      if (lm[0]%fs->bs || lm[1]%fs->bs) dERROR(PETSC_COMM_SELF,1,"Constraint matrix not a multiple of block size, should not happen");
      lm[0] /= fs->bs;
      lm[1] /= fs->bs;
      if (lm[1] != fs->nc) dERROR(PETSC_COMM_SELF,1,"Inconsistent number of closure nodes");
      lm[2] = fs->n;
      lm[3] = fs->ngh;
      err = MPI_Reduce(lm,gm,4,MPIU_INT,MPI_SUM,0,((dObject)fs)->comm);dCHK(err);
      err = PetscViewerASCIIPrintf(viewer,"On rank 0: %d/%d expanded nodes constrained against %d+%d / %d+%d real nodes, %d / %d closure\n",
                                   lm[0],gm[0], lm[2],lm[3], gm[2],gm[3], lm[1],gm[1]);dCHK(err);
      err = PetscViewerASCIIPrintf(viewer,"Block size %d: global dofs %d, ghost dofs %d, closure dofs %d\n",
                                   fs->bs,fs->bs*gm[2],fs->bs*gm[3],fs->bs*gm[1]);dCHK(err);
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

/**
Load the FS associated with a named field at the current time step
**/
dErr dFSLoadIntoFS(PetscViewer viewer,const char fieldname[],dFS fs)
{
  dErr              err;

  dFunctionBegin;
  dValidHeader(viewer,PETSC_VIEWER_CLASSID,1);
  dValidCharPointer(fieldname,2);
  dValidHeader(fs,DM_CLASSID,3);
  if (!fs->ops->loadintofs) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"FS does not support load");
  err = (*fs->ops->loadintofs)(viewer,fieldname,fs);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSDestroy(dFS fs)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  if (--((PetscObject)fs)->refct > 0) dFunctionReturn(0);
  if (fs->ops->impldestroy) {
    err = (*fs->ops->impldestroy)(fs);dCHK(err);
  }
  for (dInt i=0; i<fs->bs; i++) {err = dFree(fs->fieldname[i]);dCHK(err);}
  err = dFree(fs->fieldname);dCHK(err);
  err = VecDestroy(fs->gvec);dCHK(err);
  err = VecDestroy(fs->dcache);dCHK(err);
  err = VecScatterDestroy(fs->dscat);dCHK(err);
  err = MatDestroy(fs->E);dCHK(err);
  err = MatDestroy(fs->Ep);dCHK(err);
  err = ISLocalToGlobalMappingDestroy(fs->bmapping);dCHK(err);
  err = ISLocalToGlobalMappingDestroy(fs->mapping);dCHK(err);
  err = dFree(fs->off);dCHK(err);
  for (struct _dFSIntegrationLink *link=fs->integration,*tmp; link; link=tmp) {
    err = dFree(link->name);dCHK(err);
    err = dFree2(link->rule,link->efs);dCHK(err);
    tmp = link->next;
    err = dFree(link);dCHK(err);
  }

  /* Geometry */
  if (fs->geometry.expanded) {err = VecDestroy(fs->geometry.expanded);dCHK(err);}
  if (fs->geometry.global) {err = VecDestroy(fs->geometry.global);dCHK(err);}
  if (fs->geometry.fs) {err = dFSDestroy(fs->geometry.fs);dCHK(err);}

  /* Nodal Coordinates */
  if (fs->nodalcoord.expanded) {err = VecDestroy(fs->nodalcoord.expanded);dCHK(err);}
  if (fs->nodalcoord.global) {err = VecDestroy(fs->nodalcoord.global);dCHK(err);}
  if (fs->nodalcoord.fs) {err = dFSDestroy(fs->nodalcoord.fs);dCHK(err);}

  err = dMeshDestroy(fs->mesh);dCHK(err);
  err = dJacobiDestroy(fs->jacobi);dCHK(err);
  err = dFree(fs->ops);dCHK(err);
  { /* PetscHeaderDestroy is a macro, need it to destroy the correct ops table, and also need an lvalue */
    DM dm = (DM)fs;
    err = PetscHeaderDestroy(dm);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dFSSetOptionsPrefix(dFS fs,const char prefix[])
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  err = PetscObjectSetOptionsPrefix((PetscObject)fs,prefix);dCHK(err);
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
  dValidHeader(fs,DM_CLASSID,1);
  if (!((dObject)fs)->type_name) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_TYPENOTSET,"Cannot build space");
  if (fs->spacebuilt) dERROR(PETSC_COMM_SELF,1,"The space is already built, rebuilding is not implemented");
  if (fs->ops->buildspace) {
    err = (*fs->ops->buildspace)(fs);dCHK(err);
  }

  if (0) {
    /* Determine the number of elements in which each dof appears */
    err = dFSCreateExpandedVector(fs,&x);dCHK(err);
    err = dFSCreateGlobalVector(fs,&g);dCHK(err);
    err = VecSet(x,1);dCHK(err);
    err = VecZeroEntries(g);dCHK(err);
    err = dFSExpandedToLocal(fs,x,g,ADD_VALUES);dCHK(err);
    err = VecGhostUpdateBegin(g,ADD_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecGhostUpdateEnd(g,ADD_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecDestroy(x);dCHK(err);

    /* \todo Use g to set sparsity pattern */
    err = VecDestroy(g);dCHK(err);
  }

  fs->spacebuilt = dTRUE;
  dFunctionReturn(0);
}

dErr dFSCreateExpandedVector(dFS fs,Vec *x)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(x,2);
  err = MatGetVecs(fs->E,NULL,x);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSCreateGlobalVector(dFS fs,Vec *g)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(g,2);
  /* \todo Could give away gvec if it is only referenced once, but this make handling the composition below very
  * tricky */
  err = VecDuplicate(fs->gvec,g);dCHK(err);
  err = PetscObjectCompose((PetscObject)*g,"dFS",(PetscObject)fs);dCHK(err);
  err = VecSetOperation(*g,VECOP_VIEW,(void(*)(void))VecView_Dohp_FSCont);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSExpandedToLocal(dFS fs,Vec x,Vec l,InsertMode imode)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(x,VEC_CLASSID,2);
  dValidHeader(l,VEC_CLASSID,3);
  switch (imode) {
    case INSERT_VALUES:
      err = MatMultTranspose(fs->E,x,l);dCHK(err);
      break;
    case ADD_VALUES:
      err = MatMultTransposeAdd(fs->E,x,l,l);dCHK(err);
      break;
    default:
      dERROR(PETSC_COMM_SELF,1,"InsertMode %d not supported",imode);
  }
  dFunctionReturn(0);
}

dErr dFSLocalToExpanded(dFS fs,Vec l,Vec x,InsertMode imode)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(l,VEC_CLASSID,2);
  dValidHeader(x,VEC_CLASSID,3);
  switch (imode) {
    case INSERT_VALUES:
      err = MatMult(fs->E,l,x);dCHK(err);
      break;
    case ADD_VALUES:
      err = MatMultAdd(fs->E,l,x,x);dCHK(err);
      break;
    default:
      dERROR(PETSC_COMM_SELF,1,"InsertMode %d not supported",imode);
  }
  dFunctionReturn(0);
}

/** Take the closure vector in natural (unrotated) coordinates and cache the Dirichlet part.
*
* The closure will be returned as is, in unrotated coordinates.  It should be rotated if it's values are to be given to
* a solver component.  This function is used for setting boundary values when they are known analytically.
*
* \see dFSGetClosureCoordinates()
**/
dErr dFSInhomogeneousDirichletCommit(dFS fs,Vec gc)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(gc,VEC_CLASSID,2);
  /* \todo rotate closure vector */
  err = VecScatterBegin(fs->dscat,gc,fs->dcache,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd  (fs->dscat,gc,fs->dcache,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSGlobalToExpanded(dFS fs,Vec g,Vec x,dFSHomogeneousMode hmode,InsertMode imode)
{
  dErr err;
  Vec  gc,lf;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(g,VEC_CLASSID,2);
  dValidHeader(x,VEC_CLASSID,3);
  err = VecDohpGetClosure(g,&gc);dCHK(err);
  switch (hmode) {
    case dFS_HOMOGENEOUS: {     /* project into homogeneous space */
      dInt     n,nc;
      dScalar *a;
      err = VecGetLocalSize(g,&n);dCHK(err);
      err = VecGetLocalSize(gc,&nc);dCHK(err);
      err = VecGetArray(gc,&a);dCHK(err);
      err = dMemzero(a+n,(nc-n)*sizeof(a[0]));dCHK(err);
      err = VecRestoreArray(gc,&a);dCHK(err);
      /* \todo deal with rotations */
    } break;
    case dFS_INHOMOGENEOUS:
      err = VecScatterBegin(fs->dscat,fs->dcache,gc,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
      err = VecScatterEnd  (fs->dscat,fs->dcache,gc,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
      break;
    default: dERROR(PETSC_COMM_SELF,1,"hmode %d unsupported",hmode);
  }
  err = VecGhostUpdateBegin(gc,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecGhostUpdateEnd(gc,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecGhostGetLocalForm(gc,&lf);dCHK(err);
  err = dFSLocalToExpanded(fs,lf,x,imode);dCHK(err);
  err = VecGhostRestoreLocalForm(gc,&lf);dCHK(err);
  err = VecDohpRestoreClosure(g,&gc);dCHK(err);
  dFunctionReturn(0);
}

/** Utility function to move from expanded -> local -> closure -> global
* @param hmode Project resulting vector into this space (only matters for rotated coords because other Dirichlet conditions are in closure)
* @param imode This refers to the expanded->local operation, it does \e not refer to the ghost update which is \e always ADD_VALUES
**/
dErr dFSExpandedToGlobal(dFS fs,Vec x,Vec g,dFSHomogeneousMode hmode,InsertMode imode)
{
  dErr err;
  Vec  gc,lf;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(g,VEC_CLASSID,2);
  dValidHeader(x,VEC_CLASSID,3);
  err = VecDohpGetClosure(g,&gc);dCHK(err);
  err = VecGhostGetLocalForm(gc,&lf);dCHK(err);
  if (imode == ADD_VALUES) {    /* If we want to add, we have to kill off the ghost values otherwise they will be assembled twice */
    dInt     gstart,end;
    dScalar *a;
    err = VecGetLocalSize(gc,&gstart);dCHK(err);
    err = VecGetLocalSize(lf,&end);dCHK(err);
    err = VecGetArray(lf,&a);dCHK(err);
    err = dMemzero(a+gstart,(end-gstart)*sizeof(*a));dCHK(err);
    err = VecRestoreArray(lf,&a);dCHK(err);
  } else if (imode != INSERT_VALUES) dERROR(PETSC_COMM_SELF,1,"unsupported imode");
  err = dFSExpandedToLocal(fs,x,lf,imode);dCHK(err);
  err = VecGhostRestoreLocalForm(gc,&lf);dCHK(err);
  err = VecGhostUpdateBegin(gc,ADD_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecGhostUpdateEnd(gc,ADD_VALUES,SCATTER_REVERSE);dCHK(err);
  if (hmode == dFS_HOMOGENEOUS) { /* \todo project into homogeneous space (for rotated coords) */ }
  err = VecDohpRestoreClosure(g,&gc);dCHK(err);
  dFunctionReturn(0);
}

/** Rotate global vector to/from coordinates where components can be inforced strongly.
*
* \note We currently do not keep track of whether vectors are rotated or not.
*
* dFS_ROTATE_FORWARD: plain cartesian -> global
* dFS_ROTATE_REVERSE: global -> cartesian
*
* dFS_HOMOGENEOUS
*   with FORWARD, means do not recover cached values, enforce homogeneous conditions for these components
*   with REVERSE, means to zero homogeneous part before rotation
* dFS_INHOMOGENEOUS means do nothing special with strongly enforced part of rotated blocks
**/
dErr dFSRotateGlobal(dFS fs,Vec g,dFSRotateMode rmode,dFSHomogeneousMode hmode)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(g,VEC_CLASSID,2);
  err = dFSRotationApply(fs->rot,g,rmode,hmode);dCHK(err);
  dFunctionReturn(0);
}

static dErr dUNUSED dFSIntegrationFindLink(dFS fs,const char *name,struct _dFSIntegrationLink **found)
{
  struct _dFSIntegrationLink *link;

  dFunctionBegin;
  *found = NULL;
  for (link=fs->integration; link; link=link->next) {
    if (!strcmp(name,link->name)) {
      *found = link;
      dFunctionReturn(0);
    }
  }
  dERROR(PETSC_COMM_SELF,1,"Cannot find integration \"%s\"",name);
  dFunctionReturn(0);
}

/** dFSGetDomain - Gets the set containing all entities in the closure of the domain
 *
 */
dErr dFSGetDomain(dFS fs,dMeshESH *domain)
{

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(domain,2);
  *domain = fs->set.active;
  dFunctionReturn(0);
}

dErr dFSGetEFS(dFS fs,dRuleset rset,dInt *n,const dEFS **efs)
{
  dErr             err;
  dInt             ents_a,ents_s;
  dPolynomialOrder *order;
  dEntTopology     *topo;
  dMeshEH          *ents;
  dJacobi          jac;

  dFunctionBegin;
  err = dMeshGetNumEnts(fs->mesh,rset->set,rset->type,rset->topo,&ents_a);dCHK(err);
  err = dMallocA3(ents_a,&ents,ents_a,&topo,ents_a,&order);dCHK(err);
  err = dMeshGetEnts(fs->mesh,rset->set,rset->type,rset->topo,ents,ents_a,&ents_s);dCHK(err);
  err = dMeshGetTopo(fs->mesh,ents_s,ents,topo);dCHK(err);
  err = dMeshTagGetData(fs->mesh,fs->tag.degree,ents,ents_s,order,ents_s,dDATA_INT);dCHK(err);
  /* @bug Only correct for volume integrals. */
  *n = rset->n;
  err = dFSGetJacobi(fs,&jac);dCHK(err);
  err = dJacobiGetEFS(jac,ents_s,topo,order,rset->rules,(dEFS**)efs);dCHK(err);
  err = dFree3(ents,topo,order);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSRestoreEFS(dFS dUNUSED fs,dRuleset dUNUSED rset,dInt *n,const dEFS **efs)
{
  dErr err;

  dFunctionBegin;
  *n = 0;
  err = dFree(*efs);dCHK(err);
  dFunctionReturn(0);
}

static dErr MatGetVecs_DohpFS(Mat A,Vec *x,Vec *y)
{
  dFS fs;
  dErr err;

  dFunctionBegin;
  err = PetscObjectQuery((dObject)A,"DohpFS",(dObject*)&fs);dCHK(err);
  if (!fs) dERROR(PETSC_COMM_SELF,1,"Mat has no composed FS");
  if (x) {err = dFSCreateGlobalVector(fs,x);dCHK(err);}
  if (y) {err = dFSCreateGlobalVector(fs,y);dCHK(err);}
  dFunctionReturn(0);
}

dErr dFSGetMatrix(dFS fs,const MatType mtype,Mat *inJ)
{
  Mat    J;
  dInt   bs,n;
  dErr   err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidCharPointer(mtype,2);
  dValidPointer(inJ,3);
  *inJ = 0;
  bs = fs->bs; n = fs->n;
  err = MatCreate(((dObject)fs)->comm,&J);dCHK(err);
  err = MatSetSizes(J,bs*n,bs*n,PETSC_DETERMINE,PETSC_DETERMINE);dCHK(err);
  err = MatSetType(J,mtype);dCHK(err);
  err = MatSeqBAIJSetPreallocation(J,bs,27,NULL);dCHK(err);         /* \bug incorrect for unstructured meshes */
  err = MatMPIBAIJSetPreallocation(J,bs,27,NULL,25,NULL);dCHK(err); /* \todo this wastes a lot of space in parallel */
  err = MatSeqSBAIJSetPreallocation(J,bs,27,NULL);dCHK(err);
  err = MatMPISBAIJSetPreallocation(J,bs,27,NULL,27,NULL);dCHK(err);
  if (fs->assemblereduced) {
    err = MatSeqAIJSetPreallocation(J,27,NULL);dCHK(err);
    err = MatMPIAIJSetPreallocation(J,27,NULL,25,NULL);dCHK(err);
  } else {
    err = MatSeqAIJSetPreallocation(J,bs*27,NULL);dCHK(err);
    err = MatMPIAIJSetPreallocation(J,bs*27,NULL,bs*25,NULL);dCHK(err);
  }
  err = MatSetBlockSize(J,bs);dCHK(err);
  err = MatSetLocalToGlobalMappingBlock(J,fs->bmapping,fs->bmapping);dCHK(err);
  err = MatSetLocalToGlobalMapping(J,fs->mapping,fs->mapping);dCHK(err);

  /* We want the resulting matrices to be usable with matrix-free operations based on this FS */
  err = PetscObjectCompose((dObject)J,"DohpFS",(dObject)fs);dCHK(err);
  err = MatShellSetOperation(J,MATOP_GET_VECS,(void(*)(void))MatGetVecs_DohpFS);dCHK(err);

  *inJ = J;
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
  dBool  done;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(A,MAT_CLASSID,2);
  dValidPointer(idxm,4);
  dValidPointer(idxn,6);
  dValidPointer(v,7);
  bs = fs->bs;
  err = PetscLogEventBegin(dLOG_FSMatSetValuesExpanded,fs,A,0,0);dCHK(err);
  err = MatMAIJGetAIJ(fs->assemblefull?fs->E:fs->Ep,&E);dCHK(err); /* Does not reference so do not destroy or return E */
  err = MatGetRowIJ(E,0,dFALSE,dFALSE,&cn,&ci,&cj,&done);dCHK(err);
  if (!done) dERROR(PETSC_COMM_SELF,1,"Could not get indices");
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
  if (!done) dERROR(PETSC_COMM_SELF,1,"Failed to return indices");
  if (fs->assemblereduced) {
    dInt brow[lm],bcol[ln];
    dScalar bval[lm*ln];
    for (dInt k=0; k<bs; k++) {
      for (i=0; i<lm; i++) {
        for (j=0; j<ln; j++) {
          bval[i*ln+j] = lv[(i*bs+k)*ln*bs+(j*bs+k)];
        }
      }
      for (i=0; i<lm; i++) brow[i] = lidxm[i]*bs+k;
      for (j=0; j<ln; j++) bcol[j] = lidxn[j]*bs+k;
      err = MatSetValuesLocal(A,lm,brow,ln,bcol,bval,imode);dCHK(err);
    }
  } else {
    err = MatSetValuesBlockedLocal(A,lm,lidxm,ln,lidxn,lv,imode);dCHK(err);
  }

  if (lidxm != lidxms) {err = dFree(lidxm);dCHK(err);}
  if (lidxn != lidxns) {err = dFree(lidxn);dCHK(err);}
  if (lv != lvs)       {err = dFree(lv);dCHK(err);}
  if (lvt != lvts)     {err = dFree(lvt);dCHK(err);}
  err = PetscLogEventEnd(dLOG_FSMatSetValuesExpanded,fs,A,0,0);dCHK(err);
  dFunctionReturn(0);
}
