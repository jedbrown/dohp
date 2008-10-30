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

dErr dFSAddBdy(dFS fs,const char *name,dMeshESH entset,dMeshTag orient,dBool flip,PF constrain)
{
  dFSBoundary bdy;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  err = dNew(struct _p_dFSBoundary,&bdy);dCHK(err);
  err = PetscStrallocpy(name,&bdy->name);dCHK(err);
  bdy->entset = entset;
  bdy->orient = orient;
  bdy->fliporient = flip;
  bdy->constrain = constrain;
  bdy->next = fs->bdylist;      /* Cons with list */
  fs->bdylist = bdy;
  dFunctionReturn(0);
}

static inline dErr dBdyIntersect(dBdyType *a,dBdyType b)
{
  dFunctionBegin;
  /* FIXME: make this handle complicated cases like dBDYTYPE_NORMAL \cap dBDYTYPE_SLIP (= dBDYTYPE_DIRICHLET)*/
  *a = dMax(*a,b);
  dFunctionReturn(0);
}

/** 
* Get an array of boundary types associated with a list of entities
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
  dMesh mesh = fs->mesh;
  iMesh_Instance mi;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  /**
  * FIXME:
  *
  */
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  if (fs->ops->buildspace) {
    err = (*fs->ops->buildspace)(fs);dCHK(err);
  }
  fs->spacebuilt = true;
  dFunctionReturn(0);
}

dErr dFSCreateLocalVector(dFS fs,Vec *U)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidPointer(U,2);
  *U = 0;
  if (fs->sliced) {
    err = SlicedCreateLocalVector(fs->sliced,U);dCHK(err);
  } else {
    dERROR(1,"no sliced");
  }
  dFunctionReturn(0);
}

dErr dFSCreateGlobalVector(dFS fs,Vec *U)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(U,VEC_COOKIE,2);
  if (fs->sliced) {
    err = SlicedCreateGlobalVector(fs->sliced,U);dCHK(err);
  } else {
    dERROR(1,"no sliced");
  }
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
