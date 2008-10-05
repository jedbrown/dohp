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

dErr dFSAddBdy(dFS fs,const char *name,dMeshESH facets,dMeshTag orient,dBool flip,PF constrain)
{
  struct dFSBoundary *bdy;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  err = dNew(struct dFSBoundary,&bdy);dCHK(err);
  err = PetscStrallocpy(name,&bdy->name);dCHK(err);
  bdy->facets = facets;
  bdy->orient = orient;
  bdy->fliporient = flip;
  bdy->constrain = constrain;
  bdy->next = fs->bdylist;      /* Cons with list */
  fs->bdylist = bdy;
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
      err = PetscViewerASCIIPrintf(viewer,"General information about the mesh topology.\n");dCHK(err);
      err = PetscViewerASCIIPrintf(viewer,"number of vertices=%d edges=%d faces=%d regions=%d\n",fs->v.s,fs->e.s,fs->f.s,fs->r.s);dCHK(err);
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
