#include "private/fsimpl.h"

dErr dFSSetMesh(dFS fs,dMesh mesh,dMeshESH active,dMeshTag partition)
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
  fs->partition = partition;
  dFunctionReturn(0);
}

dErr dFSSetQuotient(dFS fs,dQuotient quot)
{
  dMesh qmesh;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(quot,dQUOTIENT_COOKIE,2);
  if (fs->mesh) {
    err = dQuotientGetMesh(quot,&qmesh);dCHK(err);
    if (fs->mesh != qmesh) {
      fs->mesh = 0;
      fs->active = 0;
      fs->partition = 0;
    }
  }
  fs->quotient = 0;
  dFunctionReturn(0);
}

dErr dFSSetDegree(dFS fs,dMeshTag deg,dJacobi jac)
{

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  dValidHeader(jac,dJACOBI_COOKIE,3);
  fs->degree = deg;
  fs->jacobi = jac;
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
  bdy->next = fs->bdy_start;    /* prepend to the list */
  fs->bdy_start = bdy;
  dFunctionReturn(0);
}

dErr dFSSetUp(dFS);

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
    if (!fs->setupcalled) {
      err = PetscViewerASCIIPrintf(viewer,"Object has not been set up.\n");dCHK(err);
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
