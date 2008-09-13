#include "cont.h"


static dErr dFSView_Cont(dFS fs,dViewer viewer)
{
  struct dFS_Cont *fsc = fs->data;
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (ascii) {
    err = PetscViewerASCIIPrintf(viewer,"Continuous Galerkin function space\n");dCHK(err);
    {                           /* print aggregate sizes */
      PetscMPIInt gm[2],lm[2];
      lm[0] = fsc->m; lm[1] = fs->n; /* set local `element' size and `local' size */
      err = MPI_Allreduce(lm,gm,2,MPI_INT,MPI_SUM,((dObject)fs)->comm);dCHK(err);
      err = PetscViewerASCIIPrintf(viewer,"%d element dofs constrained against %d local dofs assembled to %d global dofs\n",gm[0],gm[1],fs->N);dCHK(err);
    }
  }
  dFunctionReturn(0);
}

/**
* Calculate the sizes of the global and local vectors, create scatter contexts.  Assemble the constraint matrix for
* element->global maps.
*
* @param fs
*
* @return
*/
static dErr dFSSetFromOptions_Cont(dFS fs)
{
  dMesh mesh = fs->mesh;
  struct dFS_Cont *fsc = fs->data;
  dBool flg;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsHead("Continuous Galerkin options");dCHK(err);
  {
    err = PetscOptionsName("-dfs_cont_constraint_matrix","use explicit SeqAIJ constraint matrix for constraints","None",&flg);dCHK(err);
    if (flg) { fsc->usecmatrix = true; }
  }
  err = PetscOptionsTail();dCHK(err);
  dFunctionReturn(0);
}

static dErr dFSDestroy_Cont(dFS fs)
{
  dErr err;

  dFunctionBegin;
  err = dFree(fs->data);dCHK(err);
  dFunctionReturn(0);
}

static dErr dFSContPropogateDegree(dFS fs)
{
  iMesh_Instance mi;
  MeshListEH rf=MLZ;
  MeshListInt rfo=MLZ;
  dErr err;

  dFunctionBegin;
  err = dMeshGetInstance(fs->mesh,&mi);dCHK(err);
  iMesh_getEntArrAdj(mi,fs->r.v,fs->r.s,iBase_FACE,&rf.v,&rf.a,&rf.s,&rfo.v,&rfo.a,&rfo.s,&err);dICHK(mi,err);
  /* orient faces with respect to regions, propogate anisotropic degree */
#warning much to do here
  MeshListFree(rf); MeshListFree(rfo);
  dFunctionReturn(0);
}

static dErr dFSBuildSpace_Cont(dFS fs)
{
  dErr err;

  dFunctionBegin;
  err = dFSContPropogateDegree(fs);dCHK(err);
  dFunctionReturn(0);
}

/**
* Create the private structure used by a continuous galerkin function space.
*
* This function does not allocate the constraint matrices.
*
* @param fs the function space
*
* @return err
*/
dErr dFSCreate_Cont(dFS fs)
{
  struct dFS_Cont *fsc;
  dErr err;

  dFunctionBegin;
  err = dNewLog(fs,*fsc,&fsc);dCHK(err);
  fs->data = (void*)fsc;
  fs->ops->view           = dFSView_Cont;
  fs->ops->impldestroy    = dFSDestroy_Cont;
  fs->ops->setfromoptions = dFSSetFromOptions_Cont;
  fs->ops->buildspace     = dFSBuildSpace_Cont;
  dFunctionReturn(0);
}
