#include "cont.h"
#include "dohpmesh.h"

static dErr dFSView_Cont(dFS fs,dViewer viewer)
{
  dFS_Cont *fsc = fs->data;
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
      err = PetscViewerASCIIPrintf(viewer,"%d/%d element dofs constrained against %d/%d (%d/%d)\n",lm[0],gm[0],lm[1],gm[1],fsc->nowned,fs->N);dCHK(err);
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
  dFS_Cont *fsc = fs->data;
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
  dMesh mesh;
  MeshListEH adj=MLZ,conn=MLZ,aconn=MLZ;
  MeshListInt adjoff=MLZ,connoff=MLZ,aconnoff=MLZ;
  dMeshEH *ents;
  dMeshTag indexTag;
  dIInt ierr;
  dEntTopology topo;
  dInt *eind,*deg,type,cnt,total,nvtx,*aind,nents[4],estart[4];
  dMPIInt rank;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  mesh = fs->mesh;
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  err = MPI_Comm_rank(((dObject)mesh)->comm,&rank);dCHK(err);

  /* Number all entities */
  err = dMeshTagCreateTemp(mesh,"index",1,dDATA_INT,&indexTag);dCHK(err);
  err = dMeshGetNumEnts(mesh,fs->active,dTYPE_ALL,dTOPO_ALL,&total);dCHK(err);
  err = dMeshGetNumEnts(mesh,fs->active,dTYPE_VERTEX,dTOPO_ALL,&nvtx);dCHK(err);
  total -= nvtx;
  err = PetscMalloc2(total,dMeshEH,&ents,total,dInt,&eind);dCHK(err);
  /* We don't propogate to vertices so we'll ignore them, but not preserve normal indexing */
  nents[dTYPE_VERTEX] = 0; estart[dTYPE_VERTEX] = 0;
  for (type=dTYPE_EDGE,cnt=0; type<dTYPE_ALL; type++) { /* Get out all non-vertex entities */
    err = dMeshGetEnts(mesh,fs->active,type,iMesh_ALL_TOPOLOGIES,ents+cnt,total-cnt,&nents[type]);dCHK(err);
    estart[type] = cnt;
    cnt += nents[type];
  }
  if (cnt != total) {
    dMeshEH *allents;
    err = dMalloc(total*sizeof(*allents),&allents);dCHK(err);
    err = dMeshGetEnts(mesh,fs->active,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,allents,total,NULL);dCHK(err);
    for (int i=0; i<total; i++) {
      printf("%d: %p %p\n",i,(void*)ents[i],(void*)allents[i]);
      if (ents[i] != allents[i]) {
        dERROR(1,"mismatch: ents[%d]=%p  allents[%d]=%p\n",i,ents[i],i,allents[i]);
      }
    }
    dERROR(1,"count by type %d does not agree with total count %d",cnt,total);
  }
  for (dInt i=0; i<total; i++) {
    eind[i] = i;
  }
  err = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] nents=(%d %d %d %d)(%d) estart=(%d %d %d %d)\n",rank,nents[0],nents[1],nents[2],nents[3],cnt,estart[0],estart[1],estart[2],estart[3]);
  err = dMeshTagSetData(mesh,indexTag,ents,total,eind,total,dDATA_INT);dCHK(err);

  err = dMalloc(3*total*sizeof(*deg),&deg);dCHK(err);
  err = dMeshTagGetData(mesh,fs->degreetag,ents,total,deg,3*total,dDATA_INT);dCHK(err); /* Get degree everywhere */
  /* connectivity for all but vertices */
  iMesh_getEntArrAdj(mi,ents+nents[dTYPE_VERTEX],total-nents[dTYPE_VERTEX],iBase_VERTEX,MLREF(conn),MLREF(connoff),&ierr);dICHK(mi,ierr);
  for (type=iBase_REGION; type>iBase_EDGE; type--) { /* Propogate degrees downward to facets (not vertices) */
    switch (type) {
      case iBase_REGION: topo = iMesh_HEXAHEDRON; break;
      case iBase_FACE: topo = iMesh_QUADRILATERAL; break;
      case iBase_EDGE: topo = iMesh_LINE_SEGMENT; break;
      default: dERROR(1,"should not happen");
    }
    /* Get adjacent entities */
    iMesh_getEntArrAdj(mi,ents+estart[type],nents[type],type-1,MLREF(adj),MLREF(adjoff),&err);dICHK(mi,err);
    /* And adjacent connectivity */
    iMesh_getEntArrAdj(mi,adj.v,adj.s,iBase_VERTEX,MLREF(aconn),MLREF(aconnoff),&ierr);dICHK(mi,ierr);
    err = dMalloc(adj.s*sizeof(*aind),&aind);dCHK(err);
    err = dMeshTagGetData(mesh,indexTag,adj.v,adj.s,aind,adj.s,dDATA_INT);dCHK(err); /* get indices of adj ents */
    for (dInt i=0; i<nents[type]; i++) {
      err = dJacobiPropogateDown(fs->jacobi,topo,conn.v+connoff.v[estart[type]+i-nents[dTYPE_VERTEX]],deg+(estart[type]+i)*3,
                                 aconn.v+aconnoff.v[adjoff.v[i]],aind+adjoff.v[i],deg);dCHK(err);
    }
    MeshListFree(adj); MeshListFree(adjoff);
    MeshListFree(aconn); MeshListFree(aconnoff);
    err = dFree(aind);dCHK(err);
  }
  err = dMeshTagSetData(mesh,fs->degreetag,ents,total,deg,3*total,dDATA_INT);dCHK(err);
  err = dMeshTagDestroy(mesh,indexTag);dCHK(err);

  MeshListFree(conn); MeshListFree(connoff);
  err = dFree(deg);
  err = PetscFree2(ents,eind);dCHK(err);
  dFunctionReturn(0);
}


static dErr dFSBuildSpace_Cont(dFS fs)
{
  dFS_Cont *cont = fs->data;
  dMesh mesh;
  iMesh_Instance mi;
  dMeshTag indexTag,specTag;
  /* MeshListEH v=MLZ,e=MLZ,f=MLZ,r=MLZ; */
  /* MeshListInt in=MLZ,rfo=MLZ,feo=MLZ,rdegree=MLZ; */
  /* dIInt nf; */
  dMeshEH *ents;
  dInt nents,*ind,*idofs,*xdofs,*deg,cnt,estart[4];
  dEntTopology *topo;
  dEntStatus *status;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  mesh = fs->mesh;
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  err = dFSContPropogateDegree(fs);dCHK(err);

  err = dMeshGetNumEnts(mesh,fs->active,dTYPE_ALL,dTOPO_ALL,&nents);dCHK(err);
  err = dMallocA7(nents,&ents,3*nents,&deg,nents,&topo,nents,&ind,nents,&idofs,nents,&xdofs,nents,&status);dCHK(err);
  cnt = 0;
  for (dEntType type=dTYPE_VERTEX; type<dTYPE_ALL; type++) {
    err = dMeshGetEnts(mesh,fs->active,type,dTOPO_ALL,ents+cnt,nents-cnt,&cont->ne[type]);dCHK(err);
    estart[type] = cnt;
    cnt += cont->ne[type];
  }
  err = dMeshGetStatus(mesh,nents,ents,status);dCHK(err);
  err = dMeshGetTopo(mesh,nents,ents,topo);dCHK(err);
  err = dMeshTagGetData(mesh,fs->degreetag,ents,nents,deg,3*nents,dDATA_INT);
  err = dJacobiGetNodeCount(fs->jacobi,nents,topo,deg,idofs,xdofs);dCHK(err);

  /* Number all active entities with global index, determine how many dofs are local */
  err = dMeshTagCreateTemp(mesh,"local_index",1,dDATA_INT,&indexTag);dCHK(err);

  for (dInt i=0; i<nents; i++) {
    ind[i] = i;
    if (~status[i] & dSTATUS_UNOWNED) {
    }
  }
  err = dMeshTagSetData(mesh,indexTag,ents,nents,ind,nents,dDATA_INT);dCHK(err);

  /* Global index of first dof and number of interior dofs */
  err = dMeshTagCreateTemp(mesh,"field_spec",2,dDATA_INT,&specTag);dCHK(err);
  err = 0;

  // err = dFree6(...);dCHK(err);
  err = dMeshTagDestroy(mesh,indexTag);dCHK(err);
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
  dFS_Cont *fsc;
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
