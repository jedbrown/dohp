#include "cont.h"
#include "dohpmesh.h"

static dErr dFSView_Cont(dFS dUNUSED fs,dViewer viewer)
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (ascii) {
    err = PetscViewerASCIIPrintf(viewer,"Continuous Galerkin function space\n");dCHK(err);
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

static dErr dFSContPropogateDegree(dFS fs,const struct dMeshAdjacency *ma)
{
  dInt *deg;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  err = dMallocA(3*ma->nents,&deg);dCHK(err);
  err = dMeshTagGetData(fs->mesh,fs->degreetag,ma->ents,ma->nents,deg,3*ma->nents,dDATA_INT);dCHK(err); /* Get degree everywhere */
  err = dJacobiPropogateDown(fs->jacobi,ma,deg);dCHK(err);
  err = dMeshTagSetData(fs->mesh,fs->degreetag,ma->ents,ma->nents,deg,3*ma->nents,dDATA_INT);dCHK(err);
  err = dFree(deg);dCHK(err);
  dFunctionReturn(0);
}

/**
* Build a scalar continuous function space, perhaps with constraints at non-conforming nodes
*
* @param fs The space to build
*/
static dErr dFSBuildSpace_Cont(dFS fs)
{
  MPI_Comm comm = ((dObject)fs)->comm;
  /* The fact that we aren't using our context here indicates that much/all of the logic here could move up into dFS */
  dUNUSED dFS_Cont *cont = fs->data;
  struct { dInt start,ndofs; } *fieldSpec;
  struct dMeshAdjacency ma;
  dMesh mesh;
  dMeshTag specTag;
  /* MeshListEH v=MLZ,e=MLZ,f=MLZ,r=MLZ; */
  /* MeshListInt in=MLZ,rfo=MLZ,feo=MLZ,rdegree=MLZ; */
  /* dIInt nf; */
  dEntTopology *regTopo;
  dInt i,j,gcnt,ghcnt,*inodes,*xnodes,*deg,*rdeg,*ghostind,nregions;
  dInt *xind,*xstart,*istart,xcnt,*nnz,*pnnz,*regRDeg,*regBDeg;
  dMeshEH *regEnts;
  dEntStatus *status;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,dFS_COOKIE,1);
  mesh = fs->mesh;
  err = dMeshGetAdjacency(mesh,fs->active,&ma);dCHK(err);
  err = dFSContPropogateDegree(fs,&ma);dCHK(err);

  err = dMallocA5(ma.nents*3,&deg,ma.nents*3,&rdeg,ma.nents,&inodes,ma.nents,&xnodes,ma.nents,&status);dCHK(err);
  err = dMeshTagGetData(mesh,fs->degreetag,ma.ents,ma.nents,deg,3*ma.nents,dDATA_INT);dCHK(err);
  err = dMeshTagGetData(mesh,fs->ruletag,ma.ents,ma.nents,rdeg,3*ma.nents,dDATA_INT);dCHK(err);
  /* Fill the arrays \a inodes and \a xnodes with the number of interior and expanded nodes for each
  * (topology,degree) pair */
  err = dJacobiGetNodeCount(fs->jacobi,ma.nents,ma.topo,deg,inodes,xnodes);dCHK(err);
  err = dMeshGetStatus(mesh,ma.nents,ma.ents,status);dCHK(err);

  fs->n = fs->nlocal = fs->nbdofs = 0;
  for (i=0; i<ma.nents; i++) {
    if (!(status[i] & dSTATUS_UNOWNED)) {
      fs->n += inodes[i];      /* It's a global dof on this process */
    }
    fs->nlocal += inodes[i]; /* It's always a local dof */
  }
  err = MPI_Scan(&fs->n,&fs->rstart,1,MPI_INT,MPI_SUM,comm);dCHK(err);
  fs->rstart -= fs->n;
  err = MPI_Allreduce(&fs->n,&fs->N,1,MPI_INT,MPI_SUM,comm);dCHK(err);

  /* Set the global index of the first node associated with every entity (and number of nodes, for a consistency check) */
  err = dMallocA(ma.nents,&fieldSpec);dCHK(err);
  for (i=0,gcnt=fs->rstart; i<ma.nents; i++) {
    if (status[i] & dSTATUS_UNOWNED) {
      fieldSpec[i].start = -1; /* We want to know if these don't get updated */
      fieldSpec[i].ndofs = -1;
    } else {
      fieldSpec[i].start = gcnt;
      fieldSpec[i].ndofs = inodes[i];
      gcnt += inodes[i];
    }
  }

  /* Sync the field spec (make unowned entities have global offsets for their nodes) */
  err = dMeshTagCreateTemp(mesh,"field_spec",2,dDATA_INT,&specTag);dCHK(err);
  err = dMeshTagSetData(mesh,specTag,ma.ents,ma.nents,fieldSpec,ma.nents*(dInt)sizeof(fieldSpec)/(dInt)sizeof(dInt),dDATA_INT);dCHK(err);
  err = dMeshTagBcast(mesh,specTag);dCHK(err);
  err = dMeshTagGetData(mesh,specTag,ma.ents,ma.nents,fieldSpec,ma.nents*(dInt)sizeof(fieldSpec)/(dInt)sizeof(dInt),dDATA_INT);dCHK(err);
  err = dMeshTagDestroy(mesh,specTag);dCHK(err);
  for (i=0; i<ma.nents; i++) { /* Check that all unowned entities were updated */
    if (fieldSpec[i].start < 0 || fieldSpec[i].ndofs < 0) dERROR(1,"Tag exchange did not work");
    if (fieldSpec[i].ndofs != inodes[i]) dERROR(1,"Degree must not have been properly synced");
  }

  /* Set the global indices for all ghost dofs */
  err = dMallocA(fs->nlocal-fs->n,&ghostind);dCHK(err); /* global index of each ghosted dof */
  err = dMallocA(ma.nents,&istart);dCHK(err);          /* First dof of every entity, needed in this form later */
  for (i=0,ghcnt=0; i<ma.nents; i++) {
    if (status[i] & dSTATUS_UNOWNED) {
      for (j=0; j<inodes[i]; j++) ghostind[ghcnt++] = fieldSpec[i].start + j;
    }
    istart[i] = fieldSpec[i].start;
  }
  err = dFree(fieldSpec);dCHK(err);

  err = SlicedCreate(((dObject)fs)->comm,&fs->sliced);dCHK(err);
  if (fs->nlocal-fs->n == 0) {
    err = SlicedSetGhosts(fs->sliced,1,fs->n,0,NULL);dCHK(err);
  } else {
    err = SlicedSetGhosts(fs->sliced,1,fs->n,fs->nlocal-fs->n,ghostind);dCHK(err);
  }
  err = dFree(ghostind);dCHK(err);

  /**
  * At this point the local to global mapping is complete.  Now we need to assemble the constraint matrices which take
  * the local vector to an expanded vector.  If the mesh is conforming and there are no strange boundaries (i.e. slip or
  * normal) the constraint matrix will be boolean (one unit entry per row) in which case an IS would be sufficient.  In
  * the general case, there will be some non-conforming elements and some strange boundaries.  We assemble a full-order
  * constraint matrix and a low-order preconditioning constraint matrix.  The full-order matrix will be used for
  * residual evaluation and matrix-free Jacobian application.  The preconditioning constraint matrix will be used to
  * assemble the low-order preconditioner for the Jacobian (or blocks there of).
  *
  * To generate constraint matrices efficiently, we should preallocate them.  We will make the (possibly poor)
  * assumption that every element with different (must be lower!) order approximation on a downward-adjacent entity will
  * be constrained against all nodes on the adjacent entity.
  *
  * First, we define the expanded space which is just all regions in the domain (we're implementing a scalar space
  * without boundaries).
  */

  nregions = ma.toff[dTYPE_REGION+1] - ma.toff[dTYPE_REGION];
  err = dMallocA6(nregions,&xind,nregions+1,&xstart,nregions,&regTopo,nregions*3,&regRDeg,nregions*3,&regBDeg,nregions,&regEnts);dCHK(err);
  for (i=0,xcnt=0; i<nregions; i++) {
    const dInt ii = ma.toff[dTYPE_REGION] + i;
    xind[i] = ii;               /* Index of entity in full arrays */
    xstart[i] = xcnt;           /* first node on this entity */
    regTopo[i] = ma.topo[ii];
    regEnts[i] = ma.ents[ii];
    for (j=0; j<3; j++) {
      /* Perhaps use a rule which is stronger than required to be Hausdorff (see -dfs_rule_strength) */
      regRDeg[i*3+j] = dMaxInt(rdeg[ii*3+j],deg[ii*3+j] + fs->ruleStrength);
      regBDeg[i*3+j] = deg[ii*3+j];
    }
    xcnt += xnodes[ii];
  }
  fs->m = xstart[nregions] = xcnt;

  err = dMallocA2(xcnt,&nnz,xcnt,&pnnz);dCHK(err);
  err = dJacobiGetConstraintCount(fs->jacobi,nregions,xind,xstart,deg,&ma,nnz,pnnz);dCHK(err);

  /* We don't solve systems with these so it will never make sense for them to use a different format */
  err = MatCreateSeqAIJ(PETSC_COMM_SELF,fs->m,fs->nlocal,1,nnz,&fs->C);dCHK(err);
  err = MatCreateSeqAIJ(PETSC_COMM_SELF,fs->m,fs->nlocal,1,pnnz,&fs->Cp);dCHK(err);
  err = dFree2(nnz,pnnz);dCHK(err);

  err = dJacobiAddConstraints(fs->jacobi,nregions,xind,xstart,istart,deg,&ma,fs->C,fs->Cp);dCHK(err);
  err = dFree(istart);dCHK(err);
  err = dFree5(deg,rdeg,inodes,xnodes,status);dCHK(err);
  err = dMeshRestoreAdjacency(mesh,fs->active,&ma);dCHK(err);

  err = MatAssemblyBegin(fs->C,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyBegin(fs->Cp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(fs->C,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(fs->Cp,MAT_FINAL_ASSEMBLY);dCHK(err);

  /* Get Rule and EFS for domain ents. */
  fs->nelem = nregions;
  err = dMallocA3(nregions,&fs->rule,nregions,&fs->efs,nregions+1,&fs->off);dCHK(err); /* Will be freed by FS */
  err = dMemcpy(fs->off,xstart,(nregions+1)*sizeof(xstart[0]));dCHK(err);
  err = dJacobiGetRule(fs->jacobi,nregions,regTopo,regRDeg,fs->rule);dCHK(err);
  err = dJacobiGetEFS(fs->jacobi,nregions,regTopo,regBDeg,fs->rule,fs->efs);dCHK(err);
  err = dMeshGetVertexCoords(mesh,nregions,regEnts,&fs->vtxoff,&fs->vtx);dCHK(err); /* Should be restored by FS on destroy */
  err = dFree6(xind,xstart,regTopo,regRDeg,regBDeg,regEnts);dCHK(err);

  dFunctionReturn(0);
}

/**
* Create the private structure used by a continuous Galerkin function space.
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
  fs->D = 1;
  fs->data = (void*)fsc;
  fs->ops->view           = dFSView_Cont;
  fs->ops->impldestroy    = dFSDestroy_Cont;
  fs->ops->setfromoptions = dFSSetFromOptions_Cont;
  fs->ops->buildspace     = dFSBuildSpace_Cont;
  dFunctionReturn(0);
}
