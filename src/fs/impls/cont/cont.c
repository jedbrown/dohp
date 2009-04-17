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

static dErr dFSContPropogateDegree(dFS fs)
{
  dMeshAdjacency ma = fs->meshAdj;
  dInt *deg;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
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
  MPI_Comm               comm  = ((dObject)fs)->comm;
  /* \bug The fact that we aren't using our context here indicates that much/all of the logic here could move up into dFS */
  dUNUSED dFS_Cont      *cont  = fs->data;
  struct dMeshAdjacency  ma;
  dMesh                  mesh;
  iMesh_Instance         mi;
  dEntTopology          *regTopo;
  dInt                  *inodes,*xnodes,*deg,*rdeg,nregions,*bstat,ents_a,ents_s,*intdata,*idx,*ghidx;
  dInt                  *xstart,xcnt,*regRDeg,*regBDeg;
  dInt                   bs,n,ngh,ndirichlet,nghdirichlet,rstart,drstart,crstart,dsplit;
  dIInt                  ierr;
  dMeshEH               *ents;
  dEntStatus            *status;
  dErr                   err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  bs   = fs->bs;
  mesh = fs->mesh;
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  err = dMeshGetAdjacency(mesh,fs->activeSet,&fs->meshAdj);dCHK(err);
  err = dMemcpy(&ma,fs->meshAdj,sizeof ma);dCHK(err); /* To have object rather than pointer semantics in this function. */
  err = dFSContPropogateDegree(fs);dCHK(err);

  /* Allocate a workspace that's plenty big, so that we don't have to allocate memory constantly */
  ents_a = ma.nents;
  err = dMallocA4(ents_a,&ents,ents_a,&intdata,ents_a,&idx,ents_a,&ghidx);dCHK(err);

  /* Partition entities in active set into explicitly represented and Dirichlet */
  {
    dInt      nboundaries;
    dMeshESH *bdysets;
    iMesh_addEntArrToSet(mi,ma.ents,ma.nents,fs->ownedExplicitSet,&ierr);dICHK(mi,ierr); /* Start with all represented explicitly, pretent they are all owned */
    err = dMeshGetNumSubsets(mesh,fs->boundaries,1,&nboundaries);dCHK(err);
    if (!nboundaries) goto after_boundaries;
    err = dMallocA2(nboundaries,&bdysets,nboundaries,&bstat);dCHK(err);
    err = dMeshGetSubsets(mesh,fs->boundaries,1,bdysets,nboundaries,NULL);dCHK(err);
    err = dMeshTagSGetData(mesh,fs->bstatusTag,bdysets,nboundaries,bstat,nboundaries,dDATA_INT);dCHK(err);
    for (int i=0; i<nboundaries; i++) {
      if (bstat[i] & dFSBSTATUS_DIRICHLET) {
        dInt ghstart;
        iMesh_getEntitiesRec(mi,bdysets[i],dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
        iMesh_rmvEntArrFromSet(mi,ents,ents_s,fs->ownedExplicitSet,&ierr);dICHK(mi,ierr);
        err = dMeshPartitionOnOwnership(mesh,ents,ents_s,&ghstart);dCHK(err);
        iMesh_addEntArrToSet(mi,ents,ghstart,fs->ownedDirichletSet,&ierr);dICHK(mi,ierr);
        iMesh_addEntArrToSet(mi,ents+ghstart,ents_s-ghstart,fs->ghostDirichletSet,&ierr);
      }
      if (bstat[i] & dFSBSTATUS_WEAK) {
        iMesh_getEntitiesRec(mi,bdysets[i],dTYPE_FACE,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
        iMesh_addEntArrToSet(mi,ents,ents_s,fs->weakFaceSet,&ierr);dICHK(mi,ierr);
      }
    }
  }
  after_boundaries:

  /* Partition explicit entities into owned and ghost */
  {
    dInt ghstart;
    iMesh_getEntitiesRec(mi,fs->ownedExplicitSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    err = dMeshPartitionOnOwnership(mesh,ents,ents_s,&ghstart);dCHK(err);
    iMesh_rmvEntArrFromSet(mi,ents+ghstart,ents_s-ghstart,fs->ownedExplicitSet,&ierr);dICHK(mi,ierr);
    iMesh_addEntArrToSet(mi,ents+ghstart,ents_s-ghstart,fs->ghostExplicitSet,&ierr);dICHK(mi,ierr);
  }

  /* Get number of nodes for all entities, and parallel status */
  err = dMallocA5(ma.nents*3,&deg,ma.nents*3,&rdeg,ma.nents,&inodes,ma.nents,&xnodes,ma.nents,&status);dCHK(err);
  err = dMeshTagGetData(mesh,fs->degreetag,ma.ents,ma.nents,deg,3*ma.nents,dDATA_INT);dCHK(err);
  err = dMeshTagGetData(mesh,fs->ruletag,ma.ents,ma.nents,rdeg,3*ma.nents,dDATA_INT);dCHK(err);
  /* Fill the arrays \a inodes and \a xnodes with the number of interior and expanded nodes for each
  * (topology,degree) pair */
  err = dJacobiGetNodeCount(fs->jacobi,ma.nents,ma.topo,deg,inodes,xnodes);dCHK(err);
  err = dMeshGetStatus(mesh,ma.ents,ma.nents,status);dCHK(err);

  /* Count the number of nodes in each space (owned, local, owned dirichlet, local dirichlet) */
  n = ngh = ndirichlet = nghdirichlet;
  for (int i=0; i<ma.nents; i++) {
    dIInt isdirichlet,isghostdirichlet,isexplicit,isghostexplicit;
    iMesh_isEntContained(mi,fs->ownedDirichletSet,ma.ents[i],&isdirichlet,&ierr);dICHK(mi,ierr);
    iMesh_isEntContained(mi,fs->ghostDirichletSet,ma.ents[i],&isghostdirichlet,&ierr);dICHK(mi,ierr);
    iMesh_isEntContained(mi,fs->ownedExplicitSet,ma.ents[i],&isexplicit,&ierr);dICHK(mi,ierr);
    iMesh_isEntContained(mi,fs->ghostExplicitSet,ma.ents[i],&isghostexplicit,&ierr);dICHK(mi,ierr);
    if (!!isdirichlet + !!isghostdirichlet + !!isexplicit + !!isghostexplicit != 1) dERROR(1,"should not happen");
    if      (isdirichlet)      ndirichlet   += inodes[i];
    else if (isghostdirichlet) nghdirichlet += inodes[i];
    else if (isexplicit)       n            += inodes[i];
    else if (isghostexplicit)  ngh          += inodes[i];
  }
  err = MPI_Scan(&n,&rstart,1,MPIU_INT,MPI_SUM,comm);dCHK(err);
  rstart -= n;
  err = MPI_Scan(&ndirichlet,&drstart,1,MPIU_INT,MPI_SUM,comm);dCHK(err);
  drstart -= n;
  crstart = rstart + drstart;

  {                             /* Set global index of first node associated with every explicit entity */
    dInt g = rstart,gh;
    iMesh_getEntitiesRec(mi,fs->ownedExplicitSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);
    for (dInt i=0; i<ents_s; i++) {
      intdata[i] = g;
      g += inodes[idx[i]];
    }
    if (g - rstart != n) dERROR(1,"Dohp Error: g does not agree with rstart");
    err = dMeshTagSetData(mesh,fs->goffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);
    iMesh_getEntitiesRec(mi,fs->ghostExplicitSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    for (dInt i=0; i<ents_s; i++) intdata[i] = -1;
    err = dMeshTagSetData(mesh,fs->goffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);
    err = dMeshTagBcast(mesh,fs->goffsetTag);dCHK(err);
    /* Check that all ghost entities were updated */
    err = dMeshTagGetData(mesh,fs->goffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);
    for (dInt i=0; i<ents_s; i++) {
      if (intdata[i] < 0) dERROR(1,"Tag exchange did not work");
    }
    /* Set ghost indices of every node using \a ghidx */
    gh = 0;
    for (dInt i=0; i<ents_s; i++) {
      for (dInt j=0; j<inodes[idx[i]]; j++) ghidx[gh++] = intdata[i] + j;
    }
    err = SlicedCreate(((dObject)fs)->comm,&fs->slice);dCHK(err);
    err = SlicedSetGhosts(fs->slice,bs,n,gh,ghidx);dCHK(err);
  }

  {                             /* Set Dirichlet index of first node associated with every Dirichlet entity, mostly the same as above (horrible non-reuse) */
    dInt gd = drstart,gh;
    iMesh_getEntitiesRec(mi,fs->ownedDirichletSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);
    for (dInt i=0; i<ents_s; i++) {
      intdata[i] = gd;
      gd += inodes[idx[i]];
    }
    if (gd - drstart != ndirichlet) dERROR(1,"Dohp Error: gd does not agree with drstart");
    err = dMeshTagSetData(mesh,fs->gdoffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);
    iMesh_getEntitiesRec(mi,fs->ghostDirichletSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    for (dInt i=0; i<ents_s; i++) intdata[i] = -1;
    err = dMeshTagSetData(mesh,fs->gdoffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);
    err = dMeshTagBcast(mesh,fs->gdoffsetTag);dCHK(err);
    /* Check that all ghost entities were updated */
    err = dMeshTagGetData(mesh,fs->gdoffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);
    for (dInt i=0; i<ents_s; i++) {
      if (intdata[i] < 0) dERROR(1,"Tag exchange did not work");
    }
    /* Set ghost indices of every node using \a ghidx */
    gh = 0;
    for (dInt i=0; i<ents_s; i++) {
      for (dInt j=0; j<inodes[idx[i]]; j++) ghidx[gh++] = intdata[i] + j;
    }
    err = SlicedCreate(((dObject)fs)->comm,&fs->dslice);dCHK(err);
    err = SlicedSetGhosts(fs->dslice,bs,ndirichlet,gh,ghidx);dCHK(err);
    /* Create Dirichlet local and global forms (these are homogeneous for now, the user will have to fill them) */
    err = SlicedCreateGlobalVector(fs->dslice,&fs->d);dCHK(err);
    err = VecGhostGetLocalForm(fs->d,&fs->dl);dCHK(err);
  }

  {                             /* Set global closure index of first node associated with every entity */
    dMeshEH *et;
    dInt et_a,et_s,gc;
    iMesh_getEntitiesRec(mi,fs->ownedExplicitSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    et = ents + ents_s; et_a = ents_a - ents_s;
    iMesh_getEntitiesRec(mi,fs->ownedDirichletSet,dTYPE_ALL,dTOPO_ALL,1,&et,&et_a,&et_s,&ierr);dICHK(mi,ierr);
    ents_s += et_s;
    gc = crstart;
    for (dInt i=0; i<ents_s; i++) {
      intdata[i] = gc;
      gc += inodes[i];
    }
    if (gc - crstart != n + ndirichlet) dERROR(1,"Dohp Error: gc does not agree with crstart");
    err = dMeshTagSetData(mesh,fs->gcoffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);
    /* Create global closure vector */
    err = VecCreateMPI(comm,(n+ndirichlet)*bs,PETSC_DETERMINE,&fs->gc);dCHK(err);
    err = VecSetBlockSize(fs->gc,bs);dCHK(err);
  }

  {                             /* Set local index associated with every entity, global vector first, then Dirichlet */
    /* It is crucial that entities are obtained in the same order as the local and Dirichlet (global plus ghosts) vectors are created. */
    dInt li = 0;
    /* owned explicit */
    iMesh_getEntitiesRec(mi,fs->ownedExplicitSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);
    for (dInt i=0; i<ents_s; i++) {intdata[i] = li; li += inodes[idx[i]];}
    if (li != n) dERROR(1,"Dohp Error: li does not agree with n");
    err = dMeshTagSetData(mesh,fs->loffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);

    /* ghost explicit */
    iMesh_getEntitiesRec(mi,fs->ghostExplicitSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);
    for (dInt i=0; i<ents_s; i++) {intdata[i] = li; li += inodes[idx[i]];}
#if defined(dUSE_DEBUG)
    {
      Vec  gf,lf;
      dInt nl;
      err = SlicedCreateGlobalVector(fs->slice,&gf);dCHK(err);
      err = VecGhostGetLocalForm(gf,&lf);dCHK(err);
      err = VecGetSize(lf,&nl);dCHK(err);
      err = VecGhostRestoreLocalForm(gf,&lf);dCHK(err);
      err = VecDestroy(gf);dCHK(err);
      if (nl != (n+ngh)*bs) dERROR(1,"should not happen");
      if (li*bs != nl) dERROR(1,"Inconsistent sizes, should not happen");
    }
#endif
    dsplit = li;
    err = dMeshTagSetData(mesh,fs->loffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);

    /* owned Dirichlet */
    iMesh_getEntitiesRec(mi,fs->ownedDirichletSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);
    for (dInt i=0; i<ents_s; i++) {intdata[i] = li; li += inodes[idx[i]];}
    if (li-dsplit != ndirichlet) dERROR(1,"Inconsistent sizes, should not happen");
    err = dMeshTagSetData(mesh,fs->loffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);

    /* ghost Dirichlet */
    iMesh_getEntitiesRec(mi,fs->ghostDirichletSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);
    for (dInt i=0; i<ents_s; i++) {intdata[i] = li; li += inodes[idx[i]];}
    {                           /* just debugging */
      dInt nl;
      err = VecGetSize(fs->dl,&nl);dCHK(err);
      if (nl != (ndirichlet+nghdirichlet)*bs) dERROR(1,"should not happen");
      if ((li-dsplit)*bs != nl) dERROR(1,"Inconsistent sizes, should not happen");
    }
    err = dMeshTagSetData(mesh,fs->loffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);
  }

  {
    IS from,to;
    Vec vec;
    /* Create global closure to global scatter (ctog) */
    err = ISCreateStride(comm,n*bs,crstart*bs,1,&from);dCHK(err);
    err = ISCreateStride(PETSC_COMM_SELF,n*bs,0,1,&to);dCHK(err);
    err = SlicedCreateGlobalVector(fs->slice,&vec);dCHK(err);
    err = VecScatterCreate(fs->gc,from,vec,to,&fs->ctog);dCHK(err);
    err = ISDestroy(from);dCHK(err);
    err = ISDestroy(to);dCHK(err);
    err = VecDestroy(vec);dCHK(err);
    /* Create global closure to Dirichlet scatter (ctod) */
    err = ISCreateStride(comm,ndirichlet*bs,(crstart+n)*bs,1,&from);dCHK(err);
    err = ISCreateStride(PETSC_COMM_SELF,ndirichlet*bs,0,1,&to);dCHK(err);
    err = VecScatterCreate(fs->gc,from,fs->d,to,&fs->ctod);dCHK(err);
    err = ISDestroy(from);dCHK(err);
    err = ISDestroy(to);dCHK(err);
  }

  /**
  * At this point the local to global mapping is complete.  Now we need to assemble the constraint matrices which take
  * the local vector to an expanded vector and the local Dirichlet vector to an expanded.  If the mesh is conforming and
  * there are no strange boundaries (i.e. slip or normal) the constraint matrix will be boolean (one unit entry per row)
  * in which case an IS would be sufficient.  In the general case, there will be some non-conforming elements and some
  * strange boundaries.  We assemble a full-order constraint matrix and a low-order preconditioning constraint matrix.
  * The full-order matrix will be used for residual evaluation and matrix-free Jacobian application.  The
  * preconditioning constraint matrix will be used to assemble the low-order preconditioner for the Jacobian (or blocks
  * there of).  Even the full-order matrix is cheap to apply, but it's use in preconditioner assembly significantly
  * impacts sparsity.
  *
  * To generate constraint matrices efficiently, we should preallocate them.  We will make the (possibly poor)
  * assumption that every element with different (must be lower!) order approximation on a downward-adjacent entity will
  * be constrained against all nodes on the adjacent entity.
  */

  iMesh_getEntitiesRec(mi,fs->activeSet,dTYPE_REGION,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
  err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);
  nregions = ents_s;
  err = dMallocA4(nregions+1,&xstart,nregions,&regTopo,nregions*3,&regRDeg,nregions*3,&regBDeg);dCHK(err);
  xcnt = 0;
  for (dInt i=0; i<nregions; i++) {
    const dInt ii = idx[i]; /* Index in MeshAdjacency */
    dInt type;
    xstart[i] = xcnt;              /* first node on this entity */
    regTopo[i] = ma.topo[ii];
    type = iMesh_TypeFromTopology[regTopo[i]];
    for (dInt j=0; j<type && j<3; j++) {
      regRDeg[i*3+j] = dMaxInt(rdeg[ii*3+j],deg[ii*3+j]+fs->ruleStrength);
      regBDeg[i*3+j] = deg[ii*3+j];
    }
    for (dInt j=type; j<3; j++) {
      regRDeg[i*3+2] = 1;
      regBDeg[i*3+2] = 1;
    }
    xcnt += xnodes[ii];
  }
  xstart[nregions] = xcnt;

  {
    dInt *nnz,*pnnz,*dnnz;
    Mat   E,Ep,Ed;
    err = dMallocA3(xcnt,&nnz,xcnt,&pnnz,xcnt,&dnnz);dCHK(err);
    err = dMeshTagGetData(mesh,fs->loffsetTag,ma.ents,ma.nents,intdata,ma.nents,dDATA_INT);dCHK(err);
    /* To generate element assembly matrices, we need
    * \a idx the MeshAdjacency index of every region
    * \a xstart offset in expanded vector of first node associated with this region
    * \a dsplit integer which splits \a istart into local and local Dirichlet pieces
    * \a istart offset in local vectors of first dof associated with each entity (not just regions)
    *    let \c is=istart[e].  If \c is<dsplit then \c is is offset in local vector, if \c is>=dsplit
    *    then \c is-dsplit is the offset in local Dirichlet vector.  The array \a istart is the generic buffer \a intdata
    * \a deg integer array of length \c 3*ma.nents which holds the degree of every entity in MeshAdjacency
    * \a ma MeshAdjacency (array-based connectivity)
    *
    * We will create matrices
    * \a E full order element assembly matrix
    * \a Ep preconditioning element assembly matrix (as sparse as possible)
    * \a Ed Dirichlet element assembly matrix
    *
    * These are preallocated using \a nnz, \a pnnz, \a dnnz respectively.
    **/
    err = dJacobiGetConstraintCount(fs->jacobi,nregions,idx,xstart,dsplit,intdata,deg,&ma,nnz,pnnz,dnnz);dCHK(err);

    /* We don't solve systems with these so it will never make sense for them to use a different format */
    err = MatCreateSeqAIJ(PETSC_COMM_SELF,xcnt,n+ngh,1,nnz,&E);dCHK(err);
    err = MatCreateSeqAIJ(PETSC_COMM_SELF,xcnt,n+ngh,1,pnnz,&Ep);dCHK(err);
    if (ndirichlet+nghdirichlet > 0) {
      err = MatCreateSeqAIJ(PETSC_COMM_SELF,xcnt,ndirichlet+nghdirichlet,1,dnnz,&Ed);dCHK(err);
    } else {
      err = MatCreateSeqAIJ(PETSC_COMM_SELF,xcnt,0,0,NULL,&Ed);dCHK(err);
    }
    err = dFree3(nnz,pnnz,dnnz);dCHK(err);

    err = dJacobiAddConstraints(fs->jacobi,nregions,idx,xstart,dsplit,intdata,deg,&ma,E,Ep,Ed);dCHK(err);
    err = dFree5(deg,rdeg,inodes,xnodes,status);dCHK(err);
    err = dMeshRestoreAdjacency(mesh,fs->activeSet,&fs->meshAdj);dCHK(err); /* Any reason to leave this around for longer? */

    err = MatAssemblyBegin(E,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyBegin(Ep,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyBegin(Ed,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd(E,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd(Ep,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd(Ed,MAT_FINAL_ASSEMBLY);dCHK(err);

    err = MatCreateMAIJ(E,bs,&fs->E);dCHK(err);
    err = MatCreateMAIJ(Ep,bs,&fs->Ep);dCHK(err);
    err = MatCreateMAIJ(Ed,bs,&fs->Ed);dCHK(err);

    err = MatDestroy(E);dCHK(err);
    err = MatDestroy(Ep);dCHK(err);
    err = MatDestroy(Ed);dCHK(err);
  }

  /* Get Rule and EFS for domain ents. */
  fs->nelem = nregions;
  err = dMallocA3(nregions,&fs->rule,nregions,&fs->efs,nregions+1,&fs->off);dCHK(err); /* Will be freed by FS */
  err = dMemcpy(fs->off,xstart,(nregions+1)*sizeof(xstart[0]));dCHK(err);
  err = dJacobiGetRule(fs->jacobi,nregions,regTopo,regRDeg,fs->rule);dCHK(err);
  err = dJacobiGetEFS(fs->jacobi,nregions,regTopo,regBDeg,fs->rule,fs->efs);dCHK(err);
  err = dMeshGetVertexCoords(mesh,nregions,ents,&fs->vtxoff,&fs->vtx);dCHK(err); /* Should be restored by FS on destroy */
  err = dFree4(xstart,regTopo,regRDeg,regBDeg);dCHK(err);
  err = dFree4(ents,intdata,idx,ghidx);dCHK(err);

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
  fs->bs = 1;
  fs->data = (void*)fsc;
  fs->ops->view           = dFSView_Cont;
  fs->ops->impldestroy    = dFSDestroy_Cont;
  fs->ops->setfromoptions = dFSSetFromOptions_Cont;
  fs->ops->buildspace     = dFSBuildSpace_Cont;
  dFunctionReturn(0);
}
