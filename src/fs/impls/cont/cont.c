#include "cont.h"
#include <iMesh_extensions.h>
#include <dohpmesh.h>
#include <dohpvec.h>
#include <dohpviewer.h>

static dErr dFSView_Cont(dFS fs,dViewer viewer)
{
  dBool ascii,dhm;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_DHM,&dhm);dCHK(err);
  if (ascii) {
    err = PetscViewerASCIIPrintf(viewer,"Continuous Galerkin function space\n");dCHK(err);
  } else if (dhm) {
    err = dFSView_Cont_DHM(fs,viewer);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr dFSLoadIntoFS_Cont(dViewer viewer,const char fieldname[],dFS fs)
{
  dBool isdhm;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_DHM,&isdhm);dCHK(err);
  if (isdhm) {
    err = dFSLoadIntoFS_Cont_DHM(viewer,fieldname,fs);dCHK(err);
  } else dERROR(PETSC_ERR_SUP,"No support for viewer type \"%s\"",((PetscObject)viewer)->type_name);
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

/**
@note Not collective
*/
static dErr dFSContPropogateDegree(dFS fs,dMeshAdjacency ma)
{
  dPolynomialOrder *deg;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  err = dMallocA(3*ma->nents,&deg);dCHK(err);
  err = dMeshTagGetData(fs->mesh,fs->degreetag,ma->ents,ma->nents,deg,ma->nents,dDATA_INT);dCHK(err); /* Get degree everywhere */
  err = dJacobiPropogateDown(fs->jacobi,ma,deg);dCHK(err);
  err = dMeshTagSetData(fs->mesh,fs->degreetag,ma->ents,ma->nents,deg,ma->nents,dDATA_INT);dCHK(err);
  err = dFree(deg);dCHK(err);
  dFunctionReturn(0);
}

/** Compute a low-bandwidth ordering of explicit entities here (instead of [v,e,f,r])
@pre \a ents must point to an array large enough to hold all explicit ents
     \a ents_a needs the length of this array
The result of this block is that \a ents will have a good ordering
*/
static dErr dMeshPopulateOrderedSet_Private(dMesh mesh,dMeshESH orderedSet,dMeshESH explicitSet,dMeshESH dirichletSet,dMeshESH ghostSet,const MatOrderingType orderingtype)
{
  dScalar         weights[256];
  Mat             madj;
  dTruth          flg;
  dInt           *nnz,*ordering,ordering_a;
  const dInt     *newindices;
  dIInt           adj_a = 0,adj_s,*adjoff=NULL,adjoff_a=0,adjoff_s,ierr,*intdata,ents_s,ents_a,tmp;
  IS              rperm,cperm;
  dMeshEH        *ents = NULL,*adj = NULL;
  dMeshTag        orderTag;
  dErr            err;
  iMesh_Instance  mi;

  dFunctionBegin;
  err = dMeshGetInstance(mesh,&mi);dCHK(err);

  /* Determine number of entities in all sets */
  err = dMeshGetNumEnts(mesh,explicitSet,dTYPE_ALL,dTOPO_ALL,&ents_s);dCHK(err);
  ents_a = ents_s;
  err = dMeshGetNumEnts(mesh,dirichletSet,dTYPE_ALL,dTOPO_ALL,&ents_s);dCHK(err);
  ents_a += ents_s;
  err = dMeshGetNumEnts(mesh,ghostSet,dTYPE_ALL,dTOPO_ALL,&ents_s);dCHK(err);
  ents_a += ents_s;
  err = dMallocA2(ents_a,&ents,ents_a,&intdata);dCHK(err);

  err = dMeshTagCreateTemp(mesh,"ordering",1,dDATA_INT,&orderTag);dCHK(err);

  /* Get all entities packed into the \a ents array */
  err = dMeshGetEnts(mesh,explicitSet,dTYPE_ALL,dTOPO_ALL,ents,ents_a,&ents_s);dCHK(err); /* ents_s = number of explicit ents */
  err = dMeshGetEnts(mesh,dirichletSet,dTYPE_ALL,dTOPO_ALL,ents+ents_s,ents_a-ents_s,&tmp);dCHK(err);
  err = dMeshGetEnts(mesh,ghostSet,dTYPE_ALL,dTOPO_ALL,ents+ents_s+tmp,ents_a-ents_s-tmp,&tmp);dCHK(err);
  if (tmp != 0) dERROR(PETSC_ERR_LIB,"iMesh returned inconsistent count");

  /* Set default data on all ents, entities retaining this value will not be reordered */
  for (dInt i=0; i<ents_a; i++) intdata[i] = -1;
  err = dMeshTagSetData(mesh,orderTag,ents,ents_a,intdata,ents_a,dDATA_INT);dCHK(err);

  /* Get adjacencies of owned explicitly represented entities (not ghosts or dirichlet) */
  iMesh_getEntArrAdj(mi,ents,ents_s,dTYPE_ALL,&adj,&adj_a,&adj_s,&adjoff,&adjoff_a,&adjoff_s,&ierr);dICHK(mi,ierr);
  if (adj_s < ents_s) dERROR(1,"peculiar mesh");
  ordering_a = adj_s;
  err = dMallocA2(ordering_a,&ordering,ents_s,&nnz);dCHK(err); /* enough space for index on all adjacencies */
  for (dInt i=0; i<ents_s; i++) intdata[i] = i; /* Define a reference ordering on primary entities */
  err = dMeshTagSetData(mesh,orderTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);

  err = dMeshTagGetData(mesh,orderTag,adj,adj_s,ordering,ordering_a,dDATA_INT);dCHK(err);
  err = dMeshTagDestroy(mesh,orderTag);dCHK(err);

  for (dInt i=0; i<ents_s; i++) {
    dInt cnt = 0;
    for (dInt j=adjoff[i]; j<adjoff[i+1]; j++) cnt += (ordering[j] >= 0);
    nnz[i] = cnt;
  }
  /* Create matrix for number of owned entities */
  err = MatCreateSeqAIJ(PETSC_COMM_SELF,ents_s,ents_s,0,nnz,&madj);dCHK(err);
  for (dInt i=0; i<256; i++) weights[i] = 1.0; /* Maybe we should use the topology */
  for (dInt i=0; i<ents_s; i++) {
    err = MatSetValues(madj,1,&i,adjoff[i+1]-adjoff[i],ordering+adjoff[i],weights,INSERT_VALUES);dCHK(err);
  }
  err = dFree2(ordering,nnz);dCHK(err);
  free(adjoff); adjoff = NULL; adjoff_a = 0; /* iMesh allocated this memory */
  err = MatAssemblyBegin(madj,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd  (madj,MAT_FINAL_ASSEMBLY);dCHK(err);

  /* Compute reordering */
  err = MatGetOrdering(madj,orderingtype,&rperm,&cperm);dCHK(err);
  err = MatDestroy(madj);dCHK(err);
  err = ISEqual(rperm,cperm,&flg);dCHK(err);
  if (!flg) dERROR(1,"Cannot use ordering");
  err = ISGetIndices(rperm,&newindices);dCHK(err);

  /* Reuse \a adj as a buffer to apply permutation */
  for (dInt i=0; i<ents_s; i++) adj[i] = ents[i];
  for (dInt i=0; i<ents_s; i++) ents[i] = adj[newindices[i]];

  err = ISRestoreIndices(rperm,&newindices);dCHK(err);
  err = ISDestroy(rperm);dCHK(err);
  err = ISDestroy(cperm);dCHK(err);
  free(adj); adj = NULL; adj_a = 0; /* iMesh allocated this memory */

  err = dMeshSetAddEnts(mesh,orderedSet,ents,ents_a);dCHK(err);
  err = dFree2(ents,intdata);dCHK(err);
  dFunctionReturn(0);
}

/* To generate element assembly matrices, we need
* @arg fs[in] Function space, see precondition
* @arg xcnt the size of the expanded space
* @arg idx the MeshAdjacency index of every region (array of length fs->nelem)
* @arg ma MeshAdjacency (array-based connectivity)
* @arg deg integer array of length \c 3*ma.nents which holds the degree of every entity in MeshAdjacency
*
* @note The function space is assumed to be in a state where the following members are valid
*   - mesh
*   - loffsetTag : Local offset
*   - nelem      : Number of elements/regions
*   - off        : Offset of first expanded dof for each elem above (length=nelem+1)
*   - gvec       : Global vector, defines layout
*
* We will create matrices
* @arg[out] E full order element assembly matrix
* @arg[out] Ep preconditioning element assembly matrix (as sparse as possible)
*
* These are preallocated using \a nnz and \a pnnz respectively.
**/
dErr dFSBuildSpace_Cont_CreateElemAssemblyMats(dFS fs,const dInt idx[],const dMeshAdjacency ma,const dInt deg[],Mat *inE,Mat *inEp)
{
  dErr err;
  const dInt *xstart = fs->off;
  dInt nloc,bs,*nnz,*pnnz,*loffset,xcnt;
  Vec   Xclosure,Xloc;
  Mat   E,Ep;

  dFunctionBegin;
  xcnt = xstart[fs->nelem];
  err = dMallocA3(xcnt,&nnz,xcnt,&pnnz,ma->nents,&loffset);dCHK(err);
  err = dMeshTagGetData(fs->mesh,fs->loffsetTag,ma->ents,ma->nents,loffset,ma->nents,dDATA_INT);dCHK(err);

  err = VecDohpGetClosure(fs->gvec,&Xclosure);dCHK(err);
  err = VecGhostGetLocalForm(Xclosure,&Xloc);dCHK(err);
  err = VecGetLocalSize(Xloc,&nloc);dCHK(err);
  err = VecGetBlockSize(Xloc,&bs);dCHK(err);
  nloc /= bs;                   /* nloc now counts nodes */
  err = VecGhostRestoreLocalForm(Xclosure,&Xloc);dCHK(err);
  err = VecDohpRestoreClosure(fs->gvec,&Xclosure);dCHK(err);

  err = dJacobiGetConstraintCount(fs->jacobi,fs->nelem,idx,xstart,loffset,deg,ma,nnz,pnnz);dCHK(err);

  /* We don't solve systems with these so it will never make sense for them to use a different format */
  err = MatCreateSeqAIJ(PETSC_COMM_SELF,xcnt,nloc,1,nnz,&E);dCHK(err);
  err = MatCreateSeqAIJ(PETSC_COMM_SELF,xcnt,nloc,1,pnnz,&Ep);dCHK(err);

  err = dJacobiAddConstraints(fs->jacobi,fs->nelem,idx,xstart,loffset,deg,ma,E,Ep);dCHK(err);

  err = dFree3(nnz,pnnz,loffset);dCHK(err);

  err = MatAssemblyBegin(E,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyBegin(Ep,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(E,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(Ep,MAT_FINAL_ASSEMBLY);dCHK(err);

  err = MatCreateMAIJ(E,bs,inE);dCHK(err);
  err = MatCreateMAIJ(Ep,bs,inEp);dCHK(err);

  err = MatDestroy(E);dCHK(err);
  err = MatDestroy(Ep);dCHK(err);
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
  struct _p_dMeshAdjacency ma;
  dMeshAdjacency         meshAdj;
  dMesh                  mesh;
  iMesh_Instance         mi;
  dEntTopology          *regTopo;
  dInt                  *inodes,*xnodes,*deg,*rdeg,nregions,*bstat,ents_a,ents_s,ghents_s,*intdata,*idx,*ghidx;
  dInt                  *xstart,xcnt,*regBDeg;
  dInt                   bs,n,ngh,ndirichlet,nc,rstart,crstart;
  dIInt                  ierr;
  dMeshEH               *ents,*ghents;
  dEntStatus            *status;
  dErr                   err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  bs   = fs->bs;
  mesh = fs->mesh;
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  err = dMeshGetAdjacency(mesh,fs->activeSet,&meshAdj);dCHK(err);
  err = dMemcpy(&ma,meshAdj,sizeof ma);dCHK(err); /* To have object rather than pointer semantics in this function. */
  err = dFSContPropogateDegree(fs,meshAdj);dCHK(err);

  /* Allocate a workspace that's plenty big, so that we don't have to allocate memory constantly */
  ents_a = ma.nents;
  err = dMallocA4(ents_a,&ents,ents_a,&intdata,ents_a,&idx,ents_a,&ghidx);dCHK(err);

  /* Partition entities in active set into owned explicit, owned Dirichlet, and ghost */
  {
    dInt      nboundaries,ghstart;
    dMeshESH *bdysets;
    iMesh_addEntArrToSet(mi,ma.ents,ma.nents,fs->explicitSet,&ierr);dICHK(mi,ierr);
    /* Move ghost ents from \a explicitSet to \a ghostSet */
    iMesh_getEntitiesRec(mi,fs->explicitSet,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
    err = dMeshPartitionOnOwnership(mesh,ents,ents_s,&ghstart);dCHK(err);
    iMesh_rmvEntArrFromSet(mi,ents+ghstart,ents_s-ghstart,fs->explicitSet,&ierr);dICHK(mi,ierr);
    iMesh_addEntArrToSet(mi,ents+ghstart,ents_s-ghstart,fs->ghostSet,&ierr);dICHK(mi,ierr);
    /* Move owned Dirichlet ents from \a explicitSet to \a dirichletSet */
    err = dMeshGetNumSubsets(mesh,fs->boundaries,1,&nboundaries);dCHK(err);
    if (!nboundaries) goto after_boundaries;
    err = dMallocA2(nboundaries,&bdysets,nboundaries,&bstat);dCHK(err);
    err = dMeshGetSubsets(mesh,fs->boundaries,1,bdysets,nboundaries,NULL);dCHK(err);
    err = dMeshTagSGetData(mesh,fs->bstatusTag,bdysets,nboundaries,bstat,nboundaries,dDATA_INT);dCHK(err);
    for (int i=0; i<nboundaries; i++) {
      if (bstat[i] & dFSBSTATUS_DIRICHLET) {
        iMesh_getEntitiesRec(mi,bdysets[i],dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
        err = dMeshPartitionOnOwnership(mesh,ents,ents_s,&ghstart);dCHK(err);
        iMesh_rmvEntArrFromSet(mi,ents,ghstart,fs->explicitSet,&ierr);dICHK(mi,ierr);
        iMesh_addEntArrToSet(mi,ents,ghstart,fs->dirichletSet,&ierr);dICHK(mi,ierr);
      }
      if (bstat[i] & dFSBSTATUS_WEAK) {
        iMesh_getEntitiesRec(mi,bdysets[i],dTYPE_FACE,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
        iMesh_addEntArrToSet(mi,ents,ents_s,fs->weakFaceSet,&ierr);dICHK(mi,ierr);
      }
    }
    err = dFree2(bdysets,bstat);dCHK(err);
  }
  after_boundaries:

  err = dMeshPopulateOrderedSet_Private(mesh,fs->orderedSet,fs->explicitSet,fs->dirichletSet,fs->ghostSet,fs->orderingtype);dCHK(err);
  {
    char strbuf[dNAME_LEN];
    const char *fsname;
    PetscMPIInt mpirank;
    dInt rank;
    err = MPI_Comm_rank(((dObject)fs)->comm,&mpirank);dCHK(err);
    rank = mpirank;
    err = PetscObjectGetName((PetscObject)fs,&fsname);dCHK(err);
    err = PetscSNPrintf(strbuf,sizeof strbuf,"%s_%s",fsname,dTAG_PARTITION);dCHK(err);
    err = dMeshTagCreate(mesh,strbuf,1,dDATA_INT,&fs->partitionTag);dCHK(err);
    err = dMeshTagSSetData(mesh,fs->partitionTag,&fs->activeSet,1,&rank,1,dDATA_INT);dCHK(err);
    err = PetscSNPrintf(strbuf,sizeof strbuf,"%s_%s",fsname,dTAG_ORDERED_SUBDOMAIN);dCHK(err);
    err = dMeshTagCreate(mesh,strbuf,1,dDATA_INT,&fs->orderedsubTag);dCHK(err);
    err = dMeshTagSSetData(mesh,fs->orderedsubTag,&fs->orderedSet,1,&rank,1,dDATA_INT);dCHK(err);
  }
  err = dMeshGetEnts(fs->mesh,fs->orderedSet,dTYPE_ALL,dTOPO_ALL,ents,ents_a,&ents_s);dCHK(err);
  if (ents_s != ents_a) dERROR(PETSC_ERR_PLIB,"wrong set size");

  /* Get number of nodes for all entities, and parallel status */
  err = dMallocA4(ma.nents*3,&deg,ma.nents*3,&rdeg,ma.nents,&inodes,ma.nents,&status);dCHK(err);
  err = dMeshTagGetData(mesh,fs->degreetag,ma.ents,ma.nents,deg,3*ma.nents,dDATA_INT);dCHK(err);
  err = dMeshTagGetData(mesh,fs->ruletag,ma.ents,ma.nents,rdeg,3*ma.nents,dDATA_INT);dCHK(err);
  /* Fill the arrays \a inodes and \a xnodes with the number of interior and expanded nodes for each
  * (topology,degree) pair */
  err = dJacobiGetNodeCount(fs->jacobi,ma.nents,ma.topo,deg,inodes,NULL);dCHK(err);
  err = dMeshGetStatus(mesh,ma.ents,ma.nents,status);dCHK(err);

  /* Count the number of nodes in each space (explicit, dirichlet, ghost) */
  n = ndirichlet = ngh = 0;
  for (int i=0; i<ma.nents; i++) {
    dInt member;
    dMeshESH exclusive_sets[] = {fs->explicitSet,fs->dirichletSet,fs->ghostSet};
    err = dMeshEntClassifyExclusive(mesh,ma.ents[i],3,exclusive_sets,&member);dCHK(err);
    switch (member) {
      case 0: n          += inodes[i]; break;
      case 1: ndirichlet += inodes[i]; break;
      case 2: ngh        += inodes[i]; break;
      default: dERROR(PETSC_ERR_PLIB,"Should not be possible");
    }
  }
  err = MPI_Scan(&n,&rstart,1,MPIU_INT,MPI_SUM,comm);dCHK(err);
  rstart -= n;
  nc = n + ndirichlet;
  err = MPI_Scan(&nc,&crstart,1,MPIU_INT,MPI_SUM,comm);dCHK(err);
  crstart -= nc;

  fs->n = n;
  fs->nc = nc;
  fs->ngh = ngh;

  {                             /* Set offsets (global, closure, local) of first node associated with every entity */
    dInt i,scan,nentsExplicit,nentsDirichlet;

    /* We assume that orderedSet contains explicitSet+dirichletSet+ghostSet (in that order) */
    err = dMeshGetEnts(mesh,fs->orderedSet,dTYPE_ALL,dTOPO_ALL,ents,ents_a,&ents_s);dCHK(err);
    err = dMeshGetNumEnts(mesh,fs->explicitSet,dTYPE_ALL,dTOPO_ALL,&nentsExplicit);dCHK(err);
    err = dMeshGetNumEnts(mesh,fs->dirichletSet,dTYPE_ALL,dTOPO_ALL,&nentsDirichlet);dCHK(err);
    err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);

    /* global offset */
    for (i=0,scan=rstart; i<nentsExplicit; scan+=inodes[idx[i++]])
      intdata[i] = scan; /* fill \a intdata with the global offset */
    for ( ; i<ents_s; i++) intdata[i] = -1;
    err = dMeshTagSetData(mesh,fs->goffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);

    /* global closure offset */
    for (i=0,scan=crstart; i<nentsExplicit+nentsDirichlet; scan+=inodes[idx[i++]])
      intdata[i] = scan;
    for ( ; i<ents_s; i++) intdata[i] = -1;
    err = dMeshTagSetData(mesh,fs->gcoffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);

    /* local index */
    for (i=0,scan=0; i<ents_s; scan+=inodes[idx[i++]])
      intdata[i] = scan;
    err = dMeshTagSetData(mesh,fs->loffsetTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);

    /* Set a pointer to the ghost portion, we will work with that next */
    ghents   = ents + nentsExplicit + nentsDirichlet;
    ghents_s = ents_s - nentsExplicit - nentsDirichlet;
  }


  /* communicate global and closure offset for ghosts */
  err = dMeshTagBcast(mesh,fs->goffsetTag);dCHK(err);
  err = dMeshTagBcast(mesh,fs->gcoffsetTag);dCHK(err);

  /* Retrieve ghost offsets, to create localupdate. */
  err = dMeshTagGetData(mesh,fs->gcoffsetTag,ghents,ghents_s,intdata,ghents_s,dDATA_INT);dCHK(err);
  for (dInt i=0; i<ghents_s; i++) { /* Paranoia: confirm that all ghost entities were updated. */
    if (intdata[i] < 0) dERROR(1,"Tag exchange did not work");
  }

  /* Set ghost indices of every node using \a ghidx, create global vector. */
  {
    dInt gh=0;
    for (dInt i=0; i<ghents_s; i++) {
      for (dInt j=0; j<inodes[idx[i]]; j++) ghidx[gh++] = intdata[i] + j;
    }
    if (gh != fs->ngh) dERROR(1,"Ghost count inconsistent");
    err = VecCreateDohp(((dObject)fs)->comm,bs,n,nc,ngh,ghidx,&fs->gvec);dCHK(err);
  }

  /* Create fs->bmapping and fs->mapping */
  err = dFSCreateLocalToGlobal_Private(fs,n,nc,ngh,ghidx,rstart);dCHK(err);

  err = VecDohpCreateDirichletCache(fs->gvec,&fs->dcache,&fs->dscat);dCHK(err);

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

  err = dMeshGetEnts(mesh,fs->activeSet,dTYPE_REGION,dTOPO_ALL,ents,ents_a,&ents_s);dCHK(err);
  err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);
  nregions = ents_s;
  err = dMallocA4(nregions+1,&xstart,nregions,&regTopo,nregions*3,&regBDeg,nregions,&xnodes);dCHK(err);
  err = dMeshGetTopo(mesh,ents_s,ents,regTopo);dCHK(err);
  err = dMeshTagGetData(mesh,fs->degreetag,ents,ents_s,regBDeg,nregions*3,dDATA_INT);dCHK(err);
  err = dJacobiGetNodeCount(fs->jacobi,ents_s,regTopo,regBDeg,NULL,xnodes);dCHK(err);

  xcnt = xstart[0] = 0;
  for (dInt i=0; i<nregions; i++) xstart[i+1] = (xcnt += xnodes[i]);

  fs->nelem = nregions;
  err = dMallocA(nregions+1,&fs->off);dCHK(err); /* Will be freed by FS */
  err = dMemcpy(fs->off,xstart,(nregions+1)*sizeof(xstart[0]));dCHK(err);

  err = dMeshGetVertexCoords(mesh,nregions,ents,&fs->vtxoff,&fs->vtx);dCHK(err); /* Should be restored by FS on destroy */
  err = dFree4(xstart,regTopo,regBDeg,xnodes);dCHK(err);

  err = dFSBuildSpace_Cont_CreateElemAssemblyMats(fs,idx,&ma,deg,&fs->E,&fs->Ep);dCHK(err);

  err = dMeshRestoreAdjacency(fs->mesh,fs->activeSet,&meshAdj);dCHK(err); /* Any reason to leave this around for longer? */
  err = dFree4(deg,rdeg,inodes,status);dCHK(err);
  err = dFree4(ents,intdata,idx,ghidx);dCHK(err);
  dFunctionReturn(0);
}

static dErr dFSGetSubElementMeshSize_Cont(dFS fs,dInt *nelem,dInt *nvert,dInt *nconn)
{
  /* dFS_Cont *cont = fs->data; */
  dErr    err;
  dInt    nents,nsub,n;
  dMeshEH *ents;
  dPolynomialOrder *degree;
  Vec     X;

  dFunctionBegin;
  err = dMeshGetNumEnts(fs->mesh,fs->activeSet,dTYPE_REGION,dTOPO_ALL,&nents);dCHK(err);
  err = dMallocA2(nents,&ents,nents,&degree);dCHK(err);
  err = dMeshGetEnts(fs->mesh,fs->activeSet,dTYPE_REGION,dTOPO_ALL,ents,nents,NULL);dCHK(err);
  err = dMeshTagGetData(fs->mesh,fs->degreetag,ents,nents,degree,nents,dDATA_INT);dCHK(err);
  nsub = 0;
  for (dInt e=0; e<nents; e++) {
    /* Since this is for visualization, we represent everything as Q_k Gauss-Lobatto even if it was P_k or otherwise. */
    nsub += dPolynomialOrder1D(degree[e],0) * dPolynomialOrder1D(degree[e],1) * dPolynomialOrder1D(degree[e],2);
  }
  err = VecDohpGetClosure(fs->gvec,&X);dCHK(err);
  err = VecGetLocalSize(X,&n);dCHK(err);
  err = VecDohpRestoreClosure(fs->gvec,&X);dCHK(err);
  *nelem = nsub;
  *nvert = n/fs->bs;
  *nconn = nsub*8;              /* all hex */
  dFunctionReturn(0);
}

static dErr dFSGetSubElementMesh_Cont(dFS fs,dInt nsubelems,dInt nsubconn,dEntTopology subtopo[],dInt suboff[],dInt subind[])
{
  dErr       err;
  //dReal (*nx)[3],(*geom)[3];
  dInt dUNUSED n,sub,subc,*off,*geomoff,nnodes,*ai,*aj;
  dTruth     done;

  dFunctionBegin;
  err = dMemzero(subtopo,sizeof(*subtopo)*nsubelems);dCHK(err);
  err = dMemzero(suboff,sizeof(*suboff)*(nsubelems+1));dCHK(err);
  err = dMemzero(subind,sizeof(*subind)*nsubconn);dCHK(err);

  err = MatGetRowIJ(fs->E,0,PETSC_FALSE,PETSC_FALSE,&nnodes,&ai,&aj,&done);dCHK(err);
  if (!done) dERROR(PETSC_ERR_PLIB,"Element assembly matrix not gotten");
  /*
  In the following, we have to get a bit more than we actually need (topology only).  We'll change the interface if it
  becomes a performance issue.
  */
#if 1
  dERROR(1,"In flux");
#else
  err = dFSGetElements(fs,&n,&off,0,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&nx,0,0,0,0,0,0);dCHK(err); /* space for nodal coordinates (which we won't use) */
  sub = subc = 0;
  for (dInt e=0; e<n; e++) {
    dInt three,P[3];
    err = dEFSGetGlobalCoordinates(&efs[e],(const dReal(*)[3])(geom+geomoff[e]),&three,P,nx);dCHK(err);
    dASSERT(three == 3);
    for (dInt i=0; i<P[0]-1; i++) {
      for (dInt j=0; j<P[1]-1; j++) {
        for (dInt k=0; k<P[2]-1; k++,sub++) {
          dQ1CORNER_CONST_DECLARE(c,rowcol,corners,off[e],nx,P,i,j,k);
          dReal no_warn_unused;           /* \a corners is unused, but I don't want to see a warning about it until I get around */
          no_warn_unused = corners[0][0]; /* to finding a replacement for dQ1CORNER_CONST_DECLARE. */
          subtopo[sub] = dTOPO_HEX;
          suboff[sub] = subc;
          for (dInt l=0; l<8; l++,subc++) {
            dASSERT(0 <= rowcol[l] && rowcol[l] < nnodes);
            subind[subc] = aj[ai[rowcol[l]]];
            if (ai[rowcol[l]+1]-ai[rowcol[l]] != 1) dERROR(PETSC_ERR_SUP,"Element assembly matrix is not boolean");
            if (subc >= nsubconn) dERROR(1,"Insufficient preallocation for connectivity");
          }
        }
      }
    }
  }
  dASSERT(subc == nsubconn);
  suboff[nsubelems] = subc;

  err = dFSRestoreElements(fs,&n,&off,0,0,0,0);dCHK(err);
  err = dFSRestoreWorkspace(fs,__func__,&nx,0,0,0,0,0,0);dCHK(err);
#endif
  err = MatRestoreRowIJ(fs->E,0,PETSC_FALSE,PETSC_FALSE,&n,&ai,&aj,&done);dCHK(err);
  if (!done) dERROR(PETSC_ERR_PLIB,"Element assembly matrix not restored");
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
  fs->data = (void*)fsc;

  fs->ops->view                  = dFSView_Cont;
  fs->ops->impldestroy           = dFSDestroy_Cont;
  fs->ops->setfromoptions        = dFSSetFromOptions_Cont;
  fs->ops->buildspace            = dFSBuildSpace_Cont;
  fs->ops->getsubelementmeshsize = dFSGetSubElementMeshSize_Cont;
  fs->ops->getsubelementmesh     = dFSGetSubElementMesh_Cont;
  fs->ops->loadintofs            = dFSLoadIntoFS_Cont;
  dFunctionReturn(0);
}
