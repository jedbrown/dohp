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
  err = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&ascii);dCHK(err);
  err = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDHM,&dhm);dCHK(err);
  if (ascii) {
    err = PetscViewerASCIIPrintf(viewer,"Continuous Galerkin function space\n");dCHK(err);
  } else if (dhm) {
    err = dFSView_Cont_DHM(fs,viewer);dCHK(err);
  } else dERROR(((PetscObject)viewer)->comm,PETSC_ERR_SUP,"viewer type \"%s\"",((PetscObject)viewer)->type_name);
  dFunctionReturn(0);
}

static dErr dFSLoadIntoFS_Cont(dViewer viewer,const char fieldname[],dFS fs)
{
  dBool isdhm;
  dErr err;

  dFunctionBegin;
  err = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERDHM,&isdhm);dCHK(err);
  if (isdhm) {
    err = dFSLoadIntoFS_Cont_DHM(viewer,fieldname,fs);dCHK(err);
  } else dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"No support for viewer type \"%s\"",((PetscObject)viewer)->type_name);
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
    err = PetscOptionsBool("-dfs_cont_constraint_matrix","use explicit SeqAIJ constraint matrix for constraints","None",flg=dFALSE,&flg,NULL);dCHK(err);
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
static dErr dFSContPropagateDegree(dFS fs,dMeshAdjacency ma)
{
  dPolynomialOrder *deg;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  err = dMallocA(ma->nents,&deg);dCHK(err);
  err = dMeshTagGetData(fs->mesh,fs->tag.degree,ma->ents,ma->nents,deg,ma->nents,dDATA_INT);dCHK(err); /* Get degree everywhere */
  err = dJacobiPropagateDown(fs->jacobi,ma,deg);dCHK(err);
  err = dMeshTagSetData(fs->mesh,fs->tag.degree,ma->ents,ma->nents,deg,ma->nents,dDATA_INT);dCHK(err);
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
  dBool           flg;
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
  if (tmp != 0) dERROR(PETSC_COMM_SELF,PETSC_ERR_LIB,"iMesh returned inconsistent count");

  /* Set default data on all ents, entities retaining this value will not be reordered */
  for (dInt i=0; i<ents_a; i++) intdata[i] = -1;
  err = dMeshTagSetData(mesh,orderTag,ents,ents_a,intdata,ents_a,dDATA_INT);dCHK(err);

  /* Get adjacencies of owned explicitly represented entities (not ghosts or dirichlet) */
  iMesh_getEntArrAdj(mi,ents,ents_s,dTYPE_ALL,&adj,&adj_a,&adj_s,&adjoff,&adjoff_a,&adjoff_s,&ierr);dICHK(mi,ierr);
  if (adj_s < ents_s) dERROR(PETSC_COMM_SELF,1,"peculiar mesh");
  ordering_a = adj_s;
  err = dMallocA2(ordering_a,&ordering,ents_s,&nnz);dCHK(err); /* enough space for index on all adjacencies */
  for (dInt i=0; i<ents_s; i++) intdata[i] = i; /* Define a reference ordering on primary entities */
  err = dMeshTagSetData(mesh,orderTag,ents,ents_s,intdata,ents_s,dDATA_INT);dCHK(err);

  err = dMeshTagGetData(mesh,orderTag,adj,adj_s,ordering,ordering_a,dDATA_INT);dCHK(err);
  err = dMeshTagDestroy(mesh,orderTag);dCHK(err);

  for (dInt i=0; i<ents_s; i++) {
    dInt cnt = 0;
    for (dInt j=adjoff[i]; j<adjoff[i+1]; j++) cnt += (ordering[j] >= 0);
    nnz[i] = cnt + 1;           // The 1 is for the diagonal
  }
  /* Create matrix for number of owned entities */
  err = MatCreateSeqAIJ(PETSC_COMM_SELF,ents_s,ents_s,0,nnz,&madj);dCHK(err);
  for (dInt i=0; i<256; i++) weights[i] = 1.0; /* Maybe we should use the topology */
  for (dInt i=0; i<ents_s; i++) {
    err = MatSetValue(madj,i,i,1.0,INSERT_VALUES);dCHK(err); // @todo Fix MatGetRowIJ_SeqAIJ_Inode_Symmetric so the diagonal is not required
    err = MatSetValues(madj,1,&i,adjoff[i+1]-adjoff[i],ordering+adjoff[i],weights,INSERT_VALUES);dCHK(err);
  }
  err = dFree2(ordering,nnz);dCHK(err);
  free(adjoff); adjoff = NULL; adjoff_a = 0; /* iMesh allocated this memory */
  err = MatAssemblyBegin(madj,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd  (madj,MAT_FINAL_ASSEMBLY);dCHK(err);

  /* Compute reordering */
  err = MatGetOrdering(madj,orderingtype,&rperm,&cperm);dCHK(err);
  err = MatDestroy(&madj);dCHK(err);
  err = ISEqual(rperm,cperm,&flg);dCHK(err);
  if (!flg) dERROR(PETSC_COMM_SELF,1,"Cannot use ordering");
  err = ISGetIndices(rperm,&newindices);dCHK(err);

  /* Reuse \a adj as a buffer to apply permutation */
  for (dInt i=0; i<ents_s; i++) adj[i] = ents[i];
  for (dInt i=0; i<ents_s; i++) ents[i] = adj[newindices[i]];

  err = ISRestoreIndices(rperm,&newindices);dCHK(err);
  err = ISDestroy(&rperm);dCHK(err);
  err = ISDestroy(&cperm);dCHK(err);
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
dErr dFSBuildSpace_Cont_CreateElemAssemblyMats(dFS fs,const dInt idx[],const dMeshAdjacency ma,const dPolynomialOrder deg[],Mat *inE,Mat *inEp)
{
  dErr err;
  const dInt *xstart = fs->off;
  dInt nloc,bs,*nnz,*pnnz,*loffset,xcnt;
  Vec   Xclosure,Xloc;
  Mat   E,Ep;

  dFunctionBegin;
  if (!fs->off) dERROR(((dObject)fs)->comm,PETSC_ERR_ARG_WRONGSTATE,"Must set fs->off before calling this function");
  xcnt = xstart[fs->nelem];
  err = dMallocA3(xcnt,&nnz,xcnt,&pnnz,ma->nents,&loffset);dCHK(err);
  err = dMeshTagGetData(fs->mesh,fs->tag.loffset,ma->ents,ma->nents,loffset,ma->nents,dDATA_INT);dCHK(err);

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
  if (1) {
    Ep = E;
    err = PetscObjectReference((PetscObject)Ep);dCHK(err);
  }

  err = dJacobiAddConstraints(fs->jacobi,fs->nelem,idx,xstart,loffset,deg,ma,E,Ep);dCHK(err);

  err = dFree3(nnz,pnnz,loffset);dCHK(err);

  err = MatAssemblyBegin(E,MAT_FINAL_ASSEMBLY);dCHK(err);
  if (E != Ep) {err = MatAssemblyBegin(Ep,MAT_FINAL_ASSEMBLY);dCHK(err);}
  err = MatAssemblyEnd(E,MAT_FINAL_ASSEMBLY);dCHK(err);
  if (E != Ep) {err = MatAssemblyEnd(Ep,MAT_FINAL_ASSEMBLY);dCHK(err);}

  err = MatCreateMAIJ(E,bs,inE);dCHK(err);
  if (E == Ep) {
    *inEp = *inE;
    err = PetscObjectReference((PetscObject)*inEp);dCHK(err);
  } else {err = MatCreateMAIJ(Ep,bs,inEp);dCHK(err);}

  err = MatDestroy(&E);dCHK(err);
  err = MatDestroy(&Ep);dCHK(err);
  dFunctionReturn(0);
}



/**
* Build a scalar continuous function space, perhaps with constraints at non-conforming nodes
*
* @param fs The space to build
*/
dErr dFSBuildSpace_Cont(dFS fs)
{
  dErr err;
  dMeshAdjacency meshAdj;
  dMesh mesh;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  err = dFSGetMesh(fs,&mesh);dCHK(err);
  err = dMeshGetAdjacency(mesh,fs->set.active,&meshAdj);dCHK(err);
  err = dFSContPropagateDegree(fs,meshAdj);dCHK(err);

  err = dFSPopulatePartitionedSets_Private(fs,meshAdj);dCHK(err);

  err = dMeshPopulateOrderedSet_Private(mesh,fs->set.ordered,fs->set.explicit,fs->set.dirichlet,fs->set.ghost,fs->orderingtype);dCHK(err);
  {
    char strbuf[dNAME_LEN];
    const char *fsname;
    PetscMPIInt mpirank;
    dInt rank;
    err = MPI_Comm_rank(((dObject)fs)->comm,&mpirank);dCHK(err);
    rank = mpirank;
    err = PetscObjectGetName((PetscObject)fs,&fsname);dCHK(err);
    err = PetscSNPrintf(strbuf,sizeof strbuf,"%s_%s",fsname,dTAG_PARTITION);dCHK(err);
    err = dMeshTagCreate(mesh,strbuf,1,dDATA_INT,&fs->tag.partition);dCHK(err);
    err = dMeshTagSSetData(mesh,fs->tag.partition,&fs->set.active,1,&rank,1,dDATA_INT);dCHK(err);
    err = PetscSNPrintf(strbuf,sizeof strbuf,"%s_%s",fsname,dTAG_ORDERED_SUBDOMAIN);dCHK(err);
    err = dMeshTagCreate(mesh,strbuf,1,dDATA_INT,&fs->tag.orderedsub);dCHK(err);
    err = dMeshTagSSetData(mesh,fs->tag.orderedsub,&fs->set.ordered,1,&rank,1,dDATA_INT);dCHK(err);
  }
  err = dFSBuildSpaceWithOrderedSet_Private(fs,meshAdj);dCHK(err);
  err = dMeshRestoreAdjacency(fs->mesh,fs->set.active,&meshAdj);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSBuildSpaceWithOrderedSet_Private(dFS fs,dMeshAdjacency meshAdj)
{
  MPI_Comm               comm  = ((dObject)fs)->comm;
  /* \bug The fact that we aren't using our context here indicates that much/all of the logic here could move up into dFS */
  dUNUSED dFS_Cont      *cont  = fs->data;
  struct _p_dMeshAdjacency ma;
  dMesh                  mesh;
  iMesh_Instance         mi;
  dEntTopology          *regTopo;
  dPolynomialOrder      *deg,*regBDeg;
  dInt                  *inodes,*xnodes,nregions,ents_a,ents_s,ghents_s,*idx;
  dInt                  *xstart,xcnt;
  dInt                   rstart,crstart;
  dMeshEH               *ents,*ghents;
  dErr                  err;

  dFunctionBegin;
  mesh = fs->mesh;
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  err = dMemcpy(&ma,meshAdj,sizeof ma);dCHK(err); /* To have object rather than pointer semantics in this function. */

  /* Allocate a workspace that's plenty big, so that we don't have to allocate memory constantly */
  ents_a = ma.nents;
  err = dMallocA4(ents_a,&ents,ents_a,&idx,ents_a,&deg,ents_a,&inodes);dCHK(err);

  err = dMeshGetEnts(fs->mesh,fs->set.ordered,dTYPE_ALL,dTOPO_ALL,ents,ents_a,&ents_s);dCHK(err);
  if (ents_s != ents_a) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"wrong set size");

  /* Get number of nodes for all entities */
  err = dMeshTagGetData(mesh,fs->tag.degree,ma.ents,ma.nents,deg,ma.nents,dDATA_INT);dCHK(err);
  /* Fill the \a inodes array with the number of interior nodes for each (topology,degree) pair */
  err = dJacobiGetNodeCount(fs->jacobi,ma.nents,ma.topo,deg,inodes,NULL);dCHK(err);

  {
    dInt counts[3],rstarts[3];
    err = dMeshClassifyCountInt(mesh,ma.nents,ma.ents,inodes,3,(const dMeshESH[]){fs->set.explicit,fs->set.dirichlet,fs->set.ghost},counts);dCHK(err);
    err = MPI_Scan(counts,rstarts,3,MPIU_INT,MPI_SUM,comm);dCHK(err);
    for (dInt i=0; i<3; i++) rstarts[i] -= counts[i];
    fs->n   = counts[0];
    fs->nc  = counts[0] + counts[1];
    fs->ngh = counts[2];
    rstart  = rstarts[0];
    crstart = rstarts[0] + rstarts[1];
  }

  {
    dInt ghstart;
    err = dFSBuildSpaceOffsets_Private(fs,ma.indexTag,inodes,rstart,crstart,ents_s,ents,&ghstart);dCHK(err);
    ghents = ents + ghstart;
    ghents_s = ents_s - ghstart;
  }

  err = dFSBuildSpaceVectors_Private(fs,ma.indexTag,inodes,rstart,ghents_s,ghents);dCHK(err);

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

  err = dMeshGetEnts(mesh,fs->set.active,dTYPE_REGION,dTOPO_ALL,ents,ents_a,&ents_s);dCHK(err);
  err = dMeshTagGetData(mesh,ma.indexTag,ents,ents_s,idx,ents_s,dDATA_INT);dCHK(err);
  nregions = ents_s;
  err = dMallocA4(nregions+1,&xstart,nregions,&regTopo,nregions,&regBDeg,nregions,&xnodes);dCHK(err);
  err = dMeshGetTopo(mesh,ents_s,ents,regTopo);dCHK(err);
  err = dMeshTagGetData(mesh,fs->tag.degree,ents,ents_s,regBDeg,nregions,dDATA_INT);dCHK(err);
  err = dJacobiGetNodeCount(fs->jacobi,ents_s,regTopo,regBDeg,NULL,xnodes);dCHK(err);

  xcnt = xstart[0] = 0;
  for (dInt i=0; i<nregions; i++) xstart[i+1] = (xcnt += xnodes[i]);

  fs->nelem = nregions;
  err = dMallocA(nregions+1,&fs->off);dCHK(err); /* Will be freed by FS */
  err = dMemcpy(fs->off,xstart,(nregions+1)*sizeof(xstart[0]));dCHK(err);

  err = dFree4(xstart,regTopo,regBDeg,xnodes);dCHK(err);

  err = dFSBuildSpace_Cont_CreateElemAssemblyMats(fs,idx,&ma,deg,&fs->E,&fs->Ep);dCHK(err);

  err = dFree4(ents,idx,deg,inodes);dCHK(err);
  fs->spacebuilt = dTRUE;
  dFunctionReturn(0);
}

static dErr dFSGetSubElementMeshSize_Cont(dFS fs,dInt *nelem,dInt *nvert,dInt *nconn)
{
  /* dFS_Cont *cont = fs->data; */
  dErr    err;
  dInt    nents,nsub,n,bs;
  dMeshEH *ents;
  dPolynomialOrder *degree;
  Vec     X;

  dFunctionBegin;
  err = dMeshGetNumEnts(fs->mesh,fs->set.active,dTYPE_REGION,dTOPO_ALL,&nents);dCHK(err);
  err = dMallocA2(nents,&ents,nents,&degree);dCHK(err);
  err = dMeshGetEnts(fs->mesh,fs->set.active,dTYPE_REGION,dTOPO_ALL,ents,nents,NULL);dCHK(err);
  err = dMeshTagGetData(fs->mesh,fs->tag.degree,ents,nents,degree,nents,dDATA_INT);dCHK(err);
  nsub = 0;
  for (dInt e=0; e<nents; e++) {
    /* Since this is for visualization, we represent everything as Q_k Gauss-Lobatto even if it was P_k or otherwise. */
    nsub += dPolynomialOrder1D(degree[e],0) * dPolynomialOrder1D(degree[e],1) * dPolynomialOrder1D(degree[e],2);
  }
  err = VecDohpGetClosure(fs->gvec,&X);dCHK(err);
  err = VecGetLocalSize(X,&n);dCHK(err);
  err = VecDohpRestoreClosure(fs->gvec,&X);dCHK(err);
  err = dFree2(ents,degree);dCHK(err);
  err = dFSGetBlockSize(fs,&bs);dCHK(err);
  *nelem = nsub;
  *nvert = n/bs;
  *nconn = nsub*8;              /* all hex */
  dFunctionReturn(0);
}

static dErr dFSGetSubElementMesh_Cont(dFS fs,dInt nsubelems,dInt nsubconn,dEntTopology subtopo[],dInt suboff[],dInt subind[])
{
  dErr     err;
  dMeshESH domain;
  dRuleset ruleset;
  dRulesetIterator iter;
  const dInt *ai,*aj;
  dInt     sub,subc,nnz;
  dBool    done;
  Mat      E1;

  dFunctionBegin;
  err = dMemzero(subtopo,sizeof(*subtopo)*nsubelems);dCHK(err);
  err = dMemzero(suboff,sizeof(*suboff)*(nsubelems+1));dCHK(err);
  err = dMemzero(subind,sizeof(*subind)*nsubconn);dCHK(err);
  err = MatMAIJGetAIJ(fs->E,&E1);dCHK(err);
  err = MatGetRowIJ(E1,0,PETSC_FALSE,PETSC_FALSE,&nnz,&ai,&aj,&done);dCHK(err);
  if (!done) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Element assembly matrix not gotten");

  /* Traverse the domain extracting the connectivity of each sub-element.  It is okay to use the scalar FS in place of
   * the coordinate FS because we only extract topology.
   */
  err = dFSGetDomain(fs,&domain);dCHK(err);
  err = dFSGetPreferredQuadratureRuleSet(fs,domain,dTYPE_REGION,dTOPO_ALL,dQUADRATURE_METHOD_SPARSE,&ruleset);dCHK(err);
  err = dRulesetCreateIterator(ruleset,fs,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,NULL,NULL);dCHK(err);

  sub = subc = 0;
  while (dRulesetIteratorHasPatch(iter)) {
    dInt P;
    const dInt *rowcol;
    err = dRulesetIteratorSetupElement(iter);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter,&P,&rowcol,NULL,NULL);dCHK(err);
    dASSERT(P == 8);            /* Only implemented for hex elements */
    subtopo[sub] = dTOPO_HEX;
    suboff[sub] = subc;
    for (dInt i=0; i<P; i++,subc++) {
      static const dInt itaps_to_tensor[8] = {0,4,6,2,1,5,7,3}; /* Ugly to put this here */
      const dInt ti = itaps_to_tensor[i];
      subind[subc] = aj[ai[rowcol[ti]]];
      if (ai[rowcol[ti]+1]-ai[rowcol[ti]] != 1) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"Element assembly matrix is not boolean");
      if (subc >= nsubconn) dERROR(PETSC_COMM_SELF,1,"Insufficient preallocation for connectivity");
    }
    err = dRulesetIteratorRestorePatchAssembly(iter,&P,&rowcol,NULL,NULL);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
    sub++;
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = dRulesetIteratorDestroy(&iter);dCHK(err);
  err = dRulesetDestroy(&ruleset);dCHK(err);

  dASSERT(subc == nsubconn);
  suboff[nsubelems] = subc;
  err = MatRestoreRowIJ(E1,0,PETSC_FALSE,PETSC_FALSE,&nnz,&ai,&aj,&done);dCHK(err);
  if (!done) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Element assembly matrix not restored");
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
