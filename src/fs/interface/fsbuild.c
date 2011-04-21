#include <dohpfsimpl.h>
#include <dohpvec.h>
#include <iMesh_extensions.h>

/* Create fs->bmapping and fs->mapping from given layout.
*/
dErr dFSCreateLocalToGlobal_Private(dFS fs,dInt n,dInt nc,dInt ngh,dInt *ghidx,dInt rstart)
{
  dErr    err;
  Vec     g,gc,lf;
  dInt    i,bs,*globals;
  dScalar *a;

  dFunctionBegin;
  /* Create block local to global mapping */
  err = dFSCreateGlobalVector(fs,&g);dCHK(err);
  err = VecGetBlockSize(g,&bs);dCHK(err);
  err = VecDohpGetClosure(g,&gc);dCHK(err);
  err = VecSet(gc,-1);dCHK(err);
  err = VecSet(g,1);dCHK(err);
  err = VecGhostUpdateBegin(gc,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecGhostUpdateEnd(gc,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecGhostGetLocalForm(gc,&lf);dCHK(err);
  err = dMallocA(nc+ngh,&globals);dCHK(err);
  err = VecGetArray(lf,&a);dCHK(err);
  /* \a a is a mask determining whether a value is represented in the global system (1) or not (-1) */
  for (i=0; i<n; i++) {
    if (a[i*bs] != 1) dERROR(PETSC_COMM_SELF,1,"should not happen");
    globals[i] = rstart+i;
  }
  for ( ; i<nc; i++) {
    if (a[i*bs] != -1) dERROR(PETSC_COMM_SELF,1,"should not happen");
    globals[i] = -(rstart + i);
  }
  for ( ; i<nc+ngh; i++) {
    globals[i] = signbit(a[i*bs]) * ghidx[i-nc];
  }
  err = VecRestoreArray(lf,&a);dCHK(err);
  err = VecGhostRestoreLocalForm(gc,&lf);dCHK(err);
  err = VecDohpRestoreClosure(g,&gc);dCHK(err);
  err = VecDestroy(&g);dCHK(err);
  err = ISLocalToGlobalMappingCreate(((dObject)fs)->comm,nc+ngh,globals,PETSC_OWN_POINTER,&fs->bmapping);dCHK(err);
  /* Don't free \a globals because we used the no-copy variant, so the IS takes ownership. */
  {
    dInt *sglobals;        /* Scalar globals */
    err = dMallocA((nc+ngh)*bs,&sglobals);dCHK(err);
    for (i=0; i<(nc+ngh)*bs; i++) sglobals[i] = globals[i/bs]*bs + i%bs;
    err = ISLocalToGlobalMappingCreate(((dObject)fs)->comm,(nc+ngh)*bs,sglobals,PETSC_OWN_POINTER,&fs->mapping);dCHK(err);
    /* Don't free \a sglobals because we used the no-copy variant, so the IS takes ownership. */
  }
  dFunctionReturn(0);
}

/* Partition entities in active set into owned explicit, owned Dirichlet, and ghost */
dErr dFSPopulatePartitionedSets_Private(dFS fs,dMeshAdjacency meshadj)
{
  dInt      nboundaries,ghstart,ents_a=0,ents_s;
  dFSBStatus *bstat;
  dMeshESH *bdysets;
  dMesh mesh;
  iMesh_Instance mi;
  dErr err;
  dIInt ierr;
  dMeshEH *ents;

  dFunctionBegin;
  err = dFSGetMesh(fs,&mesh);dCHK(err);
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  iMesh_addEntArrToSet(mi,meshadj->ents,meshadj->nents,fs->set.explicit,&ierr);dICHK(mi,ierr);
  /* Move ghost ents from \a explicitSet to \a ghostSet */
  iMesh_getEntitiesRec(mi,fs->set.explicit,dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
  err = dMeshPartitionOnOwnership(mesh,ents,ents_s,&ghstart);dCHK(err);
  iMesh_rmvEntArrFromSet(mi,ents+ghstart,ents_s-ghstart,fs->set.explicit,&ierr);dICHK(mi,ierr);
  iMesh_addEntArrToSet(mi,ents+ghstart,ents_s-ghstart,fs->set.ghost,&ierr);dICHK(mi,ierr);
  /* Move owned Dirichlet ents from \a explicitSet to \a dirichletSet */
  err = dMeshGetNumSubsets(mesh,fs->set.boundaries,0,&nboundaries);dCHK(err);
  err = dMallocA2(nboundaries,&bdysets,nboundaries,&bstat);dCHK(err);
  err = dMeshGetSubsets(mesh,fs->set.boundaries,0,bdysets,nboundaries,NULL);dCHK(err);
  err = dMeshTagSGetData(mesh,fs->tag.bstatus,bdysets,nboundaries,bstat,nboundaries,dDATA_INT);dCHK(err);
  for (int i=0; i<nboundaries; i++) {
    if (bstat[i] & dFSBSTATUS_DIRICHLET) {
      iMesh_getEntitiesRec(mi,bdysets[i],dTYPE_ALL,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
      err = dMeshPartitionOnOwnership(mesh,ents,ents_s,&ghstart);dCHK(err);
      iMesh_rmvEntArrFromSet(mi,ents,ghstart,fs->set.explicit,&ierr);dICHK(mi,ierr);
      iMesh_addEntArrToSet(mi,ents,ghstart,fs->set.dirichlet,&ierr);dICHK(mi,ierr);
    }
    if (bstat[i] & dFSBSTATUS_WEAK) {
      iMesh_getEntitiesRec(mi,bdysets[i],dTYPE_FACE,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
      iMesh_addEntArrToSet(mi,ents,ents_s,fs->set.weakFace,&ierr);dICHK(mi,ierr);
    }
  }
  err = dFree2(bdysets,bstat);dCHK(err);
  free(ents);
  dFunctionReturn(0);
}

/* Set offsets (global, closure, local) of first node associated with every entity.
 *
 * @param indexTag index into inodes
 * @param inodes number of interior nodes associated with each entity inodes[idx[i]]
 * @param rstart relative start for global explicit dofs
 * @param crstart relative start for global closure dofs
 * @param nents size of array to hold entities
 * @param ents returns with entities in the proper ordering of the FS
 * @param ghstart returns starting offset of ghost entities in ents
 */
dErr dFSBuildSpaceOffsets_Private(dFS fs,dMeshTag indexTag,const dInt inodes[],dInt rstart,dInt crstart,dInt nents,dMeshEH ents[],dInt *ghstart)
{
  dInt i,scan,nentsExplicit,nentsDirichlet,*idx,*offset;
  dMesh mesh;
  dErr err;

  dFunctionBegin;
  err = dFSGetMesh(fs,&mesh);dCHK(err);
  err = dMallocA2(nents,&offset,nents,&idx);dCHK(err);
  /* We assume that orderedSet contains explicitSet+dirichletSet+ghostSet (in that order) */
  err = dMeshGetEnts(mesh,fs->set.ordered,dTYPE_ALL,dTOPO_ALL,ents,nents,NULL);dCHK(err);
  err = dMeshGetNumEnts(mesh,fs->set.explicit,dTYPE_ALL,dTOPO_ALL,&nentsExplicit);dCHK(err);
  err = dMeshGetNumEnts(mesh,fs->set.dirichlet,dTYPE_ALL,dTOPO_ALL,&nentsDirichlet);dCHK(err);
  err = dMeshTagGetData(mesh,indexTag,ents,nents,idx,nents,dDATA_INT);dCHK(err);

  /* global offset */
  for (i=0,scan=rstart; i<nentsExplicit; scan+=inodes[idx[i++]])
    offset[i] = scan;
  for ( ; i<nents; i++) offset[i] = -1;
  err = dMeshTagSetData(mesh,fs->tag.goffset,ents,nents,offset,nents,dDATA_INT);dCHK(err);

  /* global closure offset */
  for (i=0,scan=crstart; i<nentsExplicit+nentsDirichlet; scan+=inodes[idx[i++]])
    offset[i] = scan;
  for ( ; i<nents; i++) offset[i] = -1;
  err = dMeshTagSetData(mesh,fs->tag.gcoffset,ents,nents,offset,nents,dDATA_INT);dCHK(err);

  /* local index */
  for (i=0,scan=0; i<nents; scan+=inodes[idx[i++]])
    offset[i] = scan;
  err = dMeshTagSetData(mesh,fs->tag.loffset,ents,nents,offset,nents,dDATA_INT);dCHK(err);

  /* communicate global and closure offset for ghosts */
  err = dMeshTagBcast(mesh,fs->tag.goffset);dCHK(err);
  err = dMeshTagBcast(mesh,fs->tag.gcoffset);dCHK(err);

  err = dFree2(offset,idx);dCHK(err);
  *ghstart = nentsExplicit + nentsDirichlet;
  dFunctionReturn(0);
}

dErr dFSBuildSpaceVectors_Private(dFS fs,dMeshTag indexTag,const dInt inodes[],dInt rstart,dInt ghents_s,const dMeshEH ghents[])
{
  dErr err;
  dInt *gcoffset,*ghidx,*idx;
  dMesh mesh;

  dFunctionBegin;
  err = dMallocA3(ghents_s,&gcoffset,ghents_s,&ghidx,ghents_s,&idx);dCHK(err);
  err = dFSGetMesh(fs,&mesh);dCHK(err);
  err = dMeshTagGetData(mesh,indexTag,ghents,ghents_s,idx,ghents_s,dDATA_INT);dCHK(err);

  /* Retrieve ghost offsets, to create localupdate. */
  err = dMeshTagGetData(mesh,fs->tag.gcoffset,ghents,ghents_s,gcoffset,ghents_s,dDATA_INT);dCHK(err);
  for (dInt i=0; i<ghents_s; i++) { /* Paranoia: confirm that all ghost entities were updated. */
    if (gcoffset[i] < 0) dERROR(PETSC_COMM_SELF,1,"Tag exchange did not work");
  }

  /* Set ghost indices of every node using \a ghidx, create global vector. */
  {
    dInt gh=0;
    for (dInt i=0; i<ghents_s; i++) {
      for (dInt j=0; j<inodes[idx[i]]; j++) ghidx[gh++] = gcoffset[i] + j;
    }
    if (gh != fs->ngh) dERROR(PETSC_COMM_SELF,1,"Ghost count inconsistent");
    err = VecCreateDohp(((dObject)fs)->comm,fs->bs,fs->n,fs->nc,fs->ngh,ghidx,&fs->gvec);dCHK(err);
  }
  err = dFree3(gcoffset,ghidx,idx);dCHK(err);

  /* Create fs->bmapping and fs->mapping */
  err = dFSCreateLocalToGlobal_Private(fs,fs->n,fs->nc,fs->ngh,ghidx,rstart);dCHK(err);

  err = VecDohpCreateDirichletCache(fs->gvec,&fs->dcache,&fs->dscat);dCHK(err);
  dFunctionReturn(0);dCHK(err);
}
