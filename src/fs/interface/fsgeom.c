#include <dohpfsimpl.h>
#include <dohpvec.h>

/** dFSGetNodalCoordinateFS - Gets an FS the same size as the solution FS, but with block size 3 to hold nodal coordinates
 *
 */
dErr dFSGetNodalCoordinateFS(dFS fs,dFS *nfs)
{
  dErr err;

  dFunctionBegin;
  if (!fs->nodalcoord.fs) {err = dFSRedimension(fs,3,dFS_INTERIOR,&fs->nodalcoord.fs);dCHK(err);}
  *nfs = fs->nodalcoord.fs;
  dFunctionReturn(0);
}

dErr dFSGetNodalCoordinatesExpanded(dFS fs,Vec *inX)
{
  dErr err;
  dFS cfs;
  Vec Geom,Expanded3;
  const dScalar *geom;
  dScalar *x;
  const dEFS *cefs;
  dRuleset ruleset;
  dInt nelems;

  dFunctionBegin;
  if (!fs->nodalcoord.expanded) {
    dFS fs3;
    err = dFSGetNodalCoordinateFS(fs,&fs3);dCHK(err);
    err = dFSCreateExpandedVector(fs3,&fs->nodalcoord.expanded);dCHK(err);
  }
  Expanded3 = fs->nodalcoord.expanded;

  err = dFSGetCoordinateFS(fs,&cfs);dCHK(err);
  err = dFSGetGeometryVectorExpanded(fs,&Geom);dCHK(err);

  /* Evaluate the coordinate basis functions on the interpolation nodes of the given function space. */
  err = dFSGetPreferredQuadratureRuleSet(fs,fs->set.active,dTYPE_REGION,dTOPO_ALL,dQUADRATURE_METHOD_SELF,&ruleset);dCHK(err);
  err = dFSGetEFS(cfs,ruleset,&nelems,&cefs);dCHK(err);
  err = VecGetArrayRead(Geom,&geom);dCHK(err);
  err = VecGetArray(Expanded3,&x);dCHK(err);
  for (dInt e=0,goffset=0,offset=0; e<nelems; e++) {
    dInt step,nnodes;
    dRule rule;
    err = dEFSApply(cefs[e],NULL,3,geom+goffset,x+offset,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSGetSizes(cefs[e],NULL,NULL,&step);dCHK(err);
    err = dEFSGetRule(cefs[e],&rule);dCHK(err);
    err = dRuleGetSize(rule,NULL,&nnodes);dCHK(err);
    goffset += step*3;
    offset  += nnodes*3;
  }
  err = VecRestoreArrayRead(Geom,&geom);dCHK(err);
  err = VecRestoreArray(Expanded3,&x);dCHK(err);
  err = dFSRestoreEFS(cfs,ruleset,&nelems,&cefs);dCHK(err);

  err = dRulesetDestroy(&ruleset);dCHK(err);
  *inX = Expanded3;
  dFunctionReturn(0);
}

/** Get coordinates for every node in closure (every subelement vertex)
 *
 * @note This function cannot be implemented for all \a dFS types.  For most purposes, users should
 *       \a dFSGetGeometryVectorExpanded and evaluate (element by element) on the nodes of their choice
 *       (with a self-quadrature).
 *
 * @param fs Function space
 * @param inx the new vector with block size 3 and the same number of blocks as the closure vector
 **/
dErr dFSGetNodalCoordinatesGlobal(dFS fs,Vec *inx)
{
  dErr    err;
  Vec     Expanded,Ones,X,Count,Xclosure,Countclosure;
  dFS     fs3;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(inx,2);
  *inx = 0;

  err = dFSGetNodalCoordinateFS(fs,&fs3);dCHK(err);
  err = dFSGetNodalCoordinatesExpanded(fs,&Expanded);dCHK(err);
  if (!fs->nodalcoord.global) {err = dFSCreateGlobalVector(fs3,&fs->nodalcoord.global);dCHK(err);}
  X = fs->nodalcoord.global;

  /* Count the number of occurances of each node in the closure. */
  err = VecDuplicate(Expanded,&Ones);dCHK(err);
  err = VecDuplicate(X,&Count);dCHK(err);

  err = VecDohpZeroEntries(Count);dCHK(err);
  err = VecSet(Ones,1.);dCHK(err);
  err = dFSExpandedToGlobal(fs3,Ones,Count,dFS_INHOMOGENEOUS,ADD_VALUES);dCHK(err);
  err = VecDestroy(&Ones);dCHK(err);

  err = VecDohpZeroEntries(X);dCHK(err);
  err = dFSExpandedToGlobal(fs3,Expanded,X,dFS_INHOMOGENEOUS,ADD_VALUES);dCHK(err);

  err = VecDohpGetClosure(X,&Xclosure);dCHK(err);
  err = VecDohpGetClosure(Count,&Countclosure);dCHK(err);
  err = VecPointwiseDivide(Xclosure,Xclosure,Countclosure);dCHK(err);
  err = VecDohpRestoreClosure(X,&Xclosure);dCHK(err);
  err = VecDohpRestoreClosure(Count,&Countclosure);dCHK(err);

  err = VecDestroy(&Count);dCHK(err);
  *inx = X;
  dFunctionReturn(0);
}

static dErr dFSCreateGeometryFromMesh_Private(dFS fs)
{
  dErr err;
  dMesh mesh;
  dJacobi jacobi;
  dMeshEH *ents;
  dEntTopology *topo;
  const dInt *offset;
  const dReal *rcoords;
  dScalar *x;
  dInt nents,n;
  dFS cfs;
  Vec Expanded,Global,Count;

  dFunctionBegin;
  err = dFSGetMesh(fs,&mesh);dCHK(err);
  err = dFSGetJacobi(fs,&jacobi);dCHK(err);
  err = dFSGetCoordinateFS(fs,&cfs);dCHK(err);
  err = dFSCreateExpandedVector(cfs,&Expanded);dCHK(err);
  err = dFSCreateGlobalVector(cfs,&Global);dCHK(err);
  err = VecDuplicate(Global,&Count);dCHK(err);

  /* Count the number of elements in which each node appears */
  err = VecSet(Expanded,1.);dCHK(err);
  err = VecZeroEntries(Count);dCHK(err);
  err = dFSExpandedToGlobal(cfs,Expanded,Count,dFS_HOMOGENEOUS,ADD_VALUES);dCHK(err);

  /* Populate \a Expanded with local coordinates from mesh */
  err = dMeshGetNumEnts(mesh,fs->set.active,dTYPE_REGION,dTOPO_ALL,&nents);dCHK(err);
  err = dMallocA2(nents,&ents,nents,&topo);dCHK(err);
  err = dMeshGetEnts(mesh,fs->set.active,dTYPE_REGION,dTOPO_ALL,ents,nents,NULL);dCHK(err);
  err = dMeshGetTopo(mesh,nents,ents,topo);dCHK(err);
  err = dMeshGetAdjVertexCoords(mesh,nents,ents,&offset,&rcoords);dCHK(err);
  err = VecGetLocalSize(Expanded,&n);dCHK(err);
  if (n != 3*offset[nents]) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Size of Expanded %D does not match offset %D",n,3*offset[nents]);
  err = VecGetArray(Expanded,&x);dCHK(err);
  for (dInt e=0,k=0; e<nents; e++) {
    switch (topo[e]) {
      case dTOPO_HEX:
        for (dInt i=0; i<8; i++) {
          static const dInt connidx[8] = {0,4,3,7,1,5,2,6}; /* FIXME: generalize this */
          for (dInt j=0; j<3; j++) x[k++] = rcoords[(offset[e]+connidx[i])*3+j];
        }
        break;
      default: dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"No support for topology %s",dMeshEntTopologyName(topo[e]));
    }
  }
  err = VecRestoreArray(Expanded,&x);dCHK(err);
  err = dMeshRestoreAdjVertexCoords(mesh,nents,ents,&offset,&rcoords);dCHK(err);
  err = dFree2(ents,topo);dCHK(err);

  /* Populate \a Global with the average coordinates from all the elements whose closure contains the node */
  err = dFSExpandedToGlobal(cfs,Expanded,Global,dFS_INHOMOGENEOUS,ADD_VALUES);dCHK(err);
  err = VecPointwiseDivide(Global,Global,Count);dCHK(err);
  err = VecDestroy(&Count);dCHK(err);

  fs->geometry.expanded = Expanded;
  fs->geometry.global   = Global;
  dFunctionReturn(0);
}

/* Get a parallel dFS for representing the geometry
 *
 * @note Creates a new one if not yet set, generally based on the nodal locations in the mesh
 */
dErr dFSGetCoordinateFS(dFS fs,dFS *incfs)
{
  dErr err;

  dFunctionBegin;
  if (!fs->geometry.fs) {
    dFS cfs;
    dMesh mesh;
    dJacobi jacobi;
    dMeshTag dtag;
    char degreename[256],prefix[256];

    err = dFSGetMesh(fs,&mesh);dCHK(err);
    err = dFSGetJacobi(fs,&jacobi);dCHK(err);
    err = dFSCreate(((dObject)fs)->comm,&cfs);dCHK(err);
    err = dFSSetMesh(cfs,mesh,fs->set.active);dCHK(err);
    err = dFSSetOrderingType(cfs,fs->orderingtype);dCHK(err);
    err = PetscSNPrintf(degreename,sizeof degreename,"%scfs_degree",dNonNullElse(((dObject)fs)->prefix,""));dCHK(err);
    err = dMeshCreateRuleTagIsotropic(mesh,fs->set.active,degreename,1,&dtag);dCHK(err);
    err = dFSSetDegree(cfs,jacobi,dtag);dCHK(err);
    err = dFSSetBlockSize(cfs,3);dCHK(err);
    err = PetscSNPrintf(prefix,sizeof prefix,"%scoord_",dNonNullElse(((dObject)fs)->prefix,""));dCHK(err);
    err = dFSSetOptionsPrefix(cfs,prefix);dCHK(err);
    err = dFSSetFromOptions(cfs);dCHK(err);
    fs->geometry.fs = cfs;
  }
  *incfs = fs->geometry.fs;
  dFunctionReturn(0);
}

/** Get displacements of every node in expanded vector.
*
* @param fs Function space
* @param ingeom the new geometry vector with block size 3 and same number of blocks as points in the expanded vector.
*
* @note This is not suitable for use as a function space, 
*
* @note It is important to get geometry associated with boundary sets because it will frequently be projected against
* the boundary.
**/
dErr dFSGetGeometryVectorExpanded(dFS fs,Vec *ingeom)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(ingeom,2);
  *ingeom = 0;
  if (!fs->geometry.expanded) {err = dFSCreateGeometryFromMesh_Private(fs);dCHK(err);}
  *ingeom = fs->geometry.expanded;
  dFunctionReturn(0);
}

dErr dFSGetBoundingBox(dFS fs,dReal bbox[3][2])
{
  dErr err;
  Vec X;
  const dScalar *x;
  dInt n;

  dFunctionBegin;
  for (dInt i=0; i<3; i++) {
    bbox[i][0] = PETSC_MAX_REAL;
    bbox[i][1] = PETSC_MIN_REAL;
  }
  err = dFSGetGeometryVectorExpanded(fs,&X);dCHK(err);
  err = VecGetLocalSize(X,&n);dCHK(err);
  err = VecGetArrayRead(X,&x);dCHK(err);
  for (dInt i=0; i<n; i++) {
    dInt j = i%3;
    bbox[j][0] = dMin(bbox[j][0],x[i]);
    bbox[j][1] = dMax(bbox[j][1],x[i]);
  }
  err = VecRestoreArrayRead(X,&x);dCHK(err);
  dFunctionReturn(0);
}

/** Get the number of subelements and number of vertices of subelements.
*
* @note the number of vertices is the same as the number of local nodes in closure vertor.
**/
dErr dFSGetSubElementMeshSize(dFS fs,dInt *nelems,dInt *nverts,dInt *nconn)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidIntPointer(nelems,2);
  dValidIntPointer(nverts,3);
  dValidIntPointer(nconn,4);
  if (!fs->ops->getsubelementmeshsize) dERROR(PETSC_COMM_SELF,1,"not implemented");
  err = (*fs->ops->getsubelementmeshsize)(fs,nelems,nverts,nconn);dCHK(err);
  dFunctionReturn(0);
}

/** Get the number of subelements and number of vertices of subelements.
**/
dErr dFSGetSubElementMesh(dFS fs,dInt nelems,dInt nconn,dEntTopology topo[],dInt off[],dInt ind[])
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidIntPointer(topo,4);
  dValidIntPointer(off,5);
  dValidIntPointer(ind,6);
  if (!fs->ops->getsubelementmesh) dERROR(PETSC_COMM_SELF,1,"not implemented");
  err = (*fs->ops->getsubelementmesh)(fs,nelems,nconn,topo,off,ind);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSSubElementMeshView(dFS fs,dViewer view)
{
  dErr err;
  dInt nelem,nverts,nconn,*suboff,*subind,n,bs;
  dEntTopology *subtopo;
  Vec X,Xc;
  const dReal *x;
  const char *name;

  dFunctionBegin;
  err = dFSGetSubElementMeshSize(fs,&nelem,&nverts,&nconn);dCHK(err);
  dASSERT(nconn == 8*nelem);    /* Assume Hex */
  err = dMallocA3(nelem,&subtopo,nelem+1,&suboff,nconn,&subind);dCHK(err);
  err = dFSGetSubElementMesh(fs,nelem,nconn,subtopo,suboff,subind);dCHK(err);
  err = PetscObjectGetName((dObject)fs,&name);dCHK(err);
  err = PetscViewerASCIIPrintf(view,"SubElementMesh name=%s nelem=%D nverts=%D\n",name,nelem,nverts);dCHK(err);
  err = dIntTableView(nelem,8,subind,view,"subconn");dCHK(err);
  err = dFree3(subtopo,suboff,subind);dCHK(err);

  err = dFSGetNodalCoordinatesGlobal(fs,&X);dCHK(err);
  err = VecDohpGetClosure(X,&Xc);dCHK(err);
  err = VecGetLocalSize(Xc,&n);dCHK(err);
  err = VecGetBlockSize(Xc,&bs);dCHK(err);
  err = VecGetArrayRead(Xc,&x);dCHK(err);
  err = dRealTableView(n/bs,bs,x,view,"coords");dCHK(err);
  err = VecRestoreArrayRead(Xc,&x);dCHK(err);
  err = VecDohpRestoreClosure(X,&Xc);dCHK(err);
  dFunctionReturn(0);
}

/** dFSRedimension - Create a new function space with the same layout, but different number of dofs per node
 *
 */
dErr dFSRedimension(dFS fs,dInt bs,dFSClosureMode mode,dFS *infs)
{
  dErr err;
  dFS rfs;
  dJacobi jac;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  if (bs < 1) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Block size must be at least 1, was %D",bs);
  dValidPointer(infs,4);
  *infs = NULL;

  if (!fs->spacebuilt) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Cannot redimension a space that has not been built");
  if (!fs->mesh || !fs->set.active) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Space has been built, but does not have a mesh (should not happen)");

  err = dFSCreate(((dObject)fs)->comm,&rfs);dCHK(err);
  err = dFSSetBlockSize(rfs,bs);dCHK(err);
  err = dFSSetMesh(rfs,fs->mesh,fs->set.active);dCHK(err);
  err = dFSGetJacobi(fs,&jac);dCHK(err);
  err = dFSSetDegree(rfs,jac,fs->tag.degree);dCHK(err);
  switch (mode) {
  case dFS_CLOSURE: dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"Probably not what you want"); /* because ordering would be different, but there is nothing to do for this choice. */
    break;
  case dFS_INTERIOR: {
    dMeshESH *bsets;
    dInt nsets;
    err = dMeshGetNumSubsets(fs->mesh,fs->set.boundaries,0,&nsets);dCHK(err);
    err = dMallocA(nsets,&bsets);dCHK(err);
    err = dMeshGetSubsets(fs->mesh,fs->set.boundaries,0,bsets,nsets,NULL);dCHK(err);
    for (dInt i=0; i<nsets; i++) {
      dFSBStatus bstat;
      err = dMeshTagSGetData(fs->mesh,fs->tag.bstatus,&bsets[i],1,&bstat,sizeof(bstat),dDATA_BYTE);dCHK(err);
      err = dFSRegisterBoundarySet(rfs,bsets[i],bstat,NULL,NULL);dCHK(err); /* TODO: handle dFSConstraintCtx */
    }
    err = dFree(bsets);dCHK(err);
  } break;
  default: dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"No support for dFSClosureMode %s",dFSClosureModes[mode]);
  }
  err = dFSSetType(rfs,((dObject)fs)->type_name);dCHK(err);
  err = dFSSetOrderingType(rfs,fs->orderingtype);dCHK(err);
  rfs->set.ordered = fs->set.ordered;
  { // The FS has the layout, ordering, and boundary status tags set so we are ready to build the function space.
    dMesh          mesh;
    dMeshAdjacency meshadj;

    err = dFSGetMesh(rfs,&mesh);dCHK(err);
    err = dMeshGetAdjacency(mesh,rfs->set.ordered,&meshadj);dCHK(err);
    err = dFSPopulatePartitionedSets_Private(rfs,meshadj);dCHK(err);
    err = dFSBuildSpaceWithOrderedSet_Private(rfs,meshadj);dCHK(err);
    err = dMeshRestoreAdjacency(mesh,rfs->set.ordered,&meshadj);dCHK(err);
  }

  *infs = rfs;
  dFunctionReturn(0);
}
