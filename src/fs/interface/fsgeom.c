#include <dohpfsimpl.h>
#include <dohpvec.h>

/** Get coordinates for every node in closure (every subelement vertex)
*
* @param fs Function space
* @param inx the new vector with block size 3 and the same number of blocks as the closure vector
**/
dErr dFSGetCoordinates(dFS fs,Vec *inx)
{
  dErr    err;
  dInt    nelems;
  Vec     Expanded,Geom,X,Count;
  dScalar *x;
  const dScalar *geom;
  dFS     cfs;
  const dEFS *cefs;
  dRuleSet ruleset;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(inx,2);
  *inx = 0;

  err = dFSGetCoordinateFS(fs,&cfs);dCHK(err);
  err = dFSGetGeometryVectorExpanded(fs,&Geom);dCHK(err);
  err = VecDuplicate(Geom,&Expanded);dCHK(err);
  err = dFSCreateGlobalVector(cfs,&X);dCHK(err);
  err = VecDuplicate(X,&Count);dCHK(err);

  /* Count the number of occurances of each node in the closure. */
  err = VecSet(Expanded,1.);dCHK(err);
  err = VecZeroEntries(Count);dCHK(err);
  err = dFSExpandedToGlobal(cfs,Expanded,Count,dFS_INHOMOGENEOUS,ADD_VALUES);dCHK(err);

  /* Evaluate the coordinate basis functions on the interpolation nodes of the given function space. */
  err = dFSGetPreferredQuadratureRuleSet(fs,fs->set.active,dTYPE_REGION,dTOPO_ALL,dQUADRATURE_METHOD_SELF,&ruleset);dCHK(err);
  err = dFSGetEFS(cfs,ruleset,&nelems,&cefs);dCHK(err);
  err = VecGetArrayRead(Geom,&geom);dCHK(err);
  err = VecGetArray(Expanded,&x);dCHK(err);
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
  err = VecRestoreArray(Expanded,&x);dCHK(err);

  err = VecZeroEntries(X);dCHK(err);
  err = dFSExpandedToGlobal(cfs,Expanded,X,dFS_INHOMOGENEOUS,ADD_VALUES);dCHK(err);
  err = VecPointwiseDivide(X,X,Count);dCHK(err);

  err = VecDestroy(Expanded);dCHK(err);
  err = VecDestroy(Count);dCHK(err);
  err = VecDestroy(Geom);dCHK(err);
  *inx = X;
  dFunctionReturn(0);
}

static dErr dFSCreateGeometryFromMesh_Private(dFS fs)
{
  dErr err;
  dMesh mesh;
  dJacobi jacobi;
  dMeshEH *ents;
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
  err = dMallocA(nents,&ents);dCHK(err);
  err = dMeshGetEnts(mesh,fs->set.active,dTYPE_REGION,dTOPO_ALL,ents,nents,NULL);dCHK(err);
  err = dMeshGetAdjVertexCoords(mesh,nents,ents,&offset,&rcoords);dCHK(err);
  err = VecGetLocalSize(Expanded,&n);dCHK(err);
  if (n != 3*offset[nents]) dERROR(PETSC_ERR_PLIB,"Size of Expanded %D does not match offset %D",n,3*offset[nents]);
  err = VecGetArray(Expanded,&x);dCHK(err);
  for (dInt i=0; i<n; i++) x[i] = rcoords[i];
  err = VecRestoreArray(Expanded,&x);dCHK(err);
  err = dMeshRestoreAdjVertexCoords(mesh,nents,ents,&offset,&rcoords);dCHK(err);
  err = dFree(ents);dCHK(err);

  /* Populate \a Global with the average coordinates from all the elements whose closure contains the node */
  err = dFSExpandedToGlobal(cfs,Expanded,Global,dFS_INHOMOGENEOUS,ADD_VALUES);dCHK(err);
  err = VecPointwiseDivide(Global,Global,Count);dCHK(err);
  err = VecDestroy(Count);dCHK(err);

  fs->geometry.expanded = Expanded;
  fs->geometry.global   = Global;
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
  if (!fs->ops->getsubelementmeshsize) dERROR(1,"not implemented");
  err = (*fs->ops->getsubelementmeshsize)(fs,nelems,nverts,nconn);dCHK(err);
  dFunctionReturn(0);
}

/** Get the number of subelements and number of vertices of subelements.
*
* @note the number of vertices is the same as the number of local nodes in closure vertor.
**/
dErr dFSGetSubElementMesh(dFS fs,dInt nelems,dInt nverts,dEntTopology topo[],dInt off[],dInt ind[])
{
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidIntPointer(topo,4);
  dValidIntPointer(off,5);
  dValidIntPointer(ind,6);
  if (!fs->ops->getsubelementmesh) dERROR(1,"not implemented");
  err = (*fs->ops->getsubelementmesh)(fs,nelems,nverts,topo,off,ind);dCHK(err);
  dFunctionReturn(0);
}
