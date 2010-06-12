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
  dInt    n,*off,*geomoff,nelems;
  Mat     E1,Ebs;
  Vec     xc,xx,xc1,xx1;
  dScalar *count,*x;
  dReal   (*qx)[3],(*geom)[3];
  dEFS    *efs;
  dRuleSet ruleset;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(inx,2);
  *inx = 0;

  err = MatMAIJRedimension(fs->E,1,&E1);dCHK(err);
  err = MatGetVecs(E1,&xc1,&xx1);dCHK(err);
  err = VecSet(xx1,1);dCHK(err);
  err = MatMultTranspose(E1,xx1,xc1);dCHK(err);
  err = MatMult(E1,xc1,xx1);dCHK(err); /* xx1 now has the number of times each node occurs in the expanded vector */

  err = MatMAIJRedimension(E1,3,&Ebs);dCHK(err);
  err = MatGetVecs(Ebs,&xc,&xx);dCHK(err);

  err = dFSGetPreferredQuadratureRuleSet(fs,fs->set.active,dTYPE_REGION,dTOPO_ALL,dQUADRATURE_METHOD_FAST,&ruleset);dCHK(err);
  err = dFSGetEFS(fs,ruleset,&nelems,&efs);dCHK(err);
  dERROR(1,"In flux");
#if 0
  err = dFSGetElements(fs,&n,&off,NULL,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&qx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err);
  err = VecGetArray(xx1,&count);dCHK(err);
  err = VecGetArray(xx,&x);dCHK(err);
  for (dInt e=0; e<n; e++) {
    /* Offsets are given in terms of blocks */
    const dScalar *counte = count + off[e];
    dScalar       *xe     = x + 3*off[e];
    dInt           three,P[3];
    err = dEFSGetGlobalCoordinates(efs[e],(const dReal(*)[3])(geom+geomoff[e]),&three,P,qx);dCHK(err);
    if (three != 3) dERROR(1,"element dimension unexpected");
    for (dInt i=0; i<P[0]*P[1]*P[2]; i++) {
      for (dInt j=0; j<3; j++) xe[i*3+j] = qx[i][j] / counte[i];
    }
  }
  err = dFSRestoreElements(fs,&n,&off,NULL,&efs,&geomoff,&geom);dCHK(err);
  err = dFSRestoreWorkspace(fs,__func__,&qx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err);
#endif
  err = VecRestoreArray(xx1,&count);dCHK(err);
  err = VecRestoreArray(xx,&x);dCHK(err);

  err = MatMultTranspose(Ebs,xx,xc);dCHK(err);

  err = VecDestroy(xx1);dCHK(err);
  err = VecDestroy(xc1);dCHK(err);
  err = VecDestroy(xx);dCHK(err);
  err = MatDestroy(Ebs);dCHK(err);
  err = MatDestroy(E1);dCHK(err);
  *inx = xc;
  dFunctionReturn(0);
}

/** Get displacements of every point in geometry vector.
*
* @param fs Function space
* @param ingeom the new geometry vector with block size 3 and same number of blocks as the closure vector
*
* @note I don't yet have a clear view of how this should work.  If we have mesh motion with interesting boundary
* conditions, we probably need a bonified function space for the geometry.
*
* @note It is important to get geometry associated with boundary sets because it will frequently be projected against
* the boundary.
**/
dErr dFSGetGeometryVector(dFS dUNUSED fs,Vec *ingeom)
{

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidPointer(ingeom,2);
  *ingeom = 0;
  dERROR(1,"not implemented");
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
