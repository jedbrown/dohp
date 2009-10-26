#include <dohpfsimpl.h>
#include <dohpvec.h>

/** Get coordinates for every node in closure
*
* @param fs Function space
* @param inx the new vector with block size 3 and the same number of blocks as the closure vector
**/
dErr dFSGetCoordinates(dFS fs,Vec *inx)
{
  dErr          err;
  dInt          n,*off,*geomoff;
  Mat           E1,Ebs;
  Vec           xc,xx,xc1,xx1;
  dScalar      *count,*x;
  dReal       (*qx)[3],(*geom)[3];
  s_dEFS       *efs;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidPointer(inx,2);
  *inx = 0;

  err = MatMAIJRedimension(fs->E,1,&E1);dCHK(err);
  err = MatGetVecs(E1,&xc1,&xx1);dCHK(err);
  err = VecSet(xx1,1);dCHK(err);
  err = MatMultTranspose(E1,xx1,xc1);dCHK(err);
  err = MatMult(E1,xc1,xx1);dCHK(err); /* xx1 now has the number of times each node occurs in the expanded vector */

  err = MatMAIJRedimension(E1,3,&Ebs);dCHK(err);
  err = MatGetVecs(Ebs,&xc,&xx);dCHK(err);

  err = dFSGetElements(fs,&n,&off,NULL,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&qx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err);
  err = VecGetArray(xx1,&count);dCHK(err);
  err = VecGetArray(xx,&x);dCHK(err);
  for (dInt e=0; e<n; e++) {
    /* Offsets are given in terms of blocks */
    const dScalar *counte = count + off[e];
    dScalar       *xe     = x + 3*off[e];
    dInt           three,P[3];
    err = dEFSGetGlobalCoordinates(&efs[e],(const dReal(*)[3])(geom+geomoff[e]),&three,P,qx);dCHK(err);
    if (three != 3) dERROR(1,"element dimension unexpected");
    for (dInt i=0; i<P[0]*P[1]*P[2]; i++) {
      for (dInt j=0; j<3; j++) xe[i*3+j] = qx[i][j] / counte[i];
    }
  }
  err = dFSRestoreElements(fs,&n,&off,NULL,&efs,&geomoff,&geom);dCHK(err);
  err = dFSRestoreWorkspace(fs,__func__,&qx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err);
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
  dValidHeader(fs,DM_COOKIE,1);
  dValidPointer(ingeom,2);
  *ingeom = 0;
  dERROR(1,"not implemented");
  dFunctionReturn(0);
}
