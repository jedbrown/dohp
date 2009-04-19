#include "private/fsimpl.h"

/** Get coordinates for every node in closure
*
* @param fs Function space
* @param inx the new vector with block size 3 and the same number of blocks as the closure vector
**/
dErr dFSGetClosureCoordinates(dFS fs,Vec *inx)
{
  dErr err;
  dInt n,bs;
  const VecType type;
  Vec x;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidPointer(inx,2);
  *inx = 0;
  /* Poor man's duplicate with different block size */
  err = VecGetLocalSize(fs->gc,&n);dCHK(err);
  err = VecGetBlockSize(fs->gc,&bs);dCHK(err);
  err = VecGetType(fs->gc,&type);dCHK(err);
  err = VecCreate(((dObject)fs)->comm,&x);dCHK(err);
  err = VecSetSizes(x,3*n/bs,PETSC_DETERMINE);dCHK(err);
  err = VecSetBlockSize(x,3);dCHK(err);
  err = VecSetType(x,type);dCHK(err);
  dERROR(1,"not implemented");
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
dErr dFSGetGeometryVector(dFS fs,Vec *ingeom)
{

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  dValidPointer(ingeom,2);
  *ingeom = 0;
  dERROR(1,"not implemented");
  dFunctionReturn(0);
}
