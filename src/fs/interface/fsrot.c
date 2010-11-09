#include <dohpfsimpl.h>
#include <dohpvec.h>

/** Set rotation for certain nodes in a function space
*
* @param fs the function space
* @param is index set of nodes to rotate, sequential, with respect blocks of local vector
* @param rot Rotation matrices at all nodes in \a is.  Should have length \c bs*bs*size(is).
* @param ns number of dofs to enforce strongly at each node (every entry must have 0<=ns[i]<=bs)
* @param v Vector of values for strongly enforced dofs
*
* @example Consider 2D flow over a bed with known melt rates.  Suppose the local velocity vector is
*
*   [u0x,u0y; u1x,u1y; u2x,u2y; u3x,u3y | u4x,u4y]
*
* (4 owned blocks, one ghosted block) and nodes 1,4 are on the slip boundary with normal and tangent vectors n1,t1,n4,t4
* and melt rates r1,r4.  To enforce the melt rate strongly, use
*
* \a is = [1,4]
* \a rot = [n10,n11,t10,t11, n40,n41,t40,t41]
* \a ns = [1,1]
* \a v = [r1,r4]
*
* The rotated vector will become (. = \cdot)
*
*   [u0x,u0y; u1.n1,u1.t1; u2x,u2y; u3x,u3y | u4.n4,u4.t4]
*
* and strongly enforcing melt rate produces the global vector
*
*   [u0x,u0y; r1,u1.t1; u2x,u2y; u3x,u3y | r4,u4.t4] .
*
* This is what the solver sees, the Jacobian will always have rows and columns of the identity corresponding to the
* strongly enforced components (2,8 of the local vector) and the residual will always be 0 in these components.  Hence
* the Newton step v will always be of the form
*
*   [v0x,v0y; 0,v1y; v2x,v2y; v3x,v3y | 0,v4y] .
**/
dErr dFSRotationCreate(dFS fs,IS is,dReal rmat[],dInt ns[],Vec v,dFSRotation *inrot)
{
  dFSRotation rot;
  dInt bs,n;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_CLASSID,1);
  dValidHeader(is,IS_CLASSID,2);
  dValidRealPointer(rmat,3);
  dValidIntPointer(ns,4);
  dValidHeader(v,VEC_CLASSID,5);
  dValidPointer(inrot,6);
  *inrot = 0;
  err = PetscHeaderCreate(rot,_p_dFSRotation,struct _dFSRotationOps,dFSROT_CLASSID,0,"dFSRotation",PETSC_COMM_SELF,dFSRotationDestroy,dFSRotationView);dCHK(err);

  bs = rot->bs = fs->bs;
  err = ISGetSize(is,&n);dCHK(err);
  rot->n = n;
  err = PetscObjectReference((PetscObject)is);dCHK(err);
  rot->is = is;
  err = PetscObjectReference((PetscObject)v);dCHK(err);
  rot->strong = v;
  for (dInt i=0; i<n; i++) {
    if (ns[i] < 0 || bs < ns[i]) dERROR(PETSC_COMM_SELF,1,"Number of strong dofs must be between 0 and bs=%d (inclusive)",bs);
    /* \todo Check that every rmat is orthogonal */
  }
  err = dMallocA2(n*bs*bs,&rot->rmat,n,&rot->nstrong);dCHK(err);
  err = dMemcpy(rot->rmat,rmat,n*bs*bs*sizeof rmat[0]);dCHK(err);
  err = dMemcpy(rot->nstrong,ns,n*sizeof ns[0]);dCHK(err);
  *inrot = rot;
  dFunctionReturn(0);
}

dErr dFSRotationDestroy(dFSRotation rot)
{
  dErr err;

  dFunctionBegin;
  err = ISDestroy(rot->is);dCHK(err);
  err = VecDestroy(rot->strong);dCHK(err);
  err = dFree2(rot->rmat,rot->nstrong);dCHK(err);
  err = PetscHeaderDestroy(rot);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSRotationView(dFSRotation rot,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(rot,dFSROT_CLASSID,1);
  if (!viewer) {
    err = PetscViewerASCIIGetStdout(((dObject)rot)->comm,&viewer);dCHK(err);
  }
  dValidHeader(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(rot,1,viewer,2);

  dERROR(PETSC_COMM_SELF,1,"not implemented");
  dFunctionReturn(0);
}

/** Apply rotation to global vector
*
* @note does nothing if rotation is NULL
**/
dErr dFSRotationApply(dFSRotation rot,Vec g,dFSRotateMode rmode,dFSHomogeneousMode hmode)
{
  dErr err;
  Vec  gc,lf;

  dFunctionBegin;
  if (!rot) dFunctionReturn(0);
  dValidHeader(rot,dFSROT_CLASSID,1);
  dValidHeader(g,VEC_CLASSID,2);
  if (rot->ops->apply) {
    err = rot->ops->apply(rot,g,rmode,hmode);dCHK(err);
  } else {
    err = VecDohpGetClosure(g,&gc);dCHK(err);
    err = VecGhostGetLocalForm(gc,&lf);dCHK(err);
    err = dFSRotationApplyLocal(rot,lf,rmode,hmode);dCHK(err); /* \todo only rotate the global portion */
    err = VecGhostRestoreLocalForm(gc,&lf);dCHK(err);
    err = VecDohpRestoreClosure(g,&gc);dCHK(err);
  }
  dFunctionReturn(0);
}

/** Apply rotation to local vector
*
* @note does nothing if rotation is NULL
**/
dErr dFSRotationApplyLocal(dFSRotation rot,Vec l,dFSRotateMode rmode,dFSHomogeneousMode hmode)
{
  dErr err;

  dFunctionBegin;
  if (!rot) dFunctionReturn(0);
  dValidHeader(rot,dFSROT_CLASSID,1);
  dValidHeader(l,VEC_CLASSID,2);
  if (rot->ops->applylocal) {
    err = rot->ops->applylocal(rot,l,rmode,hmode);dCHK(err);
  } else {
    const dInt *idx,bs = rot->bs;
    const dScalar *s;
    dScalar *v,*strong;
    dScalar tmp[16];

    if (bs > 16) dERROR(PETSC_COMM_SELF,1,"large block size");
    err = ISGetIndices(rot->is,&idx);dCHK(err);
    err = VecGetArray(l,&v);dCHK(err);
    if (hmode == dFS_INHOMOGENEOUS) {
      err = VecGetArray(rot->strong,&strong);dCHK(err);
      s = strong;
    } else { s = 0; }
    for (dInt i=0; i<rot->n; i++) {
      const dInt ii = idx[i];
      const dReal *rmat = &rot->rmat[ii*bs*bs];
      for (dInt j=0; j<bs; j++) {
        tmp[j] = v[ii*bs+j];
        v[ii*bs+j] = 0;
      }
      switch (rmode) {
        case dFS_ROTATE_FORWARD:
          for (dInt j=0; j<bs; j++) {
            for (dInt k=0; k<bs; k++) {
              v[ii*bs+k] += rmat[k*bs+j] * tmp[j];
            }
          }
          switch (hmode) {
            case dFS_HOMOGENEOUS:   for (dInt j=0; j<s[i]; j++) v[ii*bs+j] = 0;    break;
            case dFS_INHOMOGENEOUS: for (dInt j=0; j<s[i]; j++) v[ii*bs+j] = *s++; break;
            default: dERROR(PETSC_COMM_SELF,1,"Invalid homogeneous mode");
          }
          break;
        case dFS_ROTATE_REVERSE:
          switch (hmode) {
            case dFS_HOMOGENEOUS:   for (dInt j=0; j<s[i]; j++) tmp[j] = 0;    break;
            case dFS_INHOMOGENEOUS: for (dInt j=0; j<s[i]; j++) tmp[j] = *s++; break;
            default: dERROR(PETSC_COMM_SELF,1,"Invalid homogeneous mode");
          }
          for (dInt j=0; j<bs; j++) {
            for (dInt k=0; k<bs; k++) {
              v[ii*bs+k] += rmat[j*bs+k] * tmp[k];
            }
          }
          break;
        default: dERROR(PETSC_COMM_SELF,1,"Invalid rotate mode");
      }
    }
    err = ISRestoreIndices(rot->is,&idx);dCHK(err);
    err = VecRestoreArray(l,&v);dCHK(err);
    if (hmode == dFS_HOMOGENEOUS) {
      dInt n;
      err = VecGetSize(rot->strong,&n);dCHK(err);
      if (s-strong != n) dERROR(PETSC_COMM_SELF,1,"should not happen");
      err = VecRestoreArray(rot->strong,&strong);dCHK(err);
    }
  }
  dFunctionReturn(0);
}
