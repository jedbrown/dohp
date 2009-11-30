#include <dohpfsimpl.h>
#include <dohpvec.h>

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
    if (a[i*bs] != 1) dERROR(1,"should not happen");
    globals[i] = rstart+i;
  }
  for ( ; i<nc; i++) {
    if (a[i*bs] != -1) dERROR(1,"should not happen");
    globals[i] = -(rstart + i);
  }
  for ( ; i<nc+ngh; i++) {
    globals[i] = signbit(a[i*bs]) * ghidx[i-nc];
  }
  err = VecRestoreArray(lf,&a);dCHK(err);
  err = VecGhostRestoreLocalForm(gc,&lf);dCHK(err);
  err = VecDohpRestoreClosure(g,&gc);dCHK(err);
  err = VecDestroy(g);dCHK(err);
  err = ISLocalToGlobalMappingCreateNC(((dObject)fs)->comm,nc+ngh,globals,&fs->bmapping);dCHK(err);
  /* Don't free \a globals because we used the no-copy variant, so the IS takes ownership. */
  {
    dInt *sglobals;        /* Scalar globals */
    err = dMallocA((nc+ngh)*bs,&sglobals);dCHK(err);
    for (i=0; i<(nc+ngh)*bs; i++) sglobals[i] = globals[i/bs]*bs + i%bs;
    err = ISLocalToGlobalMappingCreateNC(((dObject)fs)->comm,(nc+ngh)*bs,sglobals,&fs->mapping);dCHK(err);
    /* Don't free \a sglobals because we used the no-copy variant, so the IS takes ownership. */
  }
  dFunctionReturn(0);
}
