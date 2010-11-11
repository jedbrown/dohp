#include <dohpfsimpl.h>
#include <dohpvec.h>
#include <stdarg.h>

struct dRulesetIteratorLink {
  dFS fs;
  const dEFS *efs;
  Vec X,Xexp,Y,Yexp;
  dScalar *x;                   /* Entries from expanded vector */
  dScalar *y;
  dScalar *u,*du,*v,*dv;
  dInt nefs;
  dInt off;
  dInt bs;
  struct dRulesetIteratorLink *next;
};

struct dRulesetIteratorStash {
  char *patch;
  char *node;
  dInt *patchoffset;
  dInt patchbytes,nodebytes;
};

struct _n_dRulesetIterator {
  dRuleset ruleset;
  dInt curpatch;
  dInt npatches,nnodes;
  dInt maxQ;
  dInt nlinks;
  dInt Q;
  dScalar *cjinv,*jw;
  struct dRulesetIteratorStash stash;
  struct dRulesetIteratorLink *link;
};

/** Get an iterator for performing an integral on the given rule set, with coordinate vectors lying in space cfs.
 *
 * Collective on dFS
 */
dErr dRulesetCreateIterator(dRuleset rset,dFS cfs,dRulesetIterator *iter)
{
  dErr err;
  dRulesetIterator it;

  dFunctionBegin;
  *iter = 0;
  err = dCallocA(1,&it);dCHK(err);
  it->ruleset = rset;
  err = dRulesetGetSize(rset,&it->npatches);dCHK(err);
  it->nnodes = 0;
  for (dInt i=0; i<it->npatches; i++) {
    dInt nnodes;
    err = dRuleGetSize(it->ruleset->rules[i],NULL,&nnodes);dCHK(err);
    it->nnodes += nnodes;
  }
  err = dRulesetGetMaxQ(rset,&it->maxQ);dCHK(err);
  err = dRulesetIteratorAddFS(it,cfs);dCHK(err);
  *iter = it;
  dFunctionReturn(0);
}

/* Collective on dFS */
dErr dRulesetIteratorAddFS(dRulesetIterator it,dFS fs)
{
  struct dRulesetIteratorLink *link,*p;
  dErr err;

  dFunctionBegin;
  err = dCallocA(1,&link);dCHK(err);
  err = PetscObjectReference((PetscObject)fs);dCHK(err);
  if (link->fs) {err = PetscObjectDereference((PetscObject)fs);dCHK(err);}
  link->fs = fs;
  err = dFSGetBlockSize(fs,&link->bs);dCHK(err);
  if (!it->link) it->link = link;
  else {
    for (p = it->link; p->next; p=p->next) {}
    p->next = link;
  }
  it->nlinks++;
  dFunctionReturn(0);
}

static dErr dRulesetIteratorLinkCreatePatchSpace_Private(struct dRulesetIteratorLink *link,dInt maxQ)
{
  dErr err;
  dInt bs = link->bs;

  dFunctionBegin;
  if (link->u) dFunctionReturn(0);
  err = dMallocA4(maxQ*bs,&link->u,maxQ*bs*3,&link->du,maxQ*bs,&link->v,maxQ*bs*3,&link->dv);dCHK(err);
  dFunctionReturn(0);
}

static dErr dRulesetIteratorLinkFreePatchSpace_Private(struct dRulesetIteratorLink *link)
{
  dErr err;

  dFunctionBegin;
  err = dFree4(link->u,link->du,link->v,link->dv);dCHK(err);
  dFunctionReturn(0);
}


dErr dRulesetIteratorStart(dRulesetIterator it,Vec X,Vec Y,...)
{
  dErr err;
  va_list ap;
  dInt i;
  struct dRulesetIteratorLink *p;

  dFunctionBegin;
  va_start(ap,Y);
  for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
    if (i) {
      X = va_arg(ap,Vec);
      Y = va_arg(ap,Vec);
    }
    p->X = X;
    p->Y = Y;
    if (!p->Xexp) {err = dFSCreateExpandedVector(p->fs,&p->Xexp);dCHK(err);}
    if (!p->Yexp) {err = dFSCreateExpandedVector(p->fs,&p->Yexp);dCHK(err);}
    {
      dBool flg;
      err = PetscTypeCompare((PetscObject)p->X,VECDOHP,&flg);dCHK(err);
      if (flg) {err = dFSGlobalToExpanded(p->fs,p->X,p->Xexp,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);}
      else {err = VecCopy(p->X,p->Xexp);dCHK(err);}
    }
    err = VecGetArray(p->Xexp,&p->x);dCHK(err);
    err = VecGetArray(p->Yexp,&p->y);dCHK(err);
    if (!p->efs) {err = dFSGetEFS(p->fs,it->ruleset,&p->nefs,&p->efs);dCHK(err);}
    if (p->nefs != it->npatches) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Got an invalid EFS list");
    err = dRulesetIteratorLinkCreatePatchSpace_Private(p,it->maxQ);dCHK(err);
  }
  va_end(ap);
  if (!it->cjinv) {
    err = dMallocA2(it->maxQ*9,&it->cjinv,it->maxQ,&it->jw);dCHK(err);
  }
  if ((it->stash.patchbytes && !it->stash.patch) || (it->stash.nodebytes && !it->stash.node)) {
    /* Allocate space for the stash */
    err = dMallocA2(it->stash.patchbytes*it->npatches,&it->stash.patch,it->stash.nodebytes*it->nnodes,&it->stash.node);dCHK(err);
  }
  it->curpatch = 0;
  dFunctionReturn(0);
}

dErr dRulesetIteratorFinish(dRulesetIterator it)
{
  dErr err;

  dFunctionBegin;
  for (struct dRulesetIteratorLink *p=it->link; p; p=p->next) {
    err = VecRestoreArray(p->Xexp,&p->x);dCHK(err);
    err = VecRestoreArray(p->Yexp,&p->y);dCHK(err);
    if (!p->Y) continue;        /* This field is not being assembled */
    err = dFSExpandedToGlobal(p->fs,p->Yexp,p->Y,dFS_HOMOGENEOUS,ADD_VALUES);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dRulesetIteratorNextPatch(dRulesetIterator it)
{
  dErr err;

  dFunctionBegin;
  for (struct dRulesetIteratorLink *p=it->link; p; p=p->next) {
    dInt n;
    err = dEFSGetSizes(p->efs[it->curpatch],NULL,NULL,&n);dCHK(err);
    p->off += n*p->bs;
  }
  it->curpatch++;
  dFunctionReturn(0);
}

bool dRulesetIteratorHasPatch(dRulesetIterator it)
{return it->curpatch < it->npatches;}

/** dRulesetIteratorGetPatch - Get dRule, dEFSs, and fields to be evaluated on the patch
 *
 */
dErr dRulesetIteratorGetPatch(dRulesetIterator it,dRule *rule,dEFS *efs,dScalar **ex,dScalar **ey,...)
{
  dErr err;
  va_list ap;
  dInt i;
  struct dRulesetIteratorLink *p;

  dFunctionBegin;
  err = dEFSGetRule(it->link->efs[it->curpatch],rule);dCHK(err);
  va_start(ap,ey);
  for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
    if (i) {
      efs = va_arg(ap,dEFS*);
      ex = va_arg(ap,dScalar**);
      ey = va_arg(ap,dScalar**);
    }
    if (efs) *efs = p->efs[it->curpatch];
    if (ex)  *ex = &p->x[p->off];
    if (ey)  *ey = &p->y[p->off];
  }
  if (p) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"dRulesetIterator claims to have nlinks %D but linked list has more",it->nlinks);
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorGetPatchSpace - Gets space to store function values at quadrature points on the current patch
 *
 */
dErr dRulesetIteratorGetPatchSpace(dRulesetIterator it,dScalar **cjinv,dScalar **jw,dScalar **u,dScalar **du,dScalar **v,dScalar **dv,...)
{
  va_list ap;

  dFunctionBegin;
  if (cjinv) *cjinv = it->cjinv;
  if (jw)     *jw   = it->jw;
  va_start(ap,dv);
  for (struct dRulesetIteratorLink *p=it->link; p; p = p->next) {
    if (p != it->link) {
      u  = va_arg(ap,dScalar**);
      du = va_arg(ap,dScalar**);
      v  = va_arg(ap,dScalar**);
      dv = va_arg(ap,dScalar**);
    }
    if (u)  *u  = p->u;
    if (du) *du = p->du;
    if (v)  *v  = p->v;
    if (dv) *dv = p->dv;
  }
  va_end(ap);
  dFunctionReturn(0);
}

dErr dRulesetIteratorRestorePatchSpace(dRulesetIterator it,dScalar **cjinv,dScalar **jw,dScalar **u,dScalar **du,dScalar **v,dScalar **dv,...)
{
  va_list ap;

  dFunctionBegin;
  if (cjinv) *cjinv = NULL;
  if (jw)     *jw   = NULL;
  va_start(ap,dv);
  for (struct dRulesetIteratorLink *p=it->link; p; p = p->next) {
    if (p != it->link) {
      u  = va_arg(ap,dScalar**);
      du = va_arg(ap,dScalar**);
      v  = va_arg(ap,dScalar**);
      dv = va_arg(ap,dScalar**);
    }
    if (u)  *u  = NULL;
    if (du) *du = NULL;
    if (v)  *v  = NULL;
    if (dv) *dv = NULL;
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorCommitPatch - Adds contribution from the current patch to a residual vector
 *
 */
dErr dRulesetIteratorCommitPatch(dRulesetIterator it,dScalar *v,...)
{
  va_list ap;

  dFunctionBegin;
  va_start(ap,v);
  for (struct dRulesetIteratorLink *p=it->link; p; p=p->next) {
    if (p != it->link) {
      v = va_arg(ap,dScalar*);
    }
    /* Nothing to do since it goes into expanded vector, which \a v already points at */
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorGetPatchApplied - Gets a patch with function values and derivatives already evaluated on quadrature points
 *
 */
dErr dRulesetIteratorGetPatchApplied(dRulesetIterator it,dInt *Q,const dScalar **jw,dScalar **u,dScalar **du,dScalar **v,dScalar **dv,...)
{
  dErr err;
  va_list ap;
  dRule rule;
  dInt i;
  dScalar *cjinv = NULL;
  struct dRulesetIteratorLink *p;

  dFunctionBegin;
  dValidPointer(jw,3);
  err = dEFSGetRule(it->link->efs[it->curpatch],&rule);dCHK(err);
  err = dRuleGetSize(rule,0,Q);dCHK(err);
  va_start(ap,dv);
  for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
    dEFS efs;
    dScalar *ex;
    if (i) {
      u = va_arg(ap,dScalar**);
      du = va_arg(ap,dScalar**);
      v = va_arg(ap,dScalar**);
      dv = va_arg(ap,dScalar**);
    }
    efs = p->efs[it->curpatch];
    ex = &p->x[p->off];
    if (u) {
      err = dEFSApply(efs,cjinv,p->bs,ex,p->u,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
      *u = p->u;
    }
    if (du || !i) {
      err = dEFSApply(efs,cjinv,p->bs,ex,p->du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
      if (du) *du = p->du;
      if (!i) {
        cjinv = it->cjinv;
        err = dRuleComputePhysical(rule,p->du,cjinv,it->jw);dCHK(err);
        *jw = it->jw;
      }
    }
    if (v) *v = p->v;
    if (dv) *dv = p->dv;
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorCommitPatchApplied - Commits coefficients of test functions evaluated at quadrature points
 *
 */
dErr dRulesetIteratorCommitPatchApplied(dRulesetIterator it,InsertMode imode,const dScalar *v,const dScalar *dv,...)
{
  dErr err;
  va_list ap;
  dInt i;
  struct dRulesetIteratorLink *p;

  dFunctionBegin;
  va_start(ap,dv);
  for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
    dEFS efs = p->efs[it->curpatch];
    dScalar *ey = &p->y[p->off];
    if (i) {
      v = va_arg(ap,const dScalar*);
      dv = va_arg(ap,const dScalar*);
    }
    if (v) {
      err = dEFSApply(efs,it->cjinv,p->bs,v,ey,dAPPLY_INTERP_TRANSPOSE,imode);dCHK(err);
      imode = ADD_VALUES;
    }
    if (dv) {
      err = dEFSApply(efs,it->cjinv,p->bs,dv,ey,dAPPLY_GRAD_TRANSPOSE,imode);dCHK(err);
    }
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorAddStash - Attach some private memory to quadrature points and integration patches
 *
 */
dErr dRulesetIteratorAddStash(dRulesetIterator it,dInt patchbytes,dInt nodebytes)
{
  dErr err;

  dFunctionBegin;
  err = dFree(it->stash.patchoffset);dCHK(err);
  err = dFree2(it->stash.patch,it->stash.node);dCHK(err);
  it->stash.patchbytes = patchbytes;
  it->stash.nodebytes  = nodebytes;
  dFunctionReturn(0);
}

/** dRulesetIteratorGetStash - Gets a pointer to the stash for the current patch
 */
dErr dRulesetIteratorGetStash(dRulesetIterator it,void *patchstash,void *nodestash)
{
  dFunctionBegin;
  if (patchstash) *(void**)patchstash = &it->stash.patch[it->curpatch*it->stash.patchbytes];
  if (nodestash) *(void**)nodestash = &it->stash.node[it->stash.patchoffset[it->curpatch]*it->stash.nodebytes];
  dFunctionReturn(0);
}

dErr dRulesetIteratorDestroy(dRulesetIterator it)
{
  dErr err;

  dFunctionBegin;
  for (struct dRulesetIteratorLink *p=it->link,*n=p; p; p=n) {
    err = dFSRestoreEFS(p->fs,it->ruleset,&p->nefs,&p->efs);dCHK(err);
    err = dFSDestroy(p->fs);dCHK(err);
    err = VecDestroy(p->Xexp);dCHK(err);
    err = VecDestroy(p->Yexp);dCHK(err);
    err = dRulesetIteratorLinkFreePatchSpace_Private(p);dCHK(err);
    n = p->next;
    err = dFree(p);dCHK(err);
  }
  err = dFree2(it->cjinv,it->jw);dCHK(err);
  err = dFree(it);dCHK(err);
  dFunctionReturn(0);
}
