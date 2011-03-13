#include <dohpfsimpl.h>
#include <dohpvec.h>
#include <stdarg.h>

struct ValueCache {
  dInt bs,n,n_alloc;
  dScalar *u,*du,*v,*dv;
  dScalar *u_alloc,*du_alloc,*v_alloc,*dv_alloc;
};
struct ValueCachePhysical {
  dInt n_alloc;
  dReal *jw,*jw_alloc;
  dReal *cjinv,*cjinv_alloc;
};
static dErr ValueCacheSetUp(struct ValueCache *vc,dInt bs,dInt nmax);
static dErr ValueCacheReset(struct ValueCache *vc,dBool self);
static dErr ValueCacheExtract(struct ValueCache *vc,dInt n,const dInt ind[],dBool identity,dScalar su[],dScalar sdu[],dScalar sv[],dScalar sdv[]);
static dErr ValueCacheDistribute(struct ValueCache *vc,dInt n,const dInt ind[],dBool identity,dScalar su[],dScalar sdu[],dScalar sv[],dScalar sdv[]);
static dErr ValueCacheDestroy(struct ValueCache *vc);
static dErr ValueCachePhysicalSetUp(struct ValueCachePhysical *phys,dInt nmax);
static dErr ValueCachePhysicalExtract(struct ValueCachePhysical *phys,dInt n,const dInt ind[],dBool identity,dScalar *cjinv_elem,dScalar *jw_elem);
static dErr ValueCachePhysicalReset(struct ValueCachePhysical *phys);
static dErr ValueCachePhysicalDestroy(struct ValueCachePhysical *phys);

struct ExplicitCache {
  dInt maxnpatches,maxQ,maxP;
  dInt  P;                     /**< number of finite element modes on each patch */
  dInt  *bidx;                 /**< For each patch, indices of modes with support on patch. Indices are relative to expanded space. */
  dReal *interp;               /**< Interpolation: u_quad[patch*Q+i] = sum_j interp[patch*Q*P+i*P+j] * u_patch[eidx[patch*P+j]] */
  dReal *deriv;                /**< Derivatives, as above, but with gradients */
};
static dErr ExplicitCacheSetUp(struct ExplicitCache *explicit,dInt maxnpatches,dInt maxQ,dInt maxP);
static dErr ExplicitCacheDestroy(struct ExplicitCache *explicit);

struct dRulesetIteratorLink {
  dFS fs;                       /**< Function space for this link */
  const dEFS *efs;              /**< Array of element function spaces for all elements in this iterator */
  Vec X,Xexp;                   /**< Trial global and expanded vectors for residuals and Jacobian assembly */
  Vec Y,Yexp;                   /**< Test global and expanded vectors for residuals */
  dScalar *x;                   /**< Entries from trial expanded vector */
  dScalar *y;                   /**< Entries from test expanded vector */
  struct ValueCache vc_elem;    /**< Values in element-natural ordering */
  struct ValueCache vc_patch;   /**< Values in patch-natural ordering */
  dInt *rowcol;                 /**< Work array to hold row and column indices */
  dInt maxP;                    /**< Largest number of basis functions with support on this element */
  dInt nefs;                    /**< Number of dEFS */
  dInt elemstart;               /**< Offset of current element in expanded vector */
  dInt bs;                      /**< Block size */
  struct ExplicitCache explicit; /* Sparse explicit values on element */
  struct dRulesetIteratorLink *next;
};

/** User-defined storage associated with the iterator.
 *
 * Can have an arbitrary number of bytes per element and per quadrature node.
 */
struct dRulesetIteratorStash {
  char *elem;                   /**< Private data associated with the element */
  char *node;                   /**< Private data associated with each quadrature node */
  dInt *elemoffset;             /**< Node index for the first node of each element */
  dInt elembytes;               /**< Number of bytes per element (typically storing some metadata) */
  dInt nodebytes;               /**< Number of bytes per node (typically storing matrix-free Jacobian information) */
};

struct _n_dRulesetIterator {
  dRuleset ruleset;             /**< Ruleset containing quadrature rules for each patch */
  dInt curelem;                 /**< Current element number (there may be multiple patches per element) */
  dInt nelems;                  /**< Total number of elements */
  dInt curpatch_in_elem;        /**< Index of current patch in this element */
  dInt npatches_in_elem;        /**< Number of patches in current element */
  dInt patchsize;               /**< Number of nodes on each patch in current element */
  const dInt *patchind;         /**< Indices in evaluated element for each patch in current element */
  const dReal *patchweight;     /**< Patch weights for each patch in current element */
  dInt elempatch;               /**< Index of the first patch on current element with respect to expanded space */
  dInt nnodes;                  /**< Total number of nodes in expanded space */
  dInt maxnpatches;             /**< Maximum number of patches in any element in iterator */
  dInt maxQ;                    /**< Largest number of quadrature nodes in a single patch */
  dInt maxQelem;                /**< Largest number of quadrature nodes in a single element */
  dInt nlinks;                  /**< Number of dFS registered with this iterator */
  dInt Q;                       /**< Number of quadrature nodes in current patch */
  dReal *cjinv_elem;            /**< Inverse Jacobian of coordinate transformation at each quadrature node of this patch, element-natural */
  dReal *jw_elem;               /**< Physical quadrature weight (determinant of coordinate Jacobian times reference quadrature weight), element-natural */
  struct ValueCachePhysical phys; /**< Physical weights and coordinate transformation at quadrature points in patch-natural ordering */
  dScalar **Ksplit;             /**< Work space for matrix assembly between each pair of registered function spaces */
  struct dRulesetIteratorStash stash; /**< Private storage for the user */
  struct dRulesetIteratorLink *link;  /**< Links for each registered function space */
};

/** @struct dRulesetIterator
 *
 * This iterator is for efficient integration over regions and/or faces in a domain.
 *
 * Some definitions
 *
 *   An \e element is a natural unit of evaluation, such as an entire high-order element or a spline.  It is often
 *   asymptotically more efficient to evaluate a whole element than to evaluate any smaller unit.
 *
 *   A \e patch is a natural unit of integration.  It is typically the largest region with maximal order continuity and
 *   can be no larger than an element.  This is a "subelement" for composite elements and the patch between knots within
 *   a spline element.  For standard finite element methods, there is exactly one \e patch per \e element so the
 *   distinction is not normally made.
 *
 * We distinguish between \e elements and \e patches so as to maintain sparsity in the following cases:
 *
 *   \li Q_1 subelements of a nodal Q_k discretization
 *   \li splines in isogeometric analysis
 *   \li other composite elements
 *   \li reconstruction in finite volume or some discontinuous Galerkin methods
 *
 * For face patches, two elements are needed for evaluation of derivatives (provided jump terms are needed).
 *
 * This iterator hides details of the nested iteration and evaluation: the user only needs to iterate over patches.
 */

dErr dRulesetIteratorNextElement(dRulesetIterator it);
static dErr dRulesetIteratorClearElement(dRulesetIterator it);
static dErr dRulesetIteratorSetupElement(dRulesetIterator it);

#define dRulesetIteratorAssertElementSetUp(it) do { \
    if (!it->npatches_in_elem) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Element has not been set up"); \
  } while (0)

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
  err = dRulesetGetSize(rset,&it->nelems);dCHK(err);
  it->nnodes = 0;
  for (dInt i=0; i<it->nelems; i++) {
    dInt nnodes;
    err = dRuleGetSize(it->ruleset->rules[i],NULL,&nnodes);dCHK(err);
    it->nnodes += nnodes;
  }
  err = dRulesetGetMaxQ(rset,&it->maxQ,&it->maxnpatches,&it->maxQelem);dCHK(err);
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

static dErr dRulesetIteratorCreateMatrixSpace_Private(dRulesetIterator it)
{
  dErr err;
  dInt i,j;
  struct dRulesetIteratorLink *row,*col;

  dFunctionBegin;
  if (it->Ksplit) dFunctionReturn(0);
  err = dMallocA(it->nlinks*it->nlinks,&it->Ksplit);dCHK(err);
  for (i=0,row=it->link; i<it->nlinks; i++,row=row->next) {
    for (j=0,col=it->link; j<it->nlinks; j++,col=col->next) {
      err = dMallocA(row->maxP*col->maxP,&it->Ksplit[i*it->nlinks+j]);dCHK(err);
    }
  }
  dFunctionReturn(0);
}

static dErr dRulesetIteratorFreeMatrixSpace_Private(dRulesetIterator it)
{
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<dSqrInt(it->nlinks); i++) {err = dFree(it->Ksplit[i]);dCHK(err);}
  err = dFree(it->Ksplit);dCHK(err);
  dFunctionReturn(0);
}

/** Start iterating using trial vector \a X (...) and test vector \a Y (...)
 *
 * @note Collective
 */
dErr dRulesetIteratorStart(dRulesetIterator it,Vec X,Vec Y,...)
{
  dErr err;
  va_list ap;
  dInt i;
  struct dRulesetIteratorLink *p;

  dFunctionBegin;
  va_start(ap,Y);
  for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
    p->elemstart = 0;
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
    if (!p->efs) {
      err = dFSGetEFS(p->fs,it->ruleset,&p->nefs,&p->efs);dCHK(err);
      p->maxP = 0;
      for (dInt j=0; j<p->nefs; j++) {
        dInt P;
        err = dEFSGetSizes(p->efs[j],NULL,NULL,&P);dCHK(err);
        p->maxP = dMaxInt(p->maxP,P);
      }
    }
    if (p->nefs != it->nelems) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Got an invalid EFS list");
    err = ValueCacheSetUp(&p->vc_elem,p->bs,it->maxQelem);dCHK(err);
    err = ValueCacheSetUp(&p->vc_patch,p->bs,it->maxQelem);dCHK(err);
    err = ExplicitCacheSetUp(&p->explicit,it->maxnpatches,it->maxQ,p->maxP);dCHK(err);
    if (!p->rowcol) {err = dMallocA(p->maxP,&p->rowcol);dCHK(err);}
  }
  va_end(ap);
  if (!it->cjinv_elem) {
    err = dMallocA2(it->maxQelem*9,&it->cjinv_elem,it->maxQelem,&it->jw_elem);dCHK(err);
    err = ValueCachePhysicalSetUp(&it->phys,it->maxQelem);dCHK(err);
  }
  err = dRulesetIteratorCreateMatrixSpace_Private(it);dCHK(err);
  if ((it->stash.elembytes && !it->stash.elem) || (it->stash.nodebytes && !it->stash.node)) {
    /* Allocate space for the stash */
    err = dMallocA2(it->stash.elembytes*it->nelems,&it->stash.elem,it->stash.nodebytes*it->nnodes,&it->stash.node);dCHK(err);
  }
  it->curelem   = 0;
  it->elempatch = 0;
  dFunctionReturn(0);
}

/** Releases resources that were acquired for iteration and perform any collective operations to vector assembly
 *
 * @note Collective
 */
dErr dRulesetIteratorFinish(dRulesetIterator it)
{
  dErr err;

  dFunctionBegin;
  err = dRulesetIteratorClearElement(it);dCHK(err);
  for (struct dRulesetIteratorLink *p=it->link; p; p=p->next) {
    err = VecRestoreArray(p->Xexp,&p->x);dCHK(err);
    err = VecRestoreArray(p->Yexp,&p->y);dCHK(err);
    if (!p->Y) continue;        /* This field is not being assembled */
    err = dFSExpandedToGlobal(p->fs,p->Yexp,p->Y,dFS_HOMOGENEOUS,ADD_VALUES);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dRulesetIteratorNextElement(dRulesetIterator it)
{
  dErr err;

  dFunctionBegin;
  for (struct dRulesetIteratorLink *p=it->link; p; p=p->next) {
    dInt n;
    err = dEFSGetSizes(p->efs[it->curelem],NULL,NULL,&n);dCHK(err);
    p->elemstart += n*p->bs;
  }
  it->elempatch += it->npatches_in_elem;
  it->curelem++;

  err = dRulesetIteratorClearElement(it);dCHK(err);
  dFunctionReturn(0);
}

static dErr dRulesetIteratorClearElement(dRulesetIterator it)
{
  dErr err;

  dFunctionBegin;
  it->curpatch_in_elem = 0;
  it->npatches_in_elem = 0;
  it->patchsize = 0;
  it->patchind = NULL;
  it->patchweight = NULL;

  for (struct dRulesetIteratorLink *p=it->link; p; p=p->next) {
    err = ValueCacheReset(&p->vc_elem,dFALSE);dCHK(err);
    err = ValueCacheReset(&p->vc_patch,dFALSE);dCHK(err);
  }
#if defined dUSE_DEBUG
  err = dMemzero(it->cjinv_elem,3*3*it->maxQelem*sizeof(dReal));dCHK(err);
  err = dMemzero(it->jw_elem,it->maxQelem*sizeof(dReal));dCHK(err);
  err = ValueCachePhysicalReset(&it->phys);dCHK(err);
#endif
  dFunctionReturn(0);
}

static dErr dRulesetIteratorSetupElement(dRulesetIterator it)
{
  dErr  err;
  dRule rule;

  dFunctionBegin;
  err = dEFSGetRule(it->link->efs[it->curelem],&rule);dCHK(err);
  err = dRuleGetPatches(rule,&it->npatches_in_elem,&it->patchsize,&it->patchind,&it->patchweight);dCHK(err);
  err = dRuleGetSize(rule,NULL,&it->Q);dCHK(err);
  if (it->Q % it->npatches_in_elem) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Number of quadrature points %D is not divisible by number of patches %D");
  it->Q /= it->npatches_in_elem;
  dFunctionReturn(0);
}

/** Move to the next patch in iterator.
 */
dErr dRulesetIteratorNextPatch(dRulesetIterator it)
{
  dErr err;

  dFunctionBegin;
  if (it->curpatch_in_elem < it->npatches_in_elem-1) { /* We still have patches left in the current element */
    it->curpatch_in_elem++;
  } else {
    err = dRulesetIteratorNextElement(it);dCHK(err);
  }
  dFunctionReturn(0);
}

/** Check whether the iterator still has a patch.  This function cannot fail.
 *
 * @note Not collective
 */
bool dRulesetIteratorHasPatch(dRulesetIterator it)
{
  /* When patches on an element run out, curelem is advanced */
  return it->curelem < it->nelems;
}

/** dRulesetIteratorGetElement - Get dRule, dEFSs, and fields to be evaluated on the element
 *
 * @note Not collective
 * @note dRulesetIteratorStart() must have been called with non-NULL vector for each element trial function requested.
 */
dErr dRulesetIteratorGetElement(dRulesetIterator it,dRule *rule,dEFS *efs,dScalar **ex,dScalar **ey,...)
{
  dErr err;
  va_list ap;
  dInt i;
  struct dRulesetIteratorLink *p;

  dFunctionBegin;
  err = dRulesetIteratorSetupElement(it);dCHK(err);
  if (rule) {err = dEFSGetRule(it->link->efs[it->curelem],rule);dCHK(err);}

  va_start(ap,ey);
  for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
    if (i) {
      efs = va_arg(ap,dEFS*);
      ex = va_arg(ap,dScalar**);
      ey = va_arg(ap,dScalar**);
    }
    if (efs) *efs = p->efs[it->curelem];
    if (ex)  *ex = &p->x[p->elemstart];
    if (ey)  *ey = &p->y[p->elemstart];
  }
  if (p) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"dRulesetIterator claims to have nlinks %D but linked list has more",it->nlinks);
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorGetPatchSpace - Gets space to store function values at quadrature points on the current patch
 *
 * @note Not collective
 */
dErr dRulesetIteratorGetPatchSpace(dRulesetIterator it,dReal **cjinv,dReal **jw,dScalar **u,dScalar **du,dScalar **v,dScalar **dv,...)
{
  va_list ap;

  dFunctionBegin;
  if (cjinv) *cjinv = it->phys.cjinv;
  if (jw)     *jw   = it->phys.jw;
  va_start(ap,dv);
  for (struct dRulesetIteratorLink *p=it->link; p; p = p->next) {
    if (p != it->link) {
      u  = va_arg(ap,dScalar**);
      du = va_arg(ap,dScalar**);
      v  = va_arg(ap,dScalar**);
      dv = va_arg(ap,dScalar**);
    }
    if (u)  *u  = p->vc_patch.u;
    if (du) *du = p->vc_patch.du;
    if (v)  *v  = p->vc_patch.v;
    if (dv) *dv = p->vc_patch.dv;
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** Release patch space memory */
dErr dRulesetIteratorRestorePatchSpace(dRulesetIterator it,dReal **cjinv,dReal **jw,dScalar **u,dScalar **du,dScalar **v,dScalar **dv,...)
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

/** dRulesetIteratorCommitPatch - Adds contribution from coefficients of test functions on the current patch to a residual vector
 *
 * @note Not collective, but dRulesetIteratorFinish() must be called before results are guaranteed to be visible.
 * @note See dRulesetIteratorCommitPatchApplied() for a version that operators on coefficients of test functions at quadrature points
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
 * @note Not collective
 * @note dRulesetIteratorStart() must have been called with global trial vectors for any trial values that are requested here.
 */
dErr dRulesetIteratorGetPatchApplied(dRulesetIterator it,dInt *Q,const dReal **jw,dScalar **u,dScalar **du,dScalar **v,dScalar **dv,...)
{
  dErr err;
  va_list ap;
  dInt i;
  dReal *cjinv = NULL;
  struct dRulesetIteratorLink *p;

  dFunctionBegin;
  dValidPointer(jw,3);
  *jw = NULL;

  va_start(ap,dv);
  if (it->curpatch_in_elem == 0) { /* We just got to this element */
    err = dRulesetIteratorSetupElement(it);dCHK(err);
    *Q = it->patchsize;
    for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
      const dBool identity = (dBool)(it->npatches_in_elem == 1); /* Element ordering is patch ordering. Should use a better heuristic. */
      const bool has_coords = (i == 0);
      dEFS efs;
      dScalar *ex;
      if (i) {
        u = va_arg(ap,dScalar**);
        du = va_arg(ap,dScalar**);
        v = va_arg(ap,dScalar**);
        dv = va_arg(ap,dScalar**);
      }
      efs = p->efs[it->curelem];
      ex = &p->x[p->elemstart];
      err = ValueCacheReset(&p->vc_elem,dTRUE);dCHK(err);
      p->vc_elem.n = it->Q;
      if (u) {
        err = dEFSApply(efs,cjinv,p->bs,ex,p->vc_elem.u,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
      }
      if (du || has_coords) {
        err = dEFSApply(efs,cjinv,p->bs,ex,p->vc_elem.du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
        if (has_coords) {
          dRule rule;
          err = dEFSGetRule(efs,&rule);dCHK(err);
          err = dRuleComputePhysical(rule,p->vc_elem.du,it->cjinv_elem,it->jw_elem);dCHK(err);
          err = ValueCachePhysicalExtract(&it->phys,it->npatches_in_elem*it->Q,it->patchind,identity,it->cjinv_elem,it->jw_elem);dCHK(err);
          cjinv = it->cjinv_elem;      /* These coordinates are used so that the following links are mapped to the physical space */
          *jw = it->phys.jw;
        }
      }

      /* Everything is evaluated in element ordering, now get access to it in patch ordering, no-copy if identity == true */
      err = ValueCacheExtract(&p->vc_patch,it->Q*it->npatches_in_elem,it->patchind,identity,p->vc_elem.u,p->vc_elem.du,p->vc_elem.v,p->vc_elem.dv);dCHK(err);

      if (u)  *u  = p->vc_patch.u;
      if (du) *du = p->vc_patch.du;
      if (v)  *v  = p->vc_patch.v;
      if (dv) *dv = p->vc_patch.dv;
    }
  } else {                      /* Everything has already been mapped into patch-oriented storage */
    *jw = it->phys.jw + it->patchsize*it->curpatch_in_elem;
    for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
      const dInt offbs = it->patchsize*it->curpatch_in_elem*p->bs;
      if (i) {
        u = va_arg(ap,dScalar**);
        du = va_arg(ap,dScalar**);
        v = va_arg(ap,dScalar**);
        dv = va_arg(ap,dScalar**);
      }
      if (u)  *u  = p->vc_patch.u  + offbs;
      if (du) *du = p->vc_patch.du + offbs*3;
      if (v)  *v  = p->vc_patch.v  + offbs;
      if (dv) *dv = p->vc_patch.dv + offbs*3;
    }
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorCommitPatchApplied - Commits coefficients of test functions evaluated at quadrature points
 *
 * @note Not collective, but dRulesetIteratorFinish() must be called before results are guaranteed to be visible.
 */
dErr dRulesetIteratorCommitPatchApplied(dRulesetIterator it,InsertMode imode,const dScalar *v,const dScalar *dv,...)
{
  dErr err;
  va_list ap;
  dInt i;
  struct dRulesetIteratorLink *p;

  dFunctionBegin;
  if (it->curpatch_in_elem < it->npatches_in_elem-1) {
    /* Still have more patches to go in this element, no need to do anything yet. */
    dFunctionReturn(0);
  } else {
    va_start(ap,dv);
    for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
      const dBool identity = (dBool)(it->npatches_in_elem == 1); /* need better heuristic */
      dEFS efs = p->efs[it->curelem];
      dScalar *ey = &p->y[p->elemstart];
      if (i) {
        v = va_arg(ap,const dScalar*);
        dv = va_arg(ap,const dScalar*);
      }
      if (v || dv) {
        /* We don't actually reference v or dv because that memory is owned by us and sitting in vc_patch.v and
         * vc_patch.dv. Instead, we just map everything in vc_patch to vc_elem where it can be tested. */
        err = ValueCacheDistribute(&p->vc_patch,it->Q*it->npatches_in_elem,it->patchind,identity,p->vc_elem.u,p->vc_elem.du,p->vc_elem.v,p->vc_elem.dv);dCHK(err);
      }
      if (v) {
        err = dEFSApply(efs,it->cjinv_elem,p->bs,p->vc_elem.v,ey,dAPPLY_INTERP_TRANSPOSE,imode);dCHK(err);
        imode = ADD_VALUES;
      }
      if (dv) {
        err = dEFSApply(efs,it->cjinv_elem,p->bs,p->vc_elem.dv,ey,dAPPLY_GRAD_TRANSPOSE,imode);dCHK(err);
      }
    }
    va_end(ap);
  }
  dFunctionReturn(0);
}

/** Gets explicit bases for each function space on a patch
 *
 * @param[in] it iterator
 * @param[out] P number of basis functions with support on this patch (one per block)
 * @param[out] rowcol Expanded indices for each of the \a P basis functions
 * @param[out] interp Interpolation matrix for this patch, \c interp[q*P+k] is the value of basis function \c k at quadrature point \c q
 * @param[out] deriv Derivative matrix for this patch, \c deriv[(q*P+k)*3+d] is the derivative of basis function \c k at quadrature point \c q in direction \c d
 *
 * @note The vararg output tuple \c (P,rowcol,interp,deriv) is repeated once for each function space registered with dRulesetIteratorAddFS().
 *
 * @note Pass \c NULL for any of \a P, \a rowcol, \a interp, \a deriv that you are not interested in.
 */
dErr dRulesetIteratorGetPatchAssembly(dRulesetIterator it,dInt *P,const dInt **rowcol,const dReal **interp,const dReal **deriv,...)
{
  dErr err;
  va_list ap;
  dInt i;
  struct dRulesetIteratorLink *p;
  const dReal *cjinv = NULL;

  dFunctionBegin;
  dRulesetIteratorAssertElementSetUp(it);
  va_start(ap,deriv);
  if (it->curpatch_in_elem == 0) { /* We just got to this element */
    for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
      if (i) {
        P = va_arg(ap,dInt*);
        rowcol = va_arg(ap,const dInt**);
        interp = va_arg(ap,const dReal**);
        deriv = va_arg(ap,const dReal**);
      }
      if (!(P || rowcol || interp || deriv)) continue;
      err = dEFSGetExplicitSparse(p->efs[it->curelem],it->npatches_in_elem,it->Q,it->patchind,cjinv,p->elemstart,&p->explicit.P,p->explicit.bidx,p->explicit.interp,p->explicit.deriv);dCHK(err);
      if (P) *P = p->explicit.P;
      if (rowcol) *rowcol = p->explicit.bidx;
      if (interp) *interp = p->explicit.interp;
      if (deriv)  *deriv  = p->explicit.deriv;
      cjinv = it->cjinv_elem;
    }
  } else {
    for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
      if (i) {
        P = va_arg(ap,dInt*);
        rowcol = va_arg(ap,const dInt**);
        interp = va_arg(ap,const dReal**);
        deriv = va_arg(ap,const dReal**);
      }
      if (P) *P = p->explicit.P;
      if (rowcol) *rowcol = p->explicit.bidx   + p->explicit.P*it->curpatch_in_elem;
      if (interp) *interp = p->explicit.interp + it->Q*p->explicit.P*it->curpatch_in_elem;
      if (deriv)  *deriv  = p->explicit.deriv  + it->Q*p->explicit.P*it->curpatch_in_elem*3;
    }
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorGetPatchAssembly - Gets explicit bases for each function space on a patch
 *
 */
dErr dRulesetIteratorRestorePatchAssembly(dRulesetIterator it,dInt *P,const dInt **rowcol,const dReal **interp,const dReal **deriv,...)
{
  va_list ap;
  dInt i;
  struct dRulesetIteratorLink *p;
  const dReal *cjinv = NULL;

  dFunctionBegin;
  va_start(ap,deriv);
  if (it->curpatch_in_elem < it->npatches_in_elem-1) {
    for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
      if (i) {
        P = va_arg(ap,dInt*);
        rowcol = va_arg(ap,const dInt**);
        interp = va_arg(ap,const dReal**);
        deriv = va_arg(ap,const dReal**);
      }
      if (P) *P = -1;
      if (rowcol) *rowcol = NULL;
      if (interp) *interp = NULL;
      if (deriv)  *deriv  = NULL;
    }
  } else {
    for (i=0,p=it->link; i<it->nlinks; i++,p=p->next) {
      if (i) {
        P = va_arg(ap,dInt*);
        rowcol = va_arg(ap,const dInt**);
        interp = va_arg(ap,const dReal**);
        deriv = va_arg(ap,const dReal**);
      }
      if (!(P || rowcol || interp || deriv)) continue;
      /* dEFSRestoreExplicitSparse() if it existed */
      if (P) *P = -1;
      if (rowcol) *rowcol = NULL;
      if (interp) *interp = NULL;
      if (deriv)  *deriv  = NULL;
      cjinv = it->cjinv_elem;
    }
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorGetMatrixSpaceSplit - Gets storage sufficient for the constituent matrices of a coupled system
 *
 */
dErr dRulesetIteratorGetMatrixSpaceSplit(dRulesetIterator it,dScalar **K,...)
{
  va_list ap;

  dFunctionBegin;
  if (!it->Ksplit) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"No Ksplit available, must call dRulesetIteratorStart first");
  va_start(ap,K);
  for (dInt i=0; i<dSqrInt(it->nlinks); i++) {
    if (i) K = va_arg(ap,dScalar**);
    if (K) *K = it->Ksplit[i];
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetIteratorAddStash - Attach some private memory to quadrature points and integration patches
 *
 */
dErr dRulesetIteratorAddStash(dRulesetIterator it,dInt elembytes,dInt nodebytes)
{
  dErr err;

  dFunctionBegin;
  err = dFree(it->stash.elemoffset);dCHK(err);
  err = dFree2(it->stash.elem,it->stash.node);dCHK(err);
  it->stash.elembytes = elembytes;
  it->stash.nodebytes  = nodebytes;
  dFunctionReturn(0);
}

/** dRulesetIteratorGetStash - Gets a pointer to the stash for the current element
 */
dErr dRulesetIteratorGetStash(dRulesetIterator it,void *elemstash,void *nodestash)
{
  dFunctionBegin;
  if (elemstash) *(void**)elemstash = &it->stash.elem[it->curelem*it->stash.elembytes];
  if (nodestash) *(void**)nodestash = &it->stash.node[it->stash.elemoffset[it->curelem]*it->stash.nodebytes];
  dFunctionReturn(0);
}

/** Destroy an iterator when it is no longer needed.
 *
 * @note This function should only be called after dRulesetIteratorFinish().
 * @note Iterators can be reused by calling dRulesetIteratorStart() again.
 */
dErr dRulesetIteratorDestroy(dRulesetIterator it)
{
  dErr err;

  dFunctionBegin;
  for (struct dRulesetIteratorLink *p=it->link,*n=p; p; p=n) {
    err = dFSRestoreEFS(p->fs,it->ruleset,&p->nefs,&p->efs);dCHK(err);
    err = dFSDestroy(p->fs);dCHK(err);
    err = VecDestroy(p->Xexp);dCHK(err);
    err = VecDestroy(p->Yexp);dCHK(err);
    err = ValueCacheDestroy(&p->vc_elem);dCHK(err);
    err = ValueCacheDestroy(&p->vc_patch);dCHK(err);
    err = ExplicitCacheDestroy(&p->explicit);dCHK(err);
    err = dFree(p->rowcol);dCHK(err);
    n = p->next;
    err = dFree(p);dCHK(err);
  }
  err = dRulesetIteratorFreeMatrixSpace_Private(it);dCHK(err);
  err = dFree2(it->cjinv_elem,it->jw_elem);dCHK(err);
  err = ValueCachePhysicalDestroy(&it->phys);dCHK(err);
  err = dFree(it);dCHK(err);
  dFunctionReturn(0);
}


static dErr ValueCacheSetUp(struct ValueCache *vc,dInt bs,dInt nmax)
{
  dErr err;

  dFunctionBegin;
  if (bs != vc->bs || nmax > vc->n_alloc) {
    err = dFree4(vc->u_alloc,vc->du_alloc,vc->v_alloc,vc->dv_alloc);dCHK(err);
    err = dMemzero(vc,sizeof(*vc));dCHK(err);
    vc->bs = bs;
    vc->n_alloc = nmax;
    err = dMallocA4(bs*nmax,&vc->u_alloc,bs*nmax*3,&vc->du_alloc,bs*nmax,&vc->v_alloc,bs*nmax*3,&vc->dv_alloc);dCHK(err);
  }
  err = ValueCacheReset(vc,dFALSE);dCHK(err);
  dFunctionReturn(0);
}

static dErr ValueCacheReset(struct ValueCache *vc,dBool self)
{
  const dInt nbs = vc->n_alloc * vc->bs;

  dFunctionBegin;
  vc->n = 0;
  if (self) {
    vc->u  = vc->u_alloc;
    vc->du = vc->du_alloc;
    vc->v  = vc->v_alloc;
    vc->dv = vc->dv_alloc;
  } else {
    vc->u  = NULL;
    vc->du = NULL;
    vc->v  = NULL;
    vc->dv = NULL;
  }
  dMakeMemUndefined(vc->u_alloc,nbs*sizeof(dScalar));
  dMakeMemUndefined(vc->du_alloc,3*nbs*sizeof(dScalar));
  dMakeMemUndefined(vc->v_alloc,nbs*sizeof(dScalar));
  dMakeMemUndefined(vc->dv_alloc,3*nbs*sizeof(dScalar));
  dFunctionReturn(0);
}

/** Extract values from a scattered array into a contiguous one, shares memory if identity == true */
static dErr ValueCacheExtract(struct ValueCache *vc,dInt n,const dInt *ind,dBool identity,dScalar su[],dScalar sdu[],dScalar sv[],dScalar sdv[])
{
  dErr err;

  dFunctionBegin;
  if (n > vc->n_alloc) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Attempt to extract %D values, but cache size is %D",n,vc->n_alloc);
  if (identity) {
    vc->n = n;
    vc->u  = su;
    vc->du = sdu;
    vc->v  = sv;
    vc->dv = sdv;
  } else {
    const dInt bs = vc->bs;
    err = ValueCacheReset(vc,dTRUE);dCHK(err);
    vc->n = n;
    for (dInt i=0; i<n; i++) {
      for (dInt j=0; j<bs; j++) {
        vc->u[i*bs+j]        = su[ind[i]*bs+j];
        vc->du[(i*bs+j)*3+0] = sdu[(ind[i]*bs+j)*3+0];
        vc->du[(i*bs+j)*3+1] = sdu[(ind[i]*bs+j)*3+1];
        vc->du[(i*bs+j)*3+2] = sdu[(ind[i]*bs+j)*3+2];
      }
    }
  }
  dFunctionReturn(0);
}

static dErr dUNUSED ValueCacheDistribute(struct ValueCache *vc,dInt n,const dInt *ind,dBool identity,dScalar su[],dScalar sdu[],dScalar sv[],dScalar sdv[])
{

  dFunctionBegin;
  if (n > vc->n_alloc) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Attempt to extract %D values, but cache size is %D",n,vc->n_alloc);
  if (n != vc->n) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Distribute size %D different from extracted size %D",n,vc->n);
  if (!vc->v || !vc->dv) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Cache is empty");
  if (identity) {
    if (!(sv == vc->v && sdv == vc->dv)) { /* This case is recoverable, but probably a logic error */
      dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,
             ((vc->v == vc->v_alloc)
              ? "Identity was not used when the cache was filled but new vector came from a different source" /* Recoverable but probably a logic error */
              : "Cache is aliasing a different vector"));
    }
    if (!(su == vc->u && sdu == vc->du)) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Outgoing vectors match, but incoming vectors are not aliased as expected, suspect logic error");
    /* Do nothing because vectors are aliased so values already went into dest */
  } else {
    const dInt bs = vc->bs;
    for (dInt i=0; i<n; i++) {
      for (dInt j=0; j<bs; j++) {
        const dInt ibs = i*bs+j,sibs = ind[i]*bs+j;
        sv[sibs] += vc->v[ibs];
        sdv[sibs*3+0] += vc->dv[ibs*3+0];
        sdv[sibs*3+1] += vc->dv[ibs*3+1];
        sdv[sibs*3+2] += vc->dv[ibs*3+2];
      }
    }
  }
  dFunctionReturn(0);
}

static dErr ValueCacheDestroy(struct ValueCache *vc)
{
  dErr err;

  dFunctionBegin;
  err = dFree4(vc->u_alloc,vc->du_alloc,vc->v_alloc,vc->dv_alloc);dCHK(err);
  err = dMemzero(vc,sizeof(*vc));dCHK(err);
  dFunctionReturn(0);
}

static dErr ValueCachePhysicalSetUp(struct ValueCachePhysical *phys,dInt nmax)
{
  dErr err;

  dFunctionBegin;
  if (phys->n_alloc) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Cache already set up");
  phys->n_alloc = nmax;
  err = dMallocA2(nmax*9,&phys->cjinv_alloc,nmax,&phys->jw_alloc);dCHK(err);
  for (dInt i=0; i<nmax; i++) phys->jw_alloc[i] = -1;
  phys->cjinv = NULL;
  phys->jw = NULL;
  dFunctionReturn(0);
}

static dErr ValueCachePhysicalExtract(struct ValueCachePhysical *phys,dInt n,const dInt ind[],dBool identity,dScalar cjinv_elem[],dScalar jw_elem[])
{
  dFunctionBegin;
  if (n > phys->n_alloc) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Not enough space allocated in the physical map cache");
  if (identity) {
    phys->cjinv = cjinv_elem;
    phys->jw = jw_elem;
  } else {
    phys->cjinv = phys->cjinv_alloc;
    phys->jw = phys->jw_alloc;
    for (dInt i=0; i<n; i++) {
      for (dInt j=0; j<9; j++) phys->cjinv[i*9+j] = cjinv_elem[ind[i]*9+j];
      phys->jw[i] = jw_elem[ind[i]];
    }
  }
  dFunctionReturn(0);
}

static dErr ValueCachePhysicalReset(struct ValueCachePhysical *phys)
{
  dErr err;
  dInt n = phys->n_alloc;

  dFunctionBegin;
  err = dMemzero(phys->cjinv_alloc,9*n*sizeof(dReal));dCHK(err);
  err = dMemzero(phys->jw_alloc,n*sizeof(dReal));dCHK(err);
  phys->cjinv = NULL;
  phys->jw = NULL;
  dFunctionReturn(0);
}

static dErr ValueCachePhysicalDestroy(struct ValueCachePhysical *phys)
{
  dErr err;

  dFunctionBegin;
  err = dFree2(phys->cjinv_alloc,phys->jw_alloc);dCHK(err);
  err = dMemzero(phys,sizeof(*phys));dCHK(err);
  dFunctionReturn(0);
}

static dErr ExplicitCacheSetUp(struct ExplicitCache *explicit,dInt maxnpatches,dInt maxQ,dInt maxP)
{
  dErr err;

  dFunctionBegin;
  if (!(maxnpatches <= explicit->maxnpatches && maxQ <= explicit->maxQ && maxP <= explicit->maxP)) {
    err = dFree3(explicit->bidx,explicit->interp,explicit->deriv);dCHK(err);
    explicit->maxnpatches = maxnpatches;
    explicit->maxQ        = maxQ;
    explicit->maxP        = maxP;
    explicit->P = -1;
    err = dMallocA3(maxnpatches*maxP,&explicit->bidx,maxnpatches*maxQ*maxP,&explicit->interp,maxnpatches*maxQ*maxP*3,&explicit->deriv);dCHK(err);
  }
  err = dMemzero(explicit->bidx,maxnpatches*maxP*sizeof(explicit->bidx[0]));dCHK(err);
  err = dMemzero(explicit->interp,maxnpatches*maxQ*maxP*sizeof(explicit->interp[0]));dCHK(err);
  err = dMemzero(explicit->deriv,maxnpatches*maxQ*maxP*3*sizeof(explicit->deriv[0]));dCHK(err);
  dFunctionReturn(0);
}

static dErr ExplicitCacheDestroy(struct ExplicitCache *explicit)
{
  dErr err;

  dFunctionBegin;
  err = dFree3(explicit->bidx,explicit->interp,explicit->deriv);dCHK(err);
  dFunctionReturn(0);
}
