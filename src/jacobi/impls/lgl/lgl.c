/**
* @file   lgl.c
* @author Jed Brown <jed@59A2.org>
* @date   Fri Aug 22 20:39:12 2008
* 
* @brief  A nodal Legendre-Gauss-Lobatto basis and Gauss quadrature.
* 
* 
*/

#include "petsc.h"
#include "dohpjacobi.h"
#include "uthash.h"
// #include "src/util/polylib.c"

/* If pointer is not aligned to a cache line boundary, advance it so that it is. */

#define CACHE_LINE    64l       /* my cache lines are 64 bytes long */
#define DEFAULT_ALIGN 16l       /* SSE instructions require 16 byte alignment */

#define dNextCacheAligned(p) dNextAlignedAddr(CACHE_LINE,p)
#define dNextAligned(p)      dNextAlignedAddr(DEFAULT_ALIGN,p)

/** 
* Returns the next address which satisfies the given alignment.
*
* This function cannot fail.
* 
* @param alignment must be a power of 2
* @param ptr The pointer
* 
* @return aligned address
*/
static inline void *dNextAlignedAddr(size_t alignment,void *ptr)
{
  size_t base = (size_t)ptr;
  size_t mask = alignment-1;
  if (base & mask) return (void*)(base + (alignment - (base & mask)));
  return (void*)base;
}

struct HashBasis {
  dInt id;
  struct m_dBasis *basis;
  UT_hash_handle hh;
};

struct LGL {
  dErr (*hasbasis)(dJacobi,dInt,dInt,dBool*);
  struct m_dRule *rule;
  struct m_dBasis *basis;
  dReal *ar,*ab;
  dInt arrLen;
  dInt *bInd,*bRowOff;
  dReal *work;
};

static dErr dJacobiSetUp_LGL(dJacobi);
static dErr dJacobiDestroy_LGL(dJacobi);
static dErr dJacobiView_LGL(dJacobi,PetscViewer);
static dErr dJacobiGetRule_LGL(dJacobi,dInt,dRule*);
static dErr dJacobiGetBasis_LGL(dJacobi,dInt,dInt,dBasis*);

static dErr LGLNodeWeight(struct LGL*,dInt,dReal*,dReal*);
static dErr LGLBasisDeriv(struct LGL*,dInt,dInt,dReal*,dReal*);
static dErr LGLHasBasis(dJacobi,dInt,dInt,dBool*);



/** 
* Initializes the ops table.
* 
* @param jac 
* 
* @return 
*/
dErr dJacobiCreate_LGL(dJacobi jac)
{
  static const struct _dJacobiOps myops = {
    .setup = dJacobiSetUp_LGL,
    .setfromoptions = 0,
    .destroy = dJacobiDestroy_LGL,
    .view = dJacobiView_LGL,
    .getrule = dJacobiGetRule_LGL,
    .getbasis = dJacobiGetBasis_LGL
  };
  dErr err;

  dFunctionBegin;
  err = dPrintf(((PetscObject)jac)->comm,"dJacobiCreate_LGL (need to do something here)\n");dCHK(err);
  err = PetscMemcpy(jac->ops,&myops,sizeof(struct _dJacobiOps));dCHK(err);
  err = PetscNew(struct LGL,&jac->impl);dCHK(err);
  ((struct LGL*)jac->impl)->hasbasis = LGLHasBasis;
  dFunctionReturn(0);
}



static dErr dJacobiSetUp_LGL(dJacobi jac)
{
  struct LGL *impl;
  struct m_dRule *rule;
  dReal *aspace,*ar;
  dInt bsize,rsize,bframe,rspace;
  dErr err;

  dFunctionBegin;
  bsize = jac->basisdegree;
  rsize = bsize + jac->ruleexcess;
  impl = (struct LGL*)(jac->impl);
  /* Get space to store all the rule contexts as an array, also generic workspace */
  err = PetscMalloc2(rsize,struct m_dRule,&impl->rule,rsize*bsize*4,dReal,&impl->work);dCHK(err);
  err = PetscMemzero(impl->rule,rsize*sizeof(struct m_dRule));dCHK(err);
  /* Get space for the rule data (yes this is just a simple formula) */
  rspace = 0;
  for (dInt i=1; i<rsize; i++) {
    rspace += 2*i;
  }
  err = PetscMalloc(rspace*sizeof(dReal),&impl->ar);dCHK(err);
  rule = impl->rule;
  ar = impl->ar;

  bframe = 0;                   /* number of active basis/rule pairs */
  aspace = 0;                   /* amount of space that the data will take */
  for (dInt i=1, rz=0; i<rsize; i++, rz++) {
    err = LGLNodeWeight(impl,i,&ar[rz],&ar[rz+i]);dCHK(err);
    rule[i].node = &ar[rz];
    rule[i].weight = &ar[rz+i];
    rz += 2*i;
    for (dInt j=0; j<bsize; j++) {   /* Compute sparsity for basis setup */
      dBool has;
      err = impl->hasbasis(jac,i,j,&has);dCHK(err);
      if (has) {
        bframe++;
        aspace = dNextCacheAligned(aspace + i*j); /* space for basis matrix */
        aspace = dNextCacheAligned(aspace + i*j); /* space for derivative matrix */
      }
    }
  }

  /**
  * We store the basis table in CSR form (rule,basis), therefore we need:
  *   - space for basis frames \c bframe
  *   - space for column (basis) indices
  *   - space for offsets into column indices (bsize+1)
  *   - space for actual values pointed to by the basis frames \c aspace
  *
  * Note: We want the data arrays (basis and derivative matrices) to be aligned to cache line boundaries, hence the
  *   awkward declaration.
  */
  err = PetscMalloc4(bframe,struct m_dBasis,&impl->basis,bframe,dInt,&impl->bInd,rsize+1,dInt,&impl->bRowOff,
                     (size_t)(aspace-(dReal*)0)+CACHE_LINE/sizeof(dReal),dReal,&impl->ab);dCHK(err);

  {
    dReal *az = dNextCacheAligned(impl->ab); /* The current finger in the matrix buffer */
    dInt bi = 0;                           /* Index in column array (\c bInd) */
    for (dInt i=0; i<rsize; i++) {
      /* Create the rule */
      impl->bRowOff[i] = bi;
      for (dInt j=0; j<bsize; j++, bi++) {
        dReal *basis,*deriv;
        dBool has;
        err = impl->hasbasis(jac,i,j,&has);dCHK(err);
        if (!has) break;
        impl->bInd[bi] = j;
        basis = az; az = dNextCacheAligned(az+i*j);
        deriv = az; az = dNextCacheAligned(az+i*j);
        err = LGLBasisDeriv(impl,i,j,basis,deriv);dCHK(err);
        impl->basis[bi].dmat = basis;
        impl->basis[bi].dmat = deriv;
      }
    }
  }
  dFunctionReturn(0);
}

static dErr dJacobiDestroy_LGL(dJacobi jac)
{
  struct LGL *impl = (struct LGL*)jac->impl;
  dErr err;

  dFunctionBegin;
  err = PetscFree2(impl->rule,impl->work);dCHK(err);
  err = PetscFree(impl->ar);dCHK(err);
  dFunctionReturn(0);
}



static dErr LGLHasBasis(dJacobi jac,dInt rule,dInt basis,dBool *has)
{

  dFunctionBegin;
  *has = (1 < basis && basis < jac->basisdegree && basis <= rule && rule < jac->basisdegree+jac->ruleexcess);
  dFunctionReturn(0);
}


static dErr LGLNodeWeight(struct LGL *lgl,dInt n,dReal node[],dReal weight[])
{
  dReal *restrict a, *restrict b;
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<n; i++) {
    node[i] = -1.0 + 0.5*i/(n-1.0);
    weight[i] = 1.0 / n;
  }
  dFunctionReturn(0);
}

