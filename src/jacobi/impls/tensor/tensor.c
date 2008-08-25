/**
* @file   tensor.c
* @author Jed Brown <jed@59A2.org>
* @date   Fri Aug 22 20:39:12 2008
* 
* @brief  A nodal Tensor product basis.
*
* See the Tensor section for details on the storage method used.
* 
*/

#include "petsc.h"
#include "dohpjacobi.h"
// #include "private/dohpimpl.h"
#include "private/fsimpl.h"
// #include "uthash.h"
#include "inlinepoly.h"
#include "blist.h"

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

typedef enum { GAUSS, GAUSS_LOBATTO, GAUSS_RADAU } GaussFamily;

/**
* The Rule and Basis building functions need work space that is freeable later (mostly just so that valgrind won't
* complain, it's not a serious memory leak, but this design is more flexible than internal static memory management.
* 
*/
typedef struct {
  void *options;
  dInt workLength;
  dReal *work;
} TensorBuilder;
static dErr TensorBuilderInit(void*,TensorBuilder*);
static dErr TensorBuilderDestroy(TensorBuilder*);

typedef struct {
  dReal       alpha,beta;
  GaussFamily family;
} TensorRuleOptions;

/**
* We make no particular attempt to align the beginning of the structure, however we would like to align the arrays,
* especially the large derivative arrays, at least to 16-byte boundaries (to enable SSE) and preferably to cache line
* boundaries.
* 
*/
typedef struct {
  dInt  size;                     /**< number of quadrature points */
  dInt  weightOffset,coordOffset; /**< index into \c data where coordinates start  */
  dReal data[];                   /**< weights = data[0...size], nodes = data[coordOffset...coordOffset+size] */
} TensorRule;
static dErr TensorRuleCreate(TensorBuilder*,dInt,dInt,TensorRule*,dInt*);

typedef struct {
  dReal       alpha,beta;
  GaussFamily family;
} TensorBasisOptions;
typedef struct {
  dErr (*mult)(dInt,dInt,dInt,dInt,dInt,dInt,dInt,const dReal[],const dReal[],const dReal[],dReal[]);
  dErr (*multtrans)(dInt,dInt,dInt,dInt,dInt,dInt,dInt,const dReal[],const dReal[],const dReal[],dReal[]);
  dInt  ruleSize,basisSize;
  dInt  basisOffset,derivOffset,nodeOffset;
  dReal data[];
} TensorBasis;
static dErr TensorBasisCreate(TensorBuilder*,const TensorRule*,dInt,dInt,TensorBasis*,dInt*);

typedef struct m_Tensor *Tensor;
/**
* There are several factors in play which justify the storage method used.  First, we would like rapid lookup of
* TensorRule and TensorBasis objects since many lookups are needed any time the approximation order changes.  More
* importantly, we would like to be able to generate the best possible code and data alignment for the kernel operations,
* particularly multiplication of derivative and interpolation matrices across the data arrays.  For this, it is ideal to
* have these matrices aligned to a cache line boundary.  This leaves the maximum possible cache available for data and
* may allow SSE instructions, especially if we can find a way to guarantee that the data arrays are also at least
* 16-byte aligned.
*
* To achieve this, we will use raw buffers to store the TensorRule and TensorBasis.  The dBufferList is a simple way to
* manage these buffers since it is difficult to determine in advance how much room is needed.  When filling in a new
* TensorRule or TensorBasis fails, we just get a new base pointer from dBufferList and retry.
*
* For rule data, it is cheap to just generate rules for every size up to the maximum needed.  For basis data, it isn't
* necessarily so simple since there could be lots of Basis/Rule combinations when the spectral order becomes very large.
* To cope with this, we store our collection of TensorBasis is CSR format.  That is, the array \c basisOffsetByRule
* holds the start of the bases for each rule size.  To find a basis on a rule \a m with \a n functions, we search for
* order \a n in 'basisSizeFlat[]' between indices range 'basisOffsetByRule[m]' to 'basisOffsetByRule[m+1]-1'.
* 
*/
struct m_Tensor {
  TensorBuilder ruleBuilder,basisBuilder;
  dBufferList ruleData,basisData; /**< The actual storage for rules. */
  TensorRule  **rule;              /**< Array of 1D rule pointers.  The rule with \a m points is indexed as \a rule[\a m]. */
  TensorBasis **basis;             /**< Array of 1D basis pointers, length \a N.  Consider basis with \a m quadrature
                                  * points and \a n basis functions.
                                  *
                                  * If the quadrature points are \f$ q_i, i = 0,1,\dotsc,m-1 \f$
                                  *
                                  * and the nodes are \f$ x_j, j=0,1,\dotsc,n-1 \f$
                                  * 
                                  * then \f$ f(q_i) = \sum_{j=0}^n B[i*n+j] f(x_j) \f$ */
  dInt M,N;                     /**< length of arrays \a rule and \a basis */
  dInt *basisOffsetByRule;      /**< Array of length \a M+1, that is, the number of rules plus one. */
  dInt *basisSizeFlat;          /**< Array of length \a N,  */
  TensorRuleOptions ruleOpts;
  TensorBasisOptions basisOpts;
  dBufferList data;
  struct v_dRuleOps *ruleOpsLine,*ruleOpsQuad;
  struct v_dEFSOps *efsOpsLine,*efsOpsQuad;
  dBool setupcalled;
};

static dErr dJacobiSetUp_Tensor(dJacobi);
static dErr dJacobiDestroy_Tensor(dJacobi);
static dErr dJacobiView_Tensor(dJacobi,PetscViewer);
static dErr dJacobiGetRule_Tensor(dJacobi jac,dTopology top,const dInt rsize[],dInt left,dRule *rule,dInt *bytes);
static dErr dJacobiGetEFS_Tensor(dJacobi jac,dTopology top,const dInt bsize[],const dRule *rule,dInt left,dEFS *efs,dInt *bytes);

static dErr TensorGetRule(Tensor this,dInt n,TensorRule **out);
static dErr TensorGetBasis(Tensor this,dInt m,dInt n,TensorBasis **out);

static dErr TensorJacobiHasBasis(dJacobi,dInt,dInt,dBool*);

static dErr dRealTableView(dInt m,dInt n,const dReal mat[],const char name[],dViewer viewer);


/** 
* Initializes the ops table.
* 
* @param jac 
* 
* @return 
*/
dErr dJacobiCreate_Tensor(dJacobi jac)
{
  static const struct _dJacobiOps myops = {
    .setup = dJacobiSetUp_Tensor,
    .setfromoptions = 0,
    .destroy = dJacobiDestroy_Tensor,
    .view = dJacobiView_Tensor,
    .getrule = dJacobiGetRule_Tensor,
    .getefs = dJacobiGetEFS_Tensor
  };
  TensorRuleOptions *ropt;
  TensorBasisOptions *bopt;
  dErr err;

  dFunctionBegin;
  err = dPrintf(((PetscObject)jac)->comm,"dJacobiCreate_Tensor()\n");dCHK(err);
  err = dMemcpy(jac->ops,&myops,sizeof(struct _dJacobiOps));dCHK(err);
  err = dNew(struct m_Tensor,&jac->impl);dCHK(err);

  ropt         = &((Tensor)jac->impl)->ruleOpts;
  ropt->alpha  = 0.0;
  ropt->beta   = 0.0;
  ropt->family = GAUSS;

  bopt         = &((Tensor)jac->impl)->basisOpts;
  bopt->alpha  = 0.0;
  bopt->beta   = 0.0;
  bopt->family = GAUSS_LOBATTO;
  dFunctionReturn(0);
}

/** 
* Prepare the dJacobi context to return dRule and dEFS objects (with dJacobiGetRule and dJacobiGetEFS).
* 
* @param jac the context
* 
* @return err
*/
static dErr dJacobiSetUp_Tensor(dJacobi jac)
{
  Tensor impl = (Tensor)(jac->impl);
  dInt   M,N;
  dErr   err;

  dFunctionBegin;
  err = TensorBuilderInit((void*)&impl->ruleOpts,&impl->ruleBuilder);dCHK(err);
  err = TensorBuilderInit((void*)&impl->basisOpts,&impl->basisBuilder);dCHK(err);
  impl->N = N = jac->basisdegree;                     /* all valid basis degrees are < P */
  impl->M = M = N + jac->ruleexcess;                  /* all valid rule degrees are < M */

  err = dMallocM(M,TensorRule*,&impl->rule);dCHK(err); /* Get space to store all the rule pointers as an array. */
  impl->rule[0] = NULL;                               /* There is no rule with zero points. */
  {
    const size_t chunkSize = sizeof(dReal)*M*M/2 + M*sizeof(impl->rule[0]); /* Not quite enough because of alignment issues */
    TensorRule *z,*znext;
    dInt        used,left = 0;
    for (dInt i=1; i<M; i++) {
      do {
        used = 0;
        if (!left) {
          err = dBufferListMalloc(&impl->ruleData,chunkSize,(void**)&z);dCHK(err);
          left = (dInt)chunkSize;
        }
        err = TensorRuleCreate(&impl->ruleBuilder,i,left,z,&used);dCHK(err); /* If there is insufficient space 'left', used=0 */
        znext = dNextAligned((void*)((size_t)z + used));
        left -= (dInt)((size_t)znext - (size_t)z);
        if (!used) left = 0;
      } while (!used);          /* This should loop once at most. */
      impl->rule[i] = z;
      z = znext;
    }
  }

  err = dMallocM(M*N,TensorBasis*,&impl->basis);dCHK(err);
  {
    const size_t chunkSize = M*N*((2*M*N+N)*sizeof(dReal) + sizeof(TensorBasis))/2; /* Adjust this for alignment. */
    TensorBasis *z,*znext;
    dInt         used,left = 0;
    dBool        has;

    for (dInt i=0; i<M; i++) {
      TensorRule *rule = impl->rule[i];
      for (dInt j=0; j<N; j++) {
        err = TensorJacobiHasBasis(jac,i,j,&has);dCHK(err);
        if (!has || !rule) {
          impl->basis[i*N+j] = NULL;
          continue;
        }
        do {
          used = 0;
          if (!left) {
            err = dBufferListMalloc(&impl->basisData,chunkSize,(void**)&z);dCHK(err);
            left = (dInt)chunkSize;
          }
          err = TensorBasisCreate(&impl->basisBuilder,rule,j,left,z,&used);dCHK(err);
          znext = dNextAligned((void*)((size_t)z + used));
          left -= (dInt)((size_t)znext - (size_t)z);
          if (!used) left = 0;
        } while (!used);        /* should loop at most once */
        impl->basis[i*N+j] = z;
        z = znext;
      }
    }
  }
  impl->setupcalled = true;
  dFunctionReturn(0);
}

static dErr dJacobiDestroy_Tensor(dJacobi jac)
{
  Tensor impl = (Tensor)jac->impl;
  dErr err;

  dFunctionBegin;
  err = dFree(impl->rule);dCHK(err);
  err = dFree(impl->basis);dCHK(err);
  err = dBufferListFree(&impl->ruleData);dCHK(err);
  err = dBufferListFree(&impl->basisData);dCHK(err);
  err = TensorBuilderDestroy(&impl->ruleBuilder);dCHK(err);
  err = TensorBuilderDestroy(&impl->basisBuilder);dCHK(err);
  dFunctionReturn(0);
}

static dErr dJacobiView_Tensor(dJacobi jac,dViewer viewer)
{
  Tensor this = (Tensor)jac->impl;
  dBool ascii;
  const TensorRule *r;
  TensorBasis *b;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  err = PetscViewerASCIIPrintf(viewer,"Tensor based Jacobi\n");dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"TensorRule database.\n");dCHK(err);
  for (dInt i=0; i<this->M; i++) {
    r = this->rule[i];
    if (r) {
      err = PetscViewerASCIIPrintf(viewer,"TensorRule with %d nodes.\n",r->size);dCHK(err);
      err = dRealTableView(1,i,&r->data[r->coordOffset],"q",viewer);dCHK(err);
      err = dRealTableView(1,i,&r->data[r->weightOffset],"w",viewer);dCHK(err);
    }
  }
  err = PetscViewerASCIIPrintf(viewer,"TensorBasis database.\n");dCHK(err);
  for (dInt i=1; i<this->M; i++) {
    for (dInt j=1; j<this->N; j++) {
      b = 0;
      err = TensorGetBasis(this,i,j,&b);dCHK(err);
      if (b) {
        err = PetscViewerASCIIPrintf(viewer,"TensorBasis with nodes=%d basis=%d.\n");dCHK(err);
        err = dRealTableView(i,j,&b->data[b->basisOffset],"basis",viewer);dCHK(err);
      }
    }
  }
  /* view the basis functions next */
  dFunctionReturn(0);
}

static dErr dRealTableView(dInt m,dInt n,const dReal mat[],const char *name,dViewer viewer)
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  for (dInt i=0; i<m; i++) {
    if (name) {
      err = PetscViewerASCIIPrintf(viewer,"%s[%d][%d:%d] ",name,i,0,n-1);dCHK(err);
    }
    for (dInt j=0; j<n; j++) {
      err = PetscViewerASCIIPrintf(viewer,"%12.5g ",mat[i*n+j]);dCHK(err);
    }
    err = PetscViewerASCIIPrintf(viewer,"\n");dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr TensorJacobiHasBasis(dJacobi jac,dInt rule,dInt basis,dBool *has)
{
  Tensor this = jac->impl;

  dFunctionBegin;
  *has = (1 < basis && basis < this->N && basis <= rule && rule < this->M);
  dFunctionReturn(0);
}

/** 
* Get the dRule corresponding to the given sizes.
* 
* @param jac Jacobi context
* @param top topology
* @param rsize rule size in each dimension, normally
* @param left number of bytes left in the buffer which \a rule points to
* @param[in,out] rule space to write the rule
* @param bytes number of bytes taken by rule
* 
* @return 
*/
static dErr dJacobiGetRule_Tensor(dJacobi jac,dTopology top,const dInt rsize[],dInt left,dRule *rule,dInt *bytes)
{
  Tensor this = (Tensor)jac->impl;
  dInt needed;
  dErr err;

  dFunctionBegin;
  switch (top) {
    case iMesh_LINE_SEGMENT:
      needed = offsetof(struct m_dRule,data[1]);
      if (needed > left) dERROR(1,"not enough space in buffer");
      *bytes = needed;
      err = TensorGetRule(this,rsize[0],(TensorRule**)&rule->data[0]);dCHK(err);
      rule->ops = this->ruleOpsLine;
      break;
    case iMesh_QUADRILATERAL:
      needed = offsetof(struct m_dRule,data[2]);dCHK(err);
      if (needed > left) dERROR(1,"not enough space in buffer");
      *bytes = needed;
      err = TensorGetRule(this,rsize[0],(TensorRule**)&rule->data[0]);dCHK(err);
      err = TensorGetRule(this,rsize[1],(TensorRule**)&rule->data[1]);dCHK(err);
      rule->ops = this->ruleOpsQuad;
      break;
    case iMesh_HEXAHEDRON:
    default:
      dERROR(1,"no rule available for given topology");
  }
  dFunctionReturn(0);
}

static dErr dJacobiGetEFS_Tensor(dJacobi jac,dTopology top,const dInt bsize[],const dRule *rule,dInt left,dEFS *efs,dInt *bytes)
{
  Tensor this = (Tensor)jac->impl;
  dInt needed;
  dErr err;

  dFunctionBegin;
  switch (top) {
    case iMesh_LINE_SEGMENT:
      needed = offsetof(struct m_dEFS,data[1]);
      if (needed > left) dERROR(1,"not enough space in buffer");
      *bytes = needed;
      err = TensorGetBasis(this,((TensorRule*)rule->data[0])->size,bsize[0],(TensorBasis**)&efs->data[0]);dCHK(err);
      efs->ops = this->efsOpsLine;
      break;
    case iMesh_QUADRILATERAL:
      needed = offsetof(struct m_dEFS,data[2]);dCHK(err);
      if (needed > left) dERROR(1,"not enough space in buffer");
      *bytes = needed;
      err = TensorGetBasis(this,((TensorRule*)rule->data[0])->size,bsize[0],(TensorBasis**)&efs->data[0]);dCHK(err);
      err = TensorGetBasis(this,((TensorRule*)rule->data[1])->size,bsize[1],(TensorBasis**)&efs->data[1]);dCHK(err);
      efs->ops = this->efsOpsQuad;
      break;
    case iMesh_HEXAHEDRON:
    default:
      dERROR(1,"no rule available for given topology");
  }
  dFunctionReturn(0);
}
  

static dErr TensorBuilderInit(void *opts,TensorBuilder *out)
{
  dFunctionBegin;
  out->options = opts;
  out->workLength = 0;
  out->work = NULL;
  dFunctionReturn(0);
}

static dErr TensorBuilderDestroy(TensorBuilder *build)
{
  dErr err;
  
  dFunctionBegin;
  build->options = NULL;
  err = dFree(build->work);dCHK(err);
  build->workLength = 0;
  dFunctionReturn(0);
}

static dErr TensorRuleCreate(TensorBuilder *build,dInt size,dInt bytesRemaining,TensorRule *rule,dInt *bytesUsed)
{
  TensorRuleOptions *opt = (TensorRuleOptions*)build->options;
  dReal *base,*weight,*coord;
  dInt bytesNeeded;
  dErr err;

  dFunctionBegin;
  /* Determine how much output space we will use */
  base = &rule->data[0];
  weight = dNextCacheAligned(base); 
  coord = dNextCacheAligned(weight + size);
  bytesNeeded = (dInt)((size_t)(coord + size) - (size_t)rule);
  if (bytesNeeded > bytesRemaining) {
    *bytesUsed = 0;
    dFunctionReturn(0);
    /* For now, also throw an error, just to see when this fires. */
    dERROR(1,"Not enough space in buffer, %d bytes needed, %d remaining",bytesNeeded,bytesRemaining);
  }
  if (build->workLength < size) { /* make sure we have enough work space */
    err = dFree(build->work);dCHK(err);
    build->workLength = dMax(2*build->workLength,2*size);
    err = dMalloc(build->workLength*sizeof(dReal),&build->work);dCHK(err);
  }
  rule->size = size;
  rule->weightOffset = (dInt)(weight - base);
  rule->coordOffset = (dInt)(coord - base);
  /* Check options */
  if (!opt) dERROR(1,"TensorRuleOptions not set.");
  if (opt->family != GAUSS) dERROR(1,"GaussFamily %d not supported",opt->family);
  if (size == 1) {              /* Polylib function fails for this size. */
    weight[0] = 2.0;
    coord[0] = 0.0;
  } else {
    zwgj(coord,weight,size,opt->alpha,opt->beta); /* polylib function */
  }
  *bytesUsed = bytesNeeded;
  dFunctionReturn(0);
}

static dErr TensorBasisCreate(TensorBuilder *build,const TensorRule *rule,dInt basisSize,dInt bytesRemaining,TensorBasis *tbasis,dInt *bytesUsed)
{
  TensorBasisOptions *opt = (TensorBasisOptions*)build->options;
  const dInt ruleSize=rule->size,N=ruleSize*basisSize;
  dReal *base,*basis,*deriv,*node;
  dInt bytesNeeded;
  dErr err;

  dFunctionBegin;
  /* Determine how much output space we will use */
  base = &tbasis->data[0];
  basis = dNextCacheAligned(base);
  deriv = dNextCacheAligned(basis + N);
  node = dNextCacheAligned(deriv + N);
  bytesNeeded = (dInt)((size_t)(node+basisSize) - (size_t)tbasis);
  if (bytesNeeded > bytesRemaining) {
    *bytesUsed = 0;
    dFunctionReturn(0);
    dERROR(1,"Not enough space in buffer, %d bytes needed, %d remaining",bytesNeeded,bytesRemaining);
  }
  if (build->workLength < 2*N) { /* make sure we have enough work space */
    err = dFree(build->work);dCHK(err);
    build->workLength = dMax(2*build->workLength,2*N);
    err = dMalloc(build->workLength*sizeof(dReal),&build->work);dCHK(err);
  }

  tbasis->ruleSize    = ruleSize;
  tbasis->basisSize   = basisSize;
  tbasis->basisOffset = (dInt)(basis - base);
  tbasis->derivOffset = (dInt)(deriv - base);
  tbasis->nodeOffset  = (dInt)(node - base);
  /* Check options */
  if (!opt) dERROR(1,"TensorRuleOptions not set.");
  if (opt->family != GAUSS_LOBATTO) dERROR(1,"GaussFamily %d not supported",opt->family);
  if (basisSize == 1) {         /* degenerate case */
    node[0] = 0.0; basis[0] = 1.0; deriv[0] = 0.0;
    dFunctionReturn(0);
  }
  {
    const dReal alpha=opt->alpha,beta=opt->beta;
    dReal *cDeriv = build->work;
    dReal *cDerivT = build->work + dSqr(basisSize);
    dReal *coord = build->work;
    err = dMemcpy(coord,&rule->data[rule->coordOffset],rule->size);dCHK(err);
    node[0] = -1.0; node[basisSize-1] = 1.0;
    jacobz(basisSize-2,node+1,alpha+1.0,beta+1.0);
    Imglj(basis,node,coord,basisSize,ruleSize,alpha,beta);
    Dglj(cDeriv,cDerivT,node,basisSize,alpha,beta);
    for (dInt i=0; i<ruleSize; i++) {
      for (dInt j=0; j<basisSize; j++) {
        dReal z = 0;
        for (dInt k=0; k<basisSize; k++) {
          z += basis[i*basisSize+k] * deriv[k*basisSize+j];
        }
        deriv[i*basisSize+j] = z;
      }
    }
  }
  *bytesUsed = bytesNeeded;
  dFunctionReturn(0);
}

/** 
* Just an error checking indexing function.
* 
* @param this 
* @param m 
* @param out 
* 
* @return 
*/
static dErr TensorGetRule(Tensor this,dInt m,TensorRule **out)
{

  dFunctionBegin;
  if (!this->setupcalled) dERROR(1,"Attempt to get rule before Tensor setup.");
  if (!(0 < m && m < this->M)) dERROR(1,"Rule %d not less than limit %d",m,this->M);
  *out = this->rule[m];
  dFunctionReturn(0);
}

/** 
* An error checking lookup function.
* 
* @param this 
* @param m 
* @param n 
* @param out 
* 
* @return 
*/
static dErr TensorGetBasis(Tensor this,dInt m,dInt n,TensorBasis **out)
{

  dFunctionBegin;
  if (!this->setupcalled) dERROR(1,"Attempt to get basis before Tensor setup.");
  if (!(0 < m && m < this->M)) dERROR(1,"Rule size %d not less than limit %d",m,this->M);
  if (!(0 < n && n < this->N)) dERROR(1,"Basis size 5d not less than limit %d",n,this->N);
  *out = this->basis[m*this->N+n];
  dFunctionReturn(0);
}
