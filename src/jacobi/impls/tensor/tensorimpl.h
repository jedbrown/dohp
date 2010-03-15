#ifndef _TENSOR_H
#define _TENSOR_H

#include <dohpkhash.h>
#include <dohpjacimpl.h>

dEXTERN_C_BEGIN

typedef dErr (*TensorMultFunction)(dInt,const dInt[3],const dInt[3],const dReal*[],const dScalar*,dScalar[restrict],InsertMode);

typedef struct _TensorBasis *TensorBasis;
/**
* Stores a one-dimensional part of a Tensor product basis.  Includes optimized functions for 3D tensor multiplication.
*
*/
struct _TensorBasis {
  TensorMultFunction multhex[3];
  dInt  P,Q;
  dReal *interp,*deriv,*node;
  dReal *interpTranspose,*derivTranspose;
  dReal *weight;                /**< Weight for quadrature rule at interpolation nodes */
  const dReal *mscale,*lscale;  /**< Used to produce optimal scaling of sparse mass and Laplacian matrices */
};

typedef struct {
  dEFSHEADER;
  dEntTopology topo;
  TensorBasis  basis[3];
} dEFS_Tensor;

KHASH_MAP_INIT_INT(tensor,TensorBasis)

typedef struct {
  dEntTopology topo;
  dPolynomialOrder degree;
  dRule rule;
} khu_efskey_t;
static inline khint_t khu_efskey_hash_func(khu_efskey_t key)
{ return kh_int_hash_func((khint32_t)key.topo) ^ kh_int_hash_func((khint32_t)key.degree) ^ kh_int64_hash_func((khint64_t)(uintptr_t)key.rule); }
static inline bool khu_efskey_hash_equal(khu_efskey_t a,khu_efskey_t b)
{ return (a.topo == b.topo) && (a.degree == b.degree) && (a.rule == b.rule); }
KHASH_INIT(efs, khu_efskey_t, dEFS_Tensor*, 1, khu_efskey_hash_func, khu_efskey_hash_equal)

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
*/
typedef struct {
  dReal       alpha,beta;
  dGaussFamily family;
  khash_t(tensor) *tensor;  /**< Hash of bases.  Consider basis with \a m quadrature
                                  * points and \a n basis functions.
                                  *
                                  * If the quadrature points are \f$ q_i, i = 0,1,\dotsc,m-1 \f$
                                  *
                                  * and the nodes are \f$ x_j, j=0,1,\dotsc,n-1 \f$
                                  *
                                  * then \f$ f(q_i) = \sum_{j=0}^n B[i*n+j] f(x_j) \f$ */
  khash_t(efs) *efs;
  PetscTruth usemscale,uselscale,nounroll;
  /* dBufferList data; */
  struct _dEFSOps *efsOpsLine,*efsOpsQuad,*efsOpsHex;
} dJacobi_Tensor;

extern dErr dJacobiEFSOpsSetUp_Tensor(dJacobi jac);
extern dErr dJacobiEFSOpsDestroy_Tensor(dJacobi jac);

extern dErr TensorBasisView(const TensorBasis,PetscViewer);

dEXTERN_C_END

#endif  /* _TENSOR_H */
