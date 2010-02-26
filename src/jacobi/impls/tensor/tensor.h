#ifndef _TENSOR_H
#define _TENSOR_H

#include <dohpkhash.h>
#include <dohpjacimpl.h>

dEXTERN_C_BEGIN

typedef dErr (*TensorMultFunction)(dInt,const dInt[3],const dInt[3],const dReal*[],const dScalar*,dScalar[restrict],InsertMode);

typedef struct s_TensorBasis *TensorBasis;
/**
* Stores a one-dimensional part of a Tensor product basis.  Includes optimized functions for 3D tensor multiplication.
*
*/
struct s_TensorBasis {
  TensorMultFunction multhex[3];
  dInt  P,Q;
  dReal *interp,*deriv,*node;
  dReal *interpTranspose,*derivTranspose;
  dReal *weight;                /**< Weight for quadrature rule at interpolation nodes */
  const dReal *mscale,*lscale;  /**< Used to produce optimal scaling of sparse mass and Laplacian matrices */
};

KHASH_MAP_INIT_INT(basis,TensorBasis)

typedef struct s_Tensor *Tensor;
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
struct s_Tensor {
  dReal       alpha,beta;
  GaussFamily family;
  khash_t(basis) *bases;  /**< Hash of bases.  Consider basis with \a m quadrature
                                  * points and \a n basis functions.
                                  *
                                  * If the quadrature points are \f$ q_i, i = 0,1,\dotsc,m-1 \f$
                                  *
                                  * and the nodes are \f$ x_j, j=0,1,\dotsc,n-1 \f$
                                  *
                                  * then \f$ f(q_i) = \sum_{j=0}^n B[i*n+j] f(x_j) \f$ */
  dInt M;                       /**< number of rules */
  dInt N;                       /**< basis size limit */
  PetscTruth usemscale,uselscale,nounroll;
  /* dBufferList data; */
  struct _dEFSOps *efsOpsLine,*efsOpsQuad,*efsOpsHex;
  bool setupcalled;
};

typedef struct {
  dEFSHEADER;
  TensorBasis basis[3];
} dEFS_Tensor;

extern dErr dJacobiEFSOpsSetUp_Tensor(dJacobi jac);
extern dErr dJacobiEFSOpsDestroy_Tensor(dJacobi jac);

extern dErr TensorBasisView(const TensorBasis,PetscViewer);

dEXTERN_C_END

#endif  /* _TENSOR_H */
