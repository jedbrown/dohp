#ifndef _TENSOR_H
#define _TENSOR_H

#include "petscsys.h"
#include "dohpjacobi.h"
// #include "private/dohpimpl.h"
#include "private/jacimpl.h"

PETSC_EXTERN_CXX_BEGIN

typedef enum { GAUSS, GAUSS_LOBATTO, GAUSS_RADAU } GaussFamily;

/**
* Describes a traversal of a 3D array of values.
*
* @example block {1000,{5,8,9},
*
*/
struct TensorBlock {
  dInt start;                   /**< array index for [0][0][0] */
  dInt len[3];                  /**< total number of elements in each direction */
  dInt stride[3];               /**< stride in each direction */
};

/**
* The Rule and Basis building functions need work space that is freeable later (mostly just so that valgrind won't
* complain, it's not a serious memory leak, but this design is more flexible than internal static memory management).
*
*/
typedef struct s_TensorBuilder *TensorBuilder;
struct s_TensorBuilder {
  void *options;
  dInt workLength;
  dReal *work;
};

typedef struct s_TensorRuleOptions *TensorRuleOptions;
struct s_TensorRuleOptions {
  dReal       alpha,beta;
  GaussFamily family;
};

typedef struct s_TensorRule *TensorRule;
/**
* We make no particular attempt to align the beginning of the structure, however we would like to align the arrays,
* especially the large derivative arrays, at least to 16-byte boundaries (to enable SSE) and preferably to cache line
* boundaries.
*
*/
struct s_TensorRule {
  dInt  size;                      /**< number of quadrature points */
  dReal *weight,*coord;            /**< weights and nodal coordinates */
};

typedef dErr (*TensorMultFunction)(dInt,const dInt*,const dInt*,const dReal**,const dScalar*,dScalar*restrict,InsertMode);

typedef struct s_TensorBasis *TensorBasis;
/**
* Stores a one-dimensional part of a Tensor product basis.  Includes optimized functions for 3D tensor multiplication.
*
*/
struct s_TensorBasis {
  TensorMultFunction mult;
  dInt  P,Q;
  dReal *interp,*deriv,*node;
  dReal *interpTranspose,*derivTranspose;
  const dReal *mscale,*lscale;        /**< Used to produce optimal scaling of sparse mass and Laplacian matrices */
};

typedef struct s_TensorBasisOptions *TensorBasisOptions;
struct s_TensorBasisOptions {
  dReal       alpha,beta;
  GaussFamily family;
};


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
  TensorBuilder ruleBuilder,basisBuilder;
  TensorRuleOptions ruleOpts;
  TensorBasisOptions basisOpts;
  TensorRule  *rule;   /**< Array of 1D rules.  The rule with \a m points is indexed as \a rule[\a m]. */
  TensorBasis *basis;  /**< Array of 1D bases, length \f$ M*N \f$.  Consider basis with \a m quadrature
                                  * points and \a n basis functions.
                                  *
                                  * If the quadrature points are \f$ q_i, i = 0,1,\dotsc,m-1 \f$
                                  *
                                  * and the nodes are \f$ x_j, j=0,1,\dotsc,n-1 \f$
                                  *
                                  * then \f$ f(q_i) = \sum_{j=0}^n B[i*n+j] f(x_j) \f$ */
  dInt M;                       /**< number of rules */
  dInt N;                       /**< basis size limit */
  dTruth usemscale,uselscale;
  /* dBufferList data; */
  struct _dRuleOps *ruleOpsLine,*ruleOpsQuad,*ruleOpsHex;
  struct _dEFSOps *efsOpsLine,*efsOpsQuad,*efsOpsHex;
  dBool setupcalled;
};

typedef struct {
  dRuleHEADER;
  TensorRule trule[3];
} dRule_Tensor;

typedef struct {
  dEFSHEADER;
  TensorBasis basis[3];
} dEFS_Tensor;

EXTERN dErr dJacobiRuleOpsSetUp_Tensor(dJacobi jac);
EXTERN dErr dJacobiEFSOpsSetUp_Tensor(dJacobi jac);
EXTERN dErr dJacobiRuleOpsDestroy_Tensor(dJacobi jac);
EXTERN dErr dJacobiEFSOpsDestroy_Tensor(dJacobi jac);

EXTERN dErr TensorRuleView(const TensorRule,PetscViewer);
EXTERN dErr TensorBasisView(const TensorBasis,PetscViewer);

PETSC_EXTERN_CXX_END

#endif  /* _TENSOR_H */
