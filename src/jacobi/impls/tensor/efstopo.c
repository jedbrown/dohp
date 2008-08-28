#include "tensor.h"

#ifdef _F
# undef _F
#endif
#define _F(f) static dErr f(dEFS*,PetscViewer) /* dEFSView */
_F(dEFSView_Tensor_Line);
_F(dEFSView_Tensor_Quad);
_F(dEFSView_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS*,dInt*,dInt*,dInt*) /* dEFSGetSizes */
_F(dEFSGetSizes_Tensor_Line);
_F(dEFSGetSizes_Tensor_Quad);
_F(dEFSGetSizes_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS*,dInt,dInt*,dScalar**restrict,const dScalar[],dScalar[],dApplyType,InsertMode) /* dEFSApply */
_F(dEFSApply_Tensor_Line);
_F(dEFSApply_Tensor_Quad);
_F(dEFSApply_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS*,dInt,dInt,const dScalar[],dScalar[],InsertMode,ScatterMode) /* dEFSScatterInt */
_F(dEFSScatterInt_Tensor_Line);
_F(dEFSScatterInt_Tensor_Quad);
_F(dEFSScatterInt_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS,dEFS,dInt*,dScalar**restrict,const dScalar[],dScalar[],InsertMode,ScatterMode) /* dEFSScatterFacet */
_F(dEFSScatterFacet_Tensor_Line);
_F(dEFSScatterFacet_Tensor_Quad);
_F(dEFSScatterFacet_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dInt D,const dInt P[],const dInt Q[],dInt *wlen,dScalar *restrict* work,dReal *A[],dTransposeMode tpose[],const dScalar f[],dScalar g[restrict],InsertMode imode)
_F(TensorMult_Line);
_F(TensorMult_Quad);
_F(TensorMult_Hex);
#undef _F


/** 
* Set up the EFS ops table for each topology.  This is the only exported function in this file.
* 
* @param jac 
* 
* @return 
*/
dErr dJacobiEFSOpsSetUp_Tensor(dJacobi jac)
{
  static const struct v_dEFSOps efsOpsLine = { .view = dEFSView_Tensor_Line,
                                               .getSizes = dEFSGetSizes_Tensor_Line,
                                               .apply = dEFSApply_Tensor_Line,
                                               .scatterInt = 0,
                                               .scatterFacet = 0 };
  static const struct v_dEFSOps efsOpsQuad = { .view = dEFSView_Tensor_Quad,
                                               .getSizes = dEFSGetSizes_Tensor_Quad,
                                               .apply = dEFSApply_Tensor_Quad,
                                               .scatterInt = 0,
                                               .scatterFacet = 0 };
  static const struct v_dEFSOps efsOpsHex  = { .view = dEFSView_Tensor_Hex,
                                               .getSizes = dEFSGetSizes_Tensor_Hex,
                                               .apply = dEFSApply_Tensor_Hex,
                                               .scatterInt = 0,
                                               .scatterFacet = 0 };
  Tensor this = (Tensor)jac->impl;
  dErr err;

  dFunctionBegin;
  if (!this->efsOpsLine) {
    err = dMalloc(sizeof(struct v_dEFSOps),&this->efsOpsLine);dCHK(err);
    err = dMemcpy(this->efsOpsLine,&efsOpsLine,sizeof(struct v_dEFSOps));
  }
  if (!this->efsOpsQuad) {
    err = dMalloc(sizeof(struct v_dEFSOps),&this->efsOpsQuad);dCHK(err);
    err = dMemcpy(this->efsOpsQuad,&efsOpsQuad,sizeof(struct v_dEFSOps));
  }
  if (!this->efsOpsHex) {
    err = dMalloc(sizeof(struct v_dEFSOps),&this->efsOpsHex);dCHK(err);
    err = dMemcpy(this->efsOpsHex,&efsOpsHex,sizeof(struct v_dEFSOps));
  }
  dFunctionReturn(0);
}


/** 
* The core computational kernel.  Performs a tensor product operation with the matrices A,B,C.
* 
* @param[in] D number of degrees of freedom
* @param[in] P array of length 3, dimensions of input block
* @param[in] Q array of length 3, dimensions of output block
* @param[in] tpose dTRANSPOSE_YES if the arrays \a A, \a B, \a C are transposed with respect to the input and output spaces
* @param[in,out] wlen length of \a work array
* @param[in,out] work pointer to workspace for use by this function, must have been allocated with PetscMalloc.  This function will reallocate if it is not sufficient.
* @param[in] A array of pointers, length 3, transformation matrix in each direction, shape(A[i])=(Q[i],P[i]), i=0,1,2
* @param[in] f input vector with size corresponding to \in
* @param[in,out] g output vector with size corresponding to \out
* @param[in] imode ADD_VALUES or INSERT_VALUES
* 
* @return err
*/
static dErr TensorMult_Hex(dInt D,const dInt P[3],const dInt Q[3],dInt *wlen,dScalar *restrict* work,
                           dReal *A[3],dTransposeMode tpose[3],const dScalar f[],dScalar g[restrict],InsertMode imode)
{
  dInt i,j,k,l,d,idx;
  dReal *B[3];
  dScalar *restrict a,*restrict b;
  dErr err;

  dFunctionBegin;
  idx = 0;
  do {                          /* This loop will execute once if there is enough work space, otherwise it will execute
                                * a second time, allocating enough memory at the beginning. */
    if (idx > *wlen) {
      err = dFree(*work);dCHK(err);
      *wlen = idx*3/2;
      err = dMalloc((*wlen)*sizeof(dScalar),work);dCHK(err);
      idx = 0;
    }
    for (i=0; i<3; i++) {
      switch (tpose[i]) {
        case dTRANSPOSE_NO:
          B[i] = A[i];
          break;
        case dTRANSPOSE_YES:
          B[i] = &(*work)[idx]; idx += Q[i]*P[i];
          if (idx < *wlen) {    /* If we ran out of work space, the 'do {} while' loop will repeat. */
            for (j=0; j<Q[i]; j++) {
              for (k=0; k<P[i]; k++) {
                B[i][j*P[i]+k] = A[i][k*Q[i]+j];
              }
            }
          }
          break;
      }
    };
    a = &(*work)[idx]; idx += Q[0]*P[1]*P[2]*D;
    b = &(*work)[idx]; idx += Q[0]*Q[1]*P[2]*D;
  } while (idx > *wlen);

  err = dMemzero(*work,idx*sizeof(a[0]));dCHK(err);
  switch (imode) {
    case INSERT_VALUES:
      err = dMemzero(g,Q[0]*Q[1]*Q[2]*D);dCHK(err);
      break;
    case ADD_VALUES:
      break;
    default:
      dERROR(1,"Requested InsertMode %d not supported for this operation.",imode);
  }

  for (l=0; l<Q[0]; l++) {
    for (i=0; i<P[0]; i++) {
      for (j=0; j<P[1]; i++) {
        for (k=0; k<P[2]; k++) {
          for (d=0; d<D; d++) {
            a[((l*P[1]+j)*P[2]+k)*D+d] += B[0][l*P[0]+i] * f[((i*P[1]+j)*P[2]+k)*D+d];
          }
        }
      }
    }
  }
  for (i=0; i<Q[0]; i++) {
    for (l=0; l<Q[1]; l++) {
      for (j=0; j<P[1]; j++) {
        for (k=0; k<P[2]; k++) {
          for (d=0; d<D; d++) {
            b[((i*Q[1]+l)*P[2]+k)*D+d] += B[1][l*P[1]+j] * a[((i*P[1]+j)*P[2]+k)*D+d];
          }
        }
      }
    }
  }
  for (i=0; i<Q[0]; i++) {
    for (j=0; j<Q[1]; j++) {
      for (l=0; l<Q[2]; l++) {
        for (k=0; k<P[2]; k++) {
          for (d=0; d<D; d++) {
            g[((i*Q[1]+j)*Q[2]+l)*D+d] += B[2][l*P[2]+k] * b[((i*Q[1]+j)*P[2]+k)*D+d];
          }
        }
      }
    }
  }
  dFunctionReturn(0);
}
