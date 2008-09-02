#include "tensor.h"

static dErr dEFSView_Tensor_Private(const char *,dRule*,dInt,TensorBasis*,PetscViewer);
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
#define _F(f) static dErr f(dEFS*,dInt*,dInt*,dReal**restrict) /* dEFSGetTensorNodes */
_F(dEFSGetTensorNodes_Tensor_Line);
_F(dEFSGetTensorNodes_Tensor_Quad);
_F(dEFSGetTensorNodes_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS*,dInt,dInt*,dScalar**restrict,const dScalar[],dScalar[],dApplyMode,InsertMode) /* dEFSApply */
//_F(dEFSApply_Tensor_Line);
//_F(dEFSApply_Tensor_Quad);
_F(dEFSApply_Tensor_Hex);
#undef _F
#if 0
# define _F(f) static dErr f(dEFS*,dInt,dInt,const dScalar[],dScalar[],InsertMode,ScatterMode) /* dEFSScatterInt */
_F(dEFSScatterInt_Tensor_Line);
_F(dEFSScatterInt_Tensor_Quad);
_F(dEFSScatterInt_Tensor_Hex);
# undef _F
# define _F(f) static dErr f(dEFS,dEFS,dInt*,dScalar**restrict,const dScalar[],dScalar[],InsertMode,ScatterMode) /* dEFSScatterFacet */
_F(dEFSScatterFacet_Tensor_Line);
_F(dEFSScatterFacet_Tensor_Quad);
_F(dEFSScatterFacet_Tensor_Hex);
# undef _F
#endif
#define _F(f) static dErr f(dInt D,const dInt P[],const dInt Q[],dInt *wlen,dScalar **restrict work,dReal *A[],dTransposeMode tpose[],const dScalar f[],dScalar g[restrict],InsertMode imode)
//_F(TensorMult_Line);
//_F(TensorMult_Quad);
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
                                               .getTensorNodes = dEFSGetTensorNodes_Tensor_Line,
                                               .apply = 0, /* dEFSApply_Tensor_Line, */
                                               .scatterInt = 0,
                                               .scatterFacet = 0 };
  static const struct v_dEFSOps efsOpsQuad = { .view = dEFSView_Tensor_Quad,
                                               .getSizes = dEFSGetSizes_Tensor_Quad,
                                               .getTensorNodes = dEFSGetTensorNodes_Tensor_Quad,
                                               .apply = 0, /* dEFSApply_Tensor_Quad, */
                                               .scatterInt = 0,
                                               .scatterFacet = 0 };
  static const struct v_dEFSOps efsOpsHex  = { .view = dEFSView_Tensor_Hex,
                                               .getSizes = dEFSGetSizes_Tensor_Hex,
                                               .getTensorNodes = dEFSGetTensorNodes_Tensor_Hex,
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

dErr dJacobiEFSOpsDestroy_Tensor(dJacobi jac)
{
  Tensor this = (Tensor)jac->impl;
  dErr err;

  dFunctionBegin;
  if (this->efsOpsLine) { err = dFree(this->efsOpsLine);dCHK(err); }
  if (this->efsOpsQuad) { err = dFree(this->efsOpsQuad);dCHK(err); }
  if (this->efsOpsHex)  { err = dFree(this->efsOpsHex);dCHK(err); }
  dFunctionReturn(0);
}


static dErr dEFSView_Tensor_Private(const char *name,dRule *rule,dInt n,TensorBasis *b,PetscViewer viewer)
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (ascii) {
    err = PetscViewerASCIIPrintf(viewer,"dEFS type %s\n",name);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"based on dRule:\n");dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = dRuleView(rule,viewer);dCHK(err);
    err = PetscViewerASCIIPopTab(viewer);
    for (dInt i=0; i<n; i++) {
      err = TensorBasisView(b[i],viewer);
    }
  }
  dFunctionReturn(0);
}


static dErr dEFSView_Tensor_Line(dEFS *efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Tensor_Private("Tensor_Line",efs->rule,1,((struct s_dEFS_Tensor_Line*)efs->data)->basis,viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr dEFSView_Tensor_Quad(dEFS *efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Tensor_Private("Tensor_Quad",efs->rule,2,((struct s_dEFS_Tensor_Quad*)efs->data)->basis,viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr dEFSView_Tensor_Hex(dEFS *efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Tensor_Private("Tensor_Hex",efs->rule,3,((struct s_dEFS_Tensor_Hex*)efs->data)->basis,viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dEFSGetSizes_Tensor_Line(dEFS *efs,dInt *dim,dInt *inodes,dInt *total)
{
  TensorBasis *b = ((struct s_dEFS_Tensor_Line*)efs->data)->basis;

  dFunctionBegin;
  if (dim) *dim = 1;
  if (inodes) *inodes = b[0]->P - 2;
  if (total) *total = b[0]->P;
  dFunctionReturn(0);
}
static dErr dEFSGetSizes_Tensor_Quad(dEFS *efs,dInt *dim,dInt *inodes,dInt *total)
{
  TensorBasis *b = ((struct s_dEFS_Tensor_Quad*)efs->data)->basis;

  dFunctionBegin;
  if (dim) *dim = 2;
  if (inodes) *inodes = (b[0]->P - 2) * (b[1]->P - 2);
  if (total) *total = b[0]->P * b[1]->P;
  dFunctionReturn(0);
}
static dErr dEFSGetSizes_Tensor_Hex(dEFS *efs,dInt *dim,dInt *inodes,dInt *total)
{
  TensorBasis *b = ((struct s_dEFS_Tensor_Hex*)efs->data)->basis;

  dFunctionBegin;
  if (dim) *dim = 3;
  if (inodes) *inodes = (b[0]->P - 2) * (b[1]->P - 2) * (b[2]->P - 2);
  if (total) *total = b[0]->P * b[1]->P * b[2]->P;
  dFunctionReturn(0);
}

static dErr dEFSGetTensorNodes_Tensor_Line(dEFS *efs,dInt *dim,dInt tsize[restrict],dReal *x[restrict])
{
  TensorBasis *b = ((struct s_dEFS_Tensor_Line*)efs->data)->basis;

  dFunctionBegin;
  if (dim) *dim = 1;
  if (tsize) tsize[0] = b[0]->P;
  x[0] = b[0]->node;
  dFunctionReturn(0);
}

static dErr dEFSGetTensorNodes_Tensor_Quad(dEFS *efs,dInt *dim,dInt tsize[restrict],dReal *x[restrict])
{
  TensorBasis *b = ((struct s_dEFS_Tensor_Quad*)efs->data)->basis;

  dFunctionBegin;
  if (dim) *dim = 2;
  for (dInt i=0; i<2; i++) {
    if (tsize) tsize[i] = b[i]->P;
    if (x) x[i] = b[i]->node;
  }
  dFunctionReturn(0);
}

static dErr dEFSGetTensorNodes_Tensor_Hex(dEFS *efs,dInt *dim,dInt tsize[restrict],dReal *x[restrict])
{
  TensorBasis *b = ((struct s_dEFS_Tensor_Hex*)efs->data)->basis;

  dFunctionBegin;
  if (dim) *dim = 3;
  for (dInt i=0; i<3; i++) {
    if (tsize) tsize[i] = b[i]->P;
    if (x) x[i] = b[i]->node;
  }
  dFunctionReturn(0);
}

static dErr dEFSApply_Tensor_Hex(dEFS *efs,dInt D,dInt *wlen,dScalar**restrict work,const dScalar in[],dScalar out[],dApplyMode amode,InsertMode imode)
{
  TensorBasis *b = ((struct s_dEFS_Tensor_Hex*)efs->data)->basis;
  /* const dInt *P=&b->P,*Q=&b->Q; */
  dScalar *A[3];
  dInt P[3],Q[3],chunk;
  dTransposeMode tpose[3];
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<3; i++) {
    P[i] = b[i]->P;
    Q[i] = b[i]->Q;
  }
  switch (amode) {
    case dAPPLY_INTERP:
      for (dInt i=0; i<3; i++) {
        A[i] = b[i]->interp;
        tpose[i] = dTRANSPOSE_NO;
      }
      err = TensorMult_Hex(D,P,Q,wlen,work,A,tpose,in,out,imode);dCHK(err);
      break;
    case dAPPLY_INTERP_TRANSPOSE:
      for (dInt i=0; i<3; i++) {
        A[i] = b[i]->interp;
        tpose[i] = dTRANSPOSE_YES;
      }
      err = TensorMult_Hex(D,Q,P,wlen,work,A,tpose,in,out,imode);dCHK(err);
      break;
    case dAPPLY_GRAD:
      chunk = D*Q[0]*Q[1]*Q[2]; /* number of dofs for each component of the gradient */
      for (dInt i=0; i<3; i++) {
        for (dInt j=0; j<3; j++) {
          if (i == j) A[j] = b[j]->deriv;
          else A[j] = b[j]->interp;
          tpose[j] = dTRANSPOSE_NO;
        }
        err = TensorMult_Hex(D,P,Q,wlen,work,A,tpose,in,&out[i*chunk],imode);dCHK(err);
      }
      break;
    case dAPPLY_GRAD_TRANSPOSE:
    default:
      dERROR(1,"invalid dApplyMode %d specified",amode);
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
static dErr TensorMult_Hex(dInt D,const dInt P[3],const dInt Q[3],dInt *wlen,dScalar **restrict work,
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
      for (j=0; j<P[1]; j++) {
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
