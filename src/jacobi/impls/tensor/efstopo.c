#include "tensor.h"
#include "dohpmesh.h"
#include "dohpgeom.h"

static dErr dEFSView_Tensor_Private(const char *,dRule,dInt,TensorBasis*,PetscViewer);
#ifdef _F
# undef _F
#endif
#define _F(f) static dErr f(dEFS,PetscViewer) /* dEFSView */
_F(dEFSView_Tensor_Line);
_F(dEFSView_Tensor_Quad);
_F(dEFSView_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS,dInt*,dInt*,dInt*) /* dEFSGetSizes */
_F(dEFSGetSizes_Tensor_Line);
_F(dEFSGetSizes_Tensor_Quad);
_F(dEFSGetSizes_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS,dInt*,dInt*,dReal**restrict) /* dEFSGetTensorNodes */
_F(dEFSGetTensorNodes_Tensor_Line);
_F(dEFSGetTensorNodes_Tensor_Quad);
_F(dEFSGetTensorNodes_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS,const dReal[],dInt,const dScalar[],dScalar[restrict],dApplyMode,InsertMode) /* dEFSApply */
_F(dEFSApply_Tensor_Line);
_F(dEFSApply_Tensor_Quad);
_F(dEFSApply_Tensor_Hex);
#undef _F
#if 0
# define _F(f) static dErr f(dEFS,dInt,dInt,const dScalar[],dScalar[],InsertMode,ScatterMode) /* dEFSScatterInt */
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

#define _F(f) static dErr f(dInt D,const dInt P[],const dInt Q[],dReal *A[],const dScalar f[],dScalar g[restrict],InsertMode imode)
_F(TensorMult_Line);
_F(TensorMult_Quad);
_F(TensorMult_Hex);
#undef _F
#define _F(f) static dErr f(dRule_Tensor*,const dReal[],dInt,const dScalar[],dScalar[],InsertMode)
_F(dRuleMappingApply_Tensor_Line);
_F(dRuleMappingApply_Tensor_Quad);
_F(dRuleMappingApply_Tensor_Hex);
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
  static const struct _dEFSOps efsOpsLine = { .view = dEFSView_Tensor_Line,
                                              .getSizes = dEFSGetSizes_Tensor_Line,
                                              .getTensorNodes = dEFSGetTensorNodes_Tensor_Line,
                                              .apply = dEFSApply_Tensor_Line,
                                              .scatterInt = 0,
                                              .scatterFacet = 0 };
  static const struct _dEFSOps efsOpsQuad = { .view = dEFSView_Tensor_Quad,
                                              .getSizes = dEFSGetSizes_Tensor_Quad,
                                              .getTensorNodes = dEFSGetTensorNodes_Tensor_Quad,
                                              .apply = dEFSApply_Tensor_Quad,
                                              .scatterInt = 0,
                                              .scatterFacet = 0 };
  static const struct _dEFSOps efsOpsHex  = { .view = dEFSView_Tensor_Hex,
                                              .getSizes = dEFSGetSizes_Tensor_Hex,
                                              .getTensorNodes = dEFSGetTensorNodes_Tensor_Hex,
                                              .apply = dEFSApply_Tensor_Hex,
                                              .scatterInt = 0,
                                              .scatterFacet = 0 };
  Tensor this = (Tensor)jac->impl;
  dErr err;

  dFunctionBegin;
  if (!this->efsOpsLine) {
    err = dMalloc(sizeof(struct _dEFSOps),&this->efsOpsLine);dCHK(err);
    err = dMemcpy(this->efsOpsLine,&efsOpsLine,sizeof(struct _dEFSOps));
  }
  if (!this->efsOpsQuad) {
    err = dMalloc(sizeof(struct _dEFSOps),&this->efsOpsQuad);dCHK(err);
    err = dMemcpy(this->efsOpsQuad,&efsOpsQuad,sizeof(struct _dEFSOps));
  }
  if (!this->efsOpsHex) {
    err = dMalloc(sizeof(struct _dEFSOps),&this->efsOpsHex);dCHK(err);
    err = dMemcpy(this->efsOpsHex,&efsOpsHex,sizeof(struct _dEFSOps));
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


static dErr dEFSView_Tensor_Private(const char *name,dRule rule,dInt n,TensorBasis *b,PetscViewer viewer)
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


static dErr dEFSView_Tensor_Line(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Tensor_Private("Tensor_Line",efs->rule,1,((dEFS_Tensor*)efs)->basis,viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr dEFSView_Tensor_Quad(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Tensor_Private("Tensor_Quad",efs->rule,2,((dEFS_Tensor*)efs)->basis,viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr dEFSView_Tensor_Hex(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Tensor_Private("Tensor_Hex",efs->rule,3,((dEFS_Tensor*)efs)->basis,viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dEFSGetSizes_Tensor_Line(dEFS efs,dInt *dim,dInt *inodes,dInt *total)
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;

  dFunctionBegin;
  if (dim) *dim = 1;
  if (inodes) *inodes = b[0]->P - 2;
  if (total) *total = b[0]->P;
  dFunctionReturn(0);
}
static dErr dEFSGetSizes_Tensor_Quad(dEFS efs,dInt *dim,dInt *inodes,dInt *total)
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;

  dFunctionBegin;
  if (dim) *dim = 2;
  if (inodes) *inodes = (b[0]->P - 2) * (b[1]->P - 2);
  if (total) *total = b[0]->P * b[1]->P;
  dFunctionReturn(0);
}
static dErr dEFSGetSizes_Tensor_Hex(dEFS efs,dInt *dim,dInt *inodes,dInt *total)
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;

  dFunctionBegin;
  if (dim) *dim = 3;
  if (inodes) *inodes = (b[0]->P - 2) * (b[1]->P - 2) * (b[2]->P - 2);
  if (total) *total = b[0]->P * b[1]->P * b[2]->P;
  dFunctionReturn(0);
}

static dErr dEFSGetTensorNodes_Tensor_Line(dEFS efs,dInt *dim,dInt tsize[restrict],dReal *x[restrict])
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;

  dFunctionBegin;
  if (dim) *dim = 1;
  if (tsize) {
    tsize[0] = b[0]->P;
    tsize[1] = 1;
    tsize[2] = 1;
  }
  if (x) {
    x[0] = b[0]->node;
    x[1] = NULL;
    x[2] = NULL;
  }
  dFunctionReturn(0);
}

static dErr dEFSGetTensorNodes_Tensor_Quad(dEFS efs,dInt *dim,dInt tsize[restrict],dReal *x[restrict])
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;

  dFunctionBegin;
  if (dim) *dim = 2;
  if (tsize) {
    tsize[0] = b[0]->P;
    tsize[1] = b[1]->P;
    tsize[2] = 1;
  }
  if (x) {
    x[0] = b[0]->node;
    x[1] = b[1]->node;
    x[2] = NULL;
  }
  dFunctionReturn(0);
}

static dErr dEFSGetTensorNodes_Tensor_Hex(dEFS efs,dInt *dim,dInt tsize[restrict],dReal *x[restrict])
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;

  dFunctionBegin;
  if (dim) *dim = 3;
  if (tsize) {
    tsize[0] = b[0]->P;
    tsize[1] = b[1]->P;
    tsize[2] = b[2]->P;
  }
  if (x) {
    x[0] = b[0]->node;
    x[1] = b[1]->node;
    x[2] = b[2]->node;
  }
  dFunctionReturn(0);
}

static dErr dEFSApply_Tensor_Line(dEFS efs,const dReal mapdata[],dInt D,const dScalar in[],dScalar out[],dApplyMode amode,InsertMode imode)
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;
  dReal *A[1];
  dInt P[1],Q[1];
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<1; i++) {
    P[i] = b[i]->P;
    Q[i] = b[i]->Q;
  }
  switch (amode) {
    case dAPPLY_INTERP:
      A[0] = b[0]->interp;
      err = TensorMult_Line(D,P,Q,A,in,out,imode);dCHK(err);
      break;
    case dAPPLY_INTERP_TRANSPOSE:
      A[0] = b[0]->interpTranspose;
      err = TensorMult_Line(D,Q,P,A,in,out,imode);dCHK(err);
      break;
    case dAPPLY_GRAD: {
      dScalar df[1][Q[0]*D];
      A[0] = b[0]->deriv;
      err = TensorMult_Line(D,P,Q,A,in,df[0],INSERT_VALUES);dCHK(err);
      err = dRuleMappingApply_Tensor_Line((dRule_Tensor*)efs->rule,mapdata,D,&df[0][0],out,imode);dCHK(err);
    } break;
    case dAPPLY_GRAD_TRANSPOSE:
    case dAPPLY_SYMGRAD:
    case dAPPLY_SYMGRAD_TRANSPOSE:
    default:
      dERROR(1,"invalid dApplyMode %d specified",amode);
  }
  dFunctionReturn(0);
}

static dErr dEFSApply_Tensor_Quad(dEFS efs,const dReal mapdata[],dInt D,const dScalar in[],dScalar out[],dApplyMode amode,InsertMode imode)
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;
  dReal *A[2];
  dInt P[2],Q[2];
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<2; i++) {
    P[i] = b[i]->P;
    Q[i] = b[i]->Q;
  }
  switch (amode) {
    case dAPPLY_INTERP:
      A[0] = b[0]->interp;
      A[1] = b[1]->interp;
      err = TensorMult_Quad(D,P,Q,A,in,out,imode);dCHK(err);
      break;
    case dAPPLY_INTERP_TRANSPOSE:
      A[0] = b[0]->interpTranspose;
      A[1] = b[1]->interpTranspose;
      err = TensorMult_Quad(D,Q,P,A,in,out,imode);dCHK(err);
      break;
    case dAPPLY_GRAD: {
      dScalar df[2][Q[0]*Q[1]*D];
      A[0] = b[0]->deriv; A[1] = b[1]->interp;
      err = TensorMult_Quad(D,P,Q,A,in,df[0],INSERT_VALUES);dCHK(err);
      A[0] = b[0]->interp; A[1] = b[1]->deriv;
      err = TensorMult_Quad(D,P,Q,A,in,df[1],INSERT_VALUES);dCHK(err);
      err = dRuleMappingApply_Tensor_Quad((dRule_Tensor*)efs->rule,mapdata,D,&df[0][0],out,imode);dCHK(err);
    } break;
    case dAPPLY_GRAD_TRANSPOSE:
    default:
      dERROR(1,"invalid dApplyMode %d specified",amode);
  }
  dFunctionReturn(0);
}

static dErr dEFSApply_Tensor_Hex(dEFS efs,const dReal mapdata[],dInt D,const dScalar in[],dScalar out[],dApplyMode amode,InsertMode imode)
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;
  dReal *A[3];
  dInt P[3],Q[3];
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<3; i++) {
    P[i] = b[i]->P;
    Q[i] = b[i]->Q;
  }
  switch (amode) {
    case dAPPLY_INTERP:
      A[0] = b[0]->interp;
      A[1] = b[1]->interp;
      A[2] = b[2]->interp;
      err = TensorMult_Hex(D,P,Q,A,in,out,imode);dCHK(err);
      break;
    case dAPPLY_INTERP_TRANSPOSE:
      A[0] = b[0]->interpTranspose;
      A[1] = b[1]->interpTranspose;
      A[2] = b[2]->interpTranspose;
      err = TensorMult_Hex(D,Q,P,A,in,out,imode);dCHK(err);
      break;
    case dAPPLY_GRAD: {
      dScalar df[3][Q[0]*Q[1]*Q[2]*D];
      A[0] = b[0]->deriv; A[1] = b[1]->interp; A[2] = b[2]->interp;
      err = TensorMult_Hex(D,P,Q,A,in,df[0],INSERT_VALUES);dCHK(err);
      A[0] = b[0]->interp; A[1] = b[1]->deriv; A[2] = b[2]->interp;
      err = TensorMult_Hex(D,P,Q,A,in,df[1],INSERT_VALUES);dCHK(err);
      A[0] = b[0]->interp; A[1] = b[1]->interp; A[2] = b[2]->deriv;
      err = TensorMult_Hex(D,P,Q,A,in,df[2],INSERT_VALUES);dCHK(err);
      err = dRuleMappingApply_Tensor_Hex((dRule_Tensor*)efs->rule,mapdata,D,&df[0][0],out,imode);dCHK(err);
    } break;
    case dAPPLY_GRAD_TRANSPOSE:
    default:
      dERROR(1,"invalid/unimplemented dApplyMode %d specified",amode);
  }
  dFunctionReturn(0);
}

static dErr TensorMult_Line(dInt D,const dInt P[1],const dInt Q[1],dReal *A[1],const dScalar f[],dScalar g[restrict],InsertMode imode)
{
  dInt i,l,d;
  dErr err;

  dFunctionBegin;
  switch (imode) {
    case INSERT_VALUES:
      err = dMemzero(g,Q[0]*D*sizeof(g[0]));dCHK(err);
      break;
    case ADD_VALUES:
      break;
    default:
      dERROR(1,"Requested InsertMode %d not supported for this operation.",imode);
  }

  for (l=0; l<Q[0]; l++) {
    for (i=0; i<P[0]; i++) {
      for (d=0; d<D; d++) {
        g[l*D+d] += A[0][l*P[0]+i] * f[i*D+d];
      }
    }
  }
  dFunctionReturn(0);
}

static dErr TensorMult_Quad(dInt D,const dInt P[2],const dInt Q[2],dReal *A[2],const dScalar f[],dScalar g[restrict],InsertMode imode)
{
  dScalar a[Q[0]*P[1]*D];
  dInt i,j,l,d;
  dErr err;

  dFunctionBegin;
  err = dMemzero(a,D*Q[0]*P[1]*sizeof(a[0]));dCHK(err);
  switch (imode) {
    case INSERT_VALUES:
      err = dMemzero(g,D*Q[0]*Q[1]*sizeof(g[0]));dCHK(err);
      break;
    case ADD_VALUES:
      break;
    default:
      dERROR(1,"Requested InsertMode %d not supported for this operation.",imode);
  }

  for (l=0; l<Q[0]; l++) {
    for (i=0; i<P[0]; i++) {
      for (j=0; j<P[1]; j++) {
        for (d=0; d<D; d++) {
          a[(l*P[1]+j)*D+d] += A[0][l*P[0]+i] * f[(i*P[1]+j)*D+d];
        }
      }
    }
  }
  for (i=0; i<Q[0]; i++) {
    for (l=0; l<Q[1]; l++) {
      for (j=0; j<P[1]; j++) {
        for (d=0; d<D; d++) {
          g[(i*Q[1]+l)*D+d] += A[1][l*P[1]+j] * a[(i*P[1]+j)*D+d];
        }
      }
    }
  }
  dFunctionReturn(0);
}

/**
* The core computational kernel.  Performs a tensor product operation with the matrices A[0..2].
*
* @param[in] D number of degrees of freedom
* @param[in] P array of length 3, dimensions of input block
* @param[in] Q array of length 3, dimensions of output block
* @param[in] A array of pointers, length 3, transformation matrix in each direction, shape(A[i])=(Q[i],P[i]), i=0,1,2
* @param[in] f input vector with size corresponding to \in
* @param[in,out] g output vector with size corresponding to \out
* @param[in] imode ADD_VALUES or INSERT_VALUES
*
* @return err
*/
static dErr TensorMult_Hex(dInt D,const dInt P[3],const dInt Q[3],dReal *A[3],const dScalar in[],dScalar out[restrict],InsertMode imode)
{
  dScalar amem[Q[0]*P[1]*P[2]*D],bmem[Q[0]*Q[1]*P[2]*D];
  dScalar (*restrict a)[P[1]][P[2]][D] = (dScalar (*restrict)[P[1]][P[2]][D])amem;
  dScalar (*restrict b)[Q[1]][P[2]][D] = (dScalar (*restrict)[Q[1]][P[2]][D])bmem;
  dScalar (*restrict g)[Q[1]][Q[2]][D] = (dScalar (*restrict)[Q[1]][Q[2]][D])out;
  const dScalar (*restrict f)[P[1]][P[2]][D] = (const dScalar (*)[P[1]][P[2]][D])in;
  const dReal (*Ax)[P[0]] = (const dReal (*)[P[0]])A[0];
  const dReal (*Ay)[P[1]] = (const dReal (*)[P[1]])A[1];
  const dReal (*Az)[P[2]] = (const dReal (*)[P[2]])A[2];
  dInt i,j,k,l,d;
  dErr err;

  dFunctionBegin;
  if (D*Q[0]*P[1]*P[2]*sizeof(amem[0]) != Q[0]*sizeof(a[0])
      || D*Q[0]*Q[1]*Q[2]*sizeof(out[0]) != Q[0]*sizeof(g[0])) dERROR(1,"sizeof not working as expected");
  err = dMemzero(a,Q[0]*sizeof(a[0]));dCHK(err);
  err = dMemzero(b,Q[0]*sizeof(b[0]));dCHK(err);
  switch (imode) {
    case INSERT_VALUES:
      err = dMemzero(g,Q[0]*sizeof(g[0]));dCHK(err);
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
            //amem[((l*P[1]+j)*P[2]+k)*D+d] += A[0][l*P[0]+i] * in[((i*P[1]+j)*P[2]+k)*D+d];
            a[l][j][k][d] += Ax[l][i] * f[i][j][k][d];
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
            //bmem[((i*Q[1]+l)*P[2]+k)*D+d] += A[1][l*P[1]+j] * amem[((i*P[1]+j)*P[2]+k)*D+d];
            b[i][l][k][d] += Ay[l][j] * a[i][j][k][d];
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
            //out[((i*Q[1]+j)*Q[2]+l)*D+d] += A[2][l*P[2]+k] * bmem[((i*Q[1]+j)*P[2]+k)*D+d];
            g[i][j][l][d] += Az[l][k] * b[i][j][k][d];
          }
        }
      }
    }
  }
  PetscLogFlops((Q[0]*P[0]*P[1]*P[2] + Q[0]*Q[1]*P[1]*P[2] + Q[0]*Q[1]*Q[2]*P[2])*D*2);
  dFunctionReturn(0);
}


static dErr TensorRuleNoMapping(dInt Q,dInt D,const dScalar u[3][Q][D],dScalar v[Q][D][3],InsertMode imode)
{

  dFunctionBegin;
  switch (imode) {
    case INSERT_VALUES: {
      for (dInt i=0; i<Q; i++) {
        for (dInt d=0; d<D; d++) { /* component */
          for (dInt e=0; e<3; e++) { /* direction */
            v[i][d][e] = u[e][i][d];
          }
        }
      }
    } break;
    case ADD_VALUES: {
      for (dInt i=0; i<Q; i++) {
        for (dInt d=0; d<D; d++) { /* component */
          for (dInt e=0; e<3; e++) { /* direction */
            v[i][d][e] += u[e][i][d];
          }
        }
      }
    } break;
    default: dERROR(1,"InsertMode %d invalid/unimplemented",imode);
  }
  dFunctionReturn(0);
}

/**
* The following functions are defined in this file so that they can be static.  Perhaps they should become virtual and
* be put in ruletopo.c
**/
static dErr dRuleMappingApply_Tensor_Line(dRule_Tensor *rule,const dReal mapdata[],dInt D,const dScalar in[],dScalar out[],InsertMode imode)
{
  const dInt Q = rule->trule[0]->size;
  const dScalar (*u)[Q][D] = (const dScalar (*)[Q][D])in;
  dScalar (*v)[D][3] = (dScalar (*)[D][3])out;
  dErr err;

  dFunctionBegin;
  if (mapdata) {
    dERROR(1,"Not implemented");
  } else {
    err = TensorRuleNoMapping(Q,D,u,v,imode);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr dRuleMappingApply_Tensor_Quad(dRule_Tensor *rule,const dReal mapdata[],dInt D,const dScalar in[],dScalar out[],InsertMode imode)
{
  const dInt Q[2] = {rule->trule[0]->size,rule->trule[1]->size},QQ = Q[0]*Q[1];
  const dScalar (*u)[QQ][D] = (const dScalar (*)[QQ][D])in;
  dScalar (*v)[D][3] = (dScalar (*)[D][3])out;
  dErr err;

  dFunctionBegin;
  if (mapdata) {
    dERROR(1,"Not implemented");
  } else {
    err = TensorRuleNoMapping(QQ,D,u,v,imode);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr dRuleMappingApply_Tensor_Hex(dRule_Tensor *rule,const dReal mapdata[],dInt D,const dScalar in[],dScalar out[],InsertMode imode)
{
  const dInt Q[3] = {rule->trule[0]->size,rule->trule[1]->size,rule->trule[2]->size},QQ = Q[0]*Q[1]*Q[2];
  const dScalar (*u)[QQ][D] = (const dScalar (*)[QQ][D])in;
  dScalar (*restrict v)[D][3] = (dScalar(*)[D][3])out;
  dErr err;

  dFunctionBegin;
  if (mapdata) {
    dERROR(1,"Not implemented");
  } else {
    err = TensorRuleNoMapping(QQ,D,u,v,imode);dCHK(err);
  }
  dFunctionReturn(0);
}
