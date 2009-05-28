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
#define _F(f) static dErr f(dEFS,dInt*,dInt*,dReal**,dReal**,const dReal**,const dReal**) /* dEFSGetTensorNodes */
_F(dEFSGetTensorNodes_Tensor_Line);
_F(dEFSGetTensorNodes_Tensor_Quad);
_F(dEFSGetTensorNodes_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS,const dReal[restrict],dInt,const dScalar[restrict],dScalar[restrict],dApplyMode,InsertMode) /* dEFSApply */
_F(dEFSApply_Tensor_Line);
_F(dEFSApply_Tensor_Quad);
_F(dEFSApply_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dEFS,const dReal(*)[3],dInt*,dInt[],dReal(*)[3]) /* dEFSGetGlobalCoordinates */
//_F(dEFSGetGlobalCoordinates_Tensor_Line);
//_F(dEFSGetGlobalCoordinates_Tensor_Quad);
_F(dEFSGetGlobalCoordinates_Tensor_Hex);
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

#define _F(f) static dErr f(dRule_Tensor*,const dReal[],dInt,const dScalar[],dScalar[],InsertMode)
_F(dRuleMappingApply_Tensor_Line);
_F(dRuleMappingApply_Tensor_Quad);
//_F(dRuleMappingApply_Tensor_Hex);
#undef _F

static dErr TensorRuleMapping(dInt Q,const dReal jinv_flat[restrict],dInt D,const dScalar in[restrict],dScalar out[restrict],InsertMode imode);
static dErr TensorRuleMappingTranspose(dInt Q,const dReal jinv_flat[restrict],dInt D,const dScalar in[restrict],dScalar out[restrict],InsertMode imode);

#define _F(f) static dErr f(dInt D,const dInt P[],const dInt Q[],const dReal *A[],const dScalar f[],dScalar g[restrict],InsertMode imode)
_F(TensorMult_Line);
_F(TensorMult_Quad);
#undef _F

static inline dErr TensorMult_Hex(const TensorBasis b[],dInt D,const dInt P[],const dInt Q[],const dReal *A[],const dScalar f[],dScalar g[restrict],InsertMode imode)
{
#if defined(dUSE_DEBUG)
  if (D > 3) dERROR(1,"D > 3 not supported");
  if (!b[2]->multhex[D-1]) dERROR(1,"No multhex member, EFS was not set up correctly");
#endif
  return b[2]->multhex[D-1](D,P,Q,A,f,g,imode);
}


/**
* Set up the EFS ops table for each topology.  This is the only exported function in this file.
*
* @param jac
*
* @return
*/
dErr dJacobiEFSOpsSetUp_Tensor(dJacobi jac)
{
  static const struct _dEFSOps efsOpsLine = { .view                 = dEFSView_Tensor_Line,
                                              .getSizes             = dEFSGetSizes_Tensor_Line,
                                              .getTensorNodes       = dEFSGetTensorNodes_Tensor_Line,
                                              .apply                = dEFSApply_Tensor_Line,
                                              .scatterInt           = 0,
                                              .scatterFacet         = 0 };
  static const struct _dEFSOps efsOpsQuad = { .view                 = dEFSView_Tensor_Quad,
                                              .getSizes             = dEFSGetSizes_Tensor_Quad,
                                              .getTensorNodes       = dEFSGetTensorNodes_Tensor_Quad,
                                              .apply                = dEFSApply_Tensor_Quad,
                                              .scatterInt           = 0,
                                              .scatterFacet         = 0 };
  static const struct _dEFSOps efsOpsHex = { .view                 = dEFSView_Tensor_Hex,
                                             .getSizes             = dEFSGetSizes_Tensor_Hex,
                                             .getTensorNodes       = dEFSGetTensorNodes_Tensor_Hex,
                                             .apply                = dEFSApply_Tensor_Hex,
                                             .getGlobalCoordinates = dEFSGetGlobalCoordinates_Tensor_Hex,
                                             .scatterInt           = 0,
                                             .scatterFacet         = 0 };
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

static dErr dEFSGetTensorNodes_Tensor_Line(dEFS efs,dInt *dim,dInt tsize[],dReal *x[],dReal *weight[],const dReal *mscale[],const dReal *lscale[])
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
  if (weight) {
    weight[0] = b[0]->weight;
    weight[1] = NULL;
    weight[2] = NULL;
  }
  if (mscale) {
    mscale[0] = b[0]->mscale;
    mscale[1] = NULL;
    mscale[2] = NULL;
  }
  if (lscale) {
    lscale[0] = b[0]->lscale;
    lscale[1] = NULL;
    lscale[2] = NULL;
  }
  dFunctionReturn(0);
}

static dErr dEFSGetTensorNodes_Tensor_Quad(dEFS efs,dInt *dim,dInt tsize[],dReal *x[],dReal *weight[],const dReal *mscale[],const dReal *lscale[])
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
  if (weight) {
    weight[0] = b[0]->weight;
    weight[1] = b[1]->weight;
    weight[2] = NULL;
  }
  if (mscale) {
    mscale[0] = b[0]->mscale;
    mscale[1] = b[1]->mscale;
    mscale[2] = NULL;
  }
  if (lscale) {
    lscale[0] = b[0]->lscale;
    lscale[1] = b[1]->lscale;
    lscale[2] = NULL;
  }
  dFunctionReturn(0);
}

static dErr dEFSGetTensorNodes_Tensor_Hex(dEFS efs,dInt *dim,dInt tsize[],dReal *x[],dReal *weight[],const dReal *mscale[],const dReal *lscale[])
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
  if (weight) {
    weight[0] = b[0]->weight;
    weight[1] = b[1]->weight;
    weight[2] = b[2]->weight;
  }
  if (mscale) {
    mscale[0] = b[0]->mscale;
    mscale[1] = b[1]->mscale;
    mscale[2] = b[2]->lscale;
  }
  if (lscale) {
    lscale[0] = b[0]->lscale;
    lscale[1] = b[1]->lscale;
    lscale[2] = b[2]->lscale;
  }
  dFunctionReturn(0);
}

static dErr dEFSApply_Tensor_Line(dEFS efs,const dReal mapdata[],dInt D,const dScalar in[],dScalar out[],dApplyMode amode,InsertMode imode)
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;
  const dReal *A[1];
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
  const dReal *A[2];
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

static dErr dEFSApply_Tensor_Hex(dEFS efs,const dReal jinv[restrict],dInt D,const dScalar in[],dScalar out[],dApplyMode amode,InsertMode imode)
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;
  const dReal *A[3];
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
      err = TensorMult_Hex(b,D,P,Q,A,in,out,imode);dCHK(err);
      break;
    case dAPPLY_INTERP_TRANSPOSE:
      A[0] = b[0]->interpTranspose;
      A[1] = b[1]->interpTranspose;
      A[2] = b[2]->interpTranspose;
      err = TensorMult_Hex(b,D,Q,P,A,in,out,imode);dCHK(err);
      break;
    case dAPPLY_GRAD: {
      dScalar df[3][Q[0]*Q[1]*Q[2]*D];
      A[0] = b[0]->deriv; A[1] = b[1]->interp; A[2] = b[2]->interp;
      err = TensorMult_Hex(b,D,P,Q,A,in,df[0],INSERT_VALUES);dCHK(err);
      A[0] = b[0]->interp; A[1] = b[1]->deriv; A[2] = b[2]->interp;
      err = TensorMult_Hex(b,D,P,Q,A,in,df[1],INSERT_VALUES);dCHK(err);
      A[0] = b[0]->interp; A[1] = b[1]->interp; A[2] = b[2]->deriv;
      err = TensorMult_Hex(b,D,P,Q,A,in,df[2],INSERT_VALUES);dCHK(err);
      //err = dRuleMappingApply_Tensor_Hex((dRule_Tensor*)efs->rule,jinv,D,&df[0][0],out,imode);dCHK(err);
      err = TensorRuleMapping(Q[0]*Q[1]*Q[2],jinv,D,&df[0][0],out,imode);dCHK(err);
    } break;
    case dAPPLY_GRAD_TRANSPOSE: {
      dScalar df[3][Q[0]*Q[1]*Q[2]*D];
      err = TensorRuleMappingTranspose(Q[0]*Q[1]*Q[2],jinv,D,in,&df[0][0],INSERT_VALUES);dCHK(err);
      switch (imode) {
        case INSERT_VALUES: err = dMemzero(out,P[0]*P[1]*P[2]*D*sizeof(out[0]));dCHK(err); break;
        case ADD_VALUES: break;
        default: dERROR(1,"InsertMode %d invalid or unimplemented",imode);
      }
      A[0] = b[0]->derivTranspose; A[1] = b[1]->interpTranspose; A[2] = b[2]->interpTranspose;
      err = TensorMult_Hex(b,D,Q,P,A,df[0],out,ADD_VALUES);dCHK(err);
      A[0] = b[0]->interpTranspose; A[1] = b[1]->derivTranspose; A[2] = b[2]->interpTranspose;
      err = TensorMult_Hex(b,D,Q,P,A,df[1],out,ADD_VALUES);dCHK(err);
      A[0] = b[0]->interpTranspose; A[1] = b[1]->interpTranspose; A[2] = b[2]->derivTranspose;
      err = TensorMult_Hex(b,D,Q,P,A,df[2],out,ADD_VALUES);dCHK(err);
    } break;
    default:
      dERROR(1,"invalid/unimplemented dApplyMode %d specified",amode);
  }
  dFunctionReturn(0);
}

static dErr dEFSGetGlobalCoordinates_Tensor_Hex(dEFS efs,const dReal (*x)[3],dInt *dim,dInt P[static 3],dReal (*qx)[3])
{
  TensorBasis *b = ((dEFS_Tensor*)efs)->basis;

  dFunctionBegin;
  *dim = 3; P[0] = b[0]->P; P[1] = b[1]->P; P[2] = b[2]->P;
  for (dInt i=0; i<P[0]; i++) {
    const dReal q0 = b[0]->node[i],q0m = 0.125*(1-q0),q0p = 0.125*(1+q0);
    for (dInt j=0; j<P[1]; j++) {
      const dReal q1 = b[1]->node[j],q1m = 1-q1,q1p = 1+q1;
      const dReal qmm = q0m*q1m;
      const dReal qpm = q0p*q1m;
      const dReal qmp = q0m*q1p;
      const dReal qpp = q0p*q1p;
      for (dInt k=0; k<P[2]; k++) {
        const dInt p = (i*P[1]+j)*P[2]+k;                /* Index of node */
        const dReal q2 = b[2]->node[k],q2m = 1-q2,q2p = 1+q2; /* location of quadrature point in reference coordinates */
        const dReal qmmm = qmm*q2m,qmmp = qmm*q2p,qmpm = qmp*q2m,qmpp = qmp*q2p;
        const dReal qpmm = qpm*q2m,qpmp = qpm*q2p,qppm = qpp*q2m,qppp = qpp*q2p;
        for (dInt l=0; l<3; l++) {
          qx[p][l] = (+ x[0][l]*qmmm + x[1][l]*qpmm + x[2][l]*qppm + x[3][l]*qmpm + x[4][l]*qmmp + x[5][l]*qpmp + x[6][l]*qppp + x[7][l]*qmpp);
        }
      }
    }
  }
  dFunctionReturn(0);
}

static dErr TensorMult_Line(dInt D,const dInt P[1],const dInt Q[1],const dReal *A[1],const dScalar f[],dScalar g[restrict],InsertMode imode)
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

static dErr TensorMult_Quad(dInt D,const dInt P[2],const dInt Q[2],const dReal *A[2],const dScalar f[],dScalar g[restrict],InsertMode imode)
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

static dErr TensorRuleNoMapping(dInt Q,dInt D,const dScalar in[restrict],dScalar out[restrict],InsertMode imode)
{
  const dScalar (*restrict u)[Q][D] = (const dScalar(*)[Q][D])in;
  dScalar (*restrict v)[D][3] = (dScalar(*)[D][3])out;

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

/** Push gradients wrt reference space forward to gradients wrt physical space
*
* @param Q number of quadrature points
* @param jinv_flat inverse jacobian at quadrature points, logically of type 'const dReal[Q][3][3]'
* @param D number of dofs per node
* @param in gradients in reference space, logically 'const dReal[3][Q][D]'
* @param out gradients in physical space, logically 'const dReal[Q][D][3]'
* @param imode INSERT_VALUES or ADD_VALUES
*/
static dErr TensorRuleMapping(dInt Q,const dReal jinv_flat[restrict],dInt D,const dScalar in[restrict],dScalar out[restrict],InsertMode imode)
{
  const dReal (*restrict jinv)[3][3] = (const dReal(*)[3][3])jinv_flat;
  const dScalar (*restrict u)[Q][D] = (const dScalar(*)[Q][D])in;
  dScalar (*restrict v)[D][3] = (dScalar(*)[D][3])out;
  dErr err;

  dFunctionBegin;
  if (!jinv) {
    err = TensorRuleNoMapping(Q,D,in,out,imode);dCHK(err);
    dFunctionReturn(0);
  }
  switch (imode) {
    case INSERT_VALUES: {
      for (dInt i=0; i<Q; i++) {
        for (dInt d=0; d<D; d++) {
          for (dInt e=0; e<3; e++) {
            v[i][d][e] = u[0][i][d] * jinv[i][0][e] + u[1][i][d] * jinv[i][1][e] + u[2][i][d] * jinv[i][2][e];
          }
        }
      }
    } break;
    case ADD_VALUES: {
      for (dInt i=0; i<Q; i++) {
        for (dInt d=0; d<D; d++) {
          for (dInt e=0; e<3; e++) {
            v[i][d][e] += u[0][i][d] * jinv[i][0][e] + u[1][i][d] * jinv[i][1][e] + u[2][i][d] * jinv[i][2][e];
          }
        }
      }
    } break;
    default: dERROR(1,"InsertMode %d invalid or not implemented",imode);
  }
  dFunctionReturn(0);
}

static dErr TensorRuleMappingTranspose(dInt Q,const dReal jinv_flat[restrict],dInt D,const dScalar in[restrict],dScalar out[restrict],InsertMode imode)
{
  const dReal (*restrict jinv)[3][3] = (const dReal(*)[3][3])jinv_flat;
  const dScalar (*restrict u)[D][3] = (const dScalar(*)[D][3])in;
  dScalar (*restrict v)[Q][D] = (dScalar(*)[Q][D])out;

  dFunctionBegin;
  if (!jinv_flat) dERROR(1,"No Jinv, need a mapping");
  switch (imode) {
    case INSERT_VALUES: {
      for (dInt i=0; i<Q; i++) {
        for (dInt d=0; d<D; d++) {
          for (dInt e=0; e<3; e++) {
            v[e][i][d] = u[i][d][0] * jinv[i][e][0] + u[i][d][1] * jinv[i][e][1] + u[i][d][2] * jinv[i][e][2];
          }
        }
      }
    } break;
    case ADD_VALUES: {
      for (dInt i=0; i<Q; i++) {
        for (dInt d=0; d<D; d++) {
          for (dInt e=0; e<3; e++) {
            v[e][i][d] += u[i][d][0] * jinv[i][e][0] + u[i][d][1] * jinv[i][e][1] + u[i][d][2] * jinv[i][e][2];
          }
        }
      }
    } break;
    default: dERROR(1,"InsertMode %d invalid or not implemented",imode);
  }
  dFunctionReturn(0);
}


/**
* The following functions are defined in this file so that they can be static.  Perhaps they should become virtual and
* be put in ruletopo.c
**/
static dErr dRuleMappingApply_Tensor_Line(dRule_Tensor *rule,const dReal jinv[],dInt D,const dScalar in[],dScalar out[],InsertMode imode)
{
  const dInt Q = rule->trule[0]->size;
  dErr err;

  dFunctionBegin;
  if (jinv) {
    dERROR(1,"Not implemented");
  } else {
    err = TensorRuleNoMapping(Q,D,in,out,imode);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr dRuleMappingApply_Tensor_Quad(dRule_Tensor *rule,const dReal jinv[],dInt D,const dScalar in[],dScalar out[],InsertMode imode)
{
  const dInt Q = rule->trule[0]->size * rule->trule[1]->size;
  dErr err;

  dFunctionBegin;
  if (jinv) {
    dERROR(1,"Not implemented");
  } else {
    err = TensorRuleNoMapping(Q,D,in,out,imode);dCHK(err);
  }
  dFunctionReturn(0);
}
