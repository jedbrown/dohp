#include <dohpgeom.h>
#include "tensorquad.h"

static dErr dRuleView_Tensor_Private(const char*,dInt,TensorRule*,PetscViewer);
static dErr dRuleGetPatches_Tensor_All(dRule,dInt*,dInt*,const dInt **,const dReal**);

#define _F(f) static dErr f(dRule,PetscViewer)
_F(dRuleView_Tensor_Line);
_F(dRuleView_Tensor_Quad);
_F(dRuleView_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dRule,dInt*,dInt*)
_F(dRuleGetSize_Tensor_Line);
_F(dRuleGetSize_Tensor_Quad);
_F(dRuleGetSize_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dRule,dInt*,dInt[],const dReal*[],const dReal*[])
_F(dRuleGetTensorNodeWeight_Tensor_Line);
_F(dRuleGetTensorNodeWeight_Tensor_Quad);
_F(dRuleGetTensorNodeWeight_Tensor_Hex);
#undef _F
#define _F(f) static dErr f(dRule,const dReal[restrict][3],dReal[restrict][3],dReal[restrict][3][3],dReal[restrict])
/* _F(dRuleComputeGeometry_Tensor_Line); */
/* _F(dRuleComputeGeometry_Tensor_Quad); */
_F(dRuleComputeGeometry_Tensor_Hex);
#undef _F

static dErr dRuleView_Tensor_Private(const char *name,dInt n,TensorRule *tr,PetscViewer viewer)
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&ascii);dCHK(err);
  if (ascii) {
    err = PetscViewerASCIIPrintf(viewer,"dRule type %s\n",name);dCHK(err);
    for (dInt i=0; i<n; i++) {
      err = TensorRuleView(tr[i],viewer);
    }
  }
  dFunctionReturn(0);
}

static dErr dRuleView_Tensor_Line(dRule rule,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dRuleView_Tensor_Private("Tensor_Line",1,((dRule_Tensor*)rule)->trule,viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dRuleView_Tensor_Quad(dRule rule,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dRuleView_Tensor_Private("Tensor_Quad",2,((dRule_Tensor*)rule)->trule,viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dRuleView_Tensor_Hex(dRule rule,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dRuleView_Tensor_Private("Tensor_Hex",3,((dRule_Tensor*)rule)->trule,viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dRuleGetSize_Tensor_Line(dRule rule,dInt *dim,dInt *size)
{
  dFunctionBegin;
  if (dim) *dim = 1;
  if (size) *size = ((dRule_Tensor*)rule)->trule[0]->size;
  dFunctionReturn(0);
}

static dErr dRuleGetSize_Tensor_Quad(dRule rule,dInt *dim,dInt *size)
{
  TensorRule *r = ((dRule_Tensor*)rule)->trule;
  dFunctionBegin;
  if (dim) *dim = 2;
  if (size) *size = r[0]->size * r[1]->size;
  dFunctionReturn(0);
}

static dErr dRuleGetSize_Tensor_Hex(dRule rule,dInt *dim,dInt *size)
{
  TensorRule *r = ((dRule_Tensor*)rule)->trule;
  dFunctionBegin;
  if (dim) *dim = 3;
  if (size) *size = r[0]->size * r[1]->size * r[2]->size;
  dFunctionReturn(0);
}

static dErr dRuleGetTensorNodeWeight_Tensor_Line(dRule rule,dInt *dim,dInt nnodes[],const dReal *coord[],const dReal *weight[])
{
  TensorRule *r = ((dRule_Tensor*)rule)->trule;
  dFunctionBegin;
  if (dim) *dim = 1;
  if (nnodes) {
    nnodes[0] = r[0]->size;
    nnodes[1] = 0;
    nnodes[2] = 0;
  }
  if (coord) {
    coord[0] = r[0]->coord;
    coord[1] = NULL;
    coord[2] = NULL;
  }
  if (weight) {
    weight[0] = r[0]->weight;
    weight[1] = NULL;
    weight[2] = NULL;
  }
  dFunctionReturn(0);
}

static dErr dRuleGetTensorNodeWeight_Tensor_Quad(dRule rule,dInt *dim,dInt nnodes[],const dReal *coord[],const dReal *weight[])
{
  TensorRule *r = ((dRule_Tensor*)rule)->trule;
  dFunctionBegin;
  if (dim) *dim = 2;
  if (nnodes) {
    nnodes[0] = r[0]->size;
    nnodes[1] = r[1]->size;
    nnodes[2] = 0;
  }
  if (coord) {
    coord[0] = r[0]->coord;
    coord[1] = r[1]->coord;
    coord[2] = NULL;
  }
  if (weight) {
    weight[0] = r[0]->weight;
    weight[1] = r[1]->weight;
    weight[2] = NULL;
  }
  dFunctionReturn(0);
}

static dErr dRuleGetTensorNodeWeight_Tensor_Hex(dRule rule,dInt *dim,dInt nnodes[],const dReal *coord[],const dReal *weight[])
{
  TensorRule *r = ((dRule_Tensor*)rule)->trule;
  dFunctionBegin;
  if (dim) *dim = 3;
  if (nnodes) {
    nnodes[0] = r[0]->size;
    nnodes[1] = r[1]->size;
    nnodes[2] = r[2]->size;
  }
  if (coord) {
    coord[0] = r[0]->coord;
    coord[1] = r[1]->coord;
    coord[2] = r[2]->coord;
  }
  if (weight) {
    weight[0] = r[0]->weight;
    weight[1] = r[1]->weight;
    weight[2] = r[2]->weight;
  }
  dFunctionReturn(0);
}

static dErr dRuleComputeGeometry_Tensor_Hex(dRule rule,const dReal x[restrict][3],dReal qg[restrict][3],dReal jinv[restrict][3][3],dReal jw[restrict])
{
  const TensorRule *r = ((dRule_Tensor*)rule)->trule;
  const dInt Q[3] = {r[0]->size,r[1]->size,r[2]->size},QQ = Q[0]*Q[1]*Q[2];
  const dReal *qx[3] = {r[0]->coord,r[1]->coord,r[2]->coord};
  const dReal *qw[3] = {r[0]->weight,r[1]->weight,r[2]->weight};
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<Q[0]; i++) {
    const dReal q0 = qx[0][i],q0m = 0.125*(1-q0),q0p = 0.125*(1+q0),qmd = q0m,qpd = q0p;
    for (dInt j=0; j<Q[1]; j++) {
      const dReal q1 = qx[1][j],q1m = 1-q1,q1p = 1+q1;
      const dReal qmm = q0m*q1m;
      const dReal qpm = q0p*q1m;
      const dReal qmp = q0m*q1p;
      const dReal qpp = q0p*q1p;
      const dReal qdm = 0.125*q1m,qdp = 0.125*q1p;
      for (dInt k=0; k<Q[2]; k++) {
        const dInt p = (i*Q[1]+j)*Q[2]+k;                /* Index of quadrature point */
        const dReal q2 = qx[2][k],q2m = 1-q2,q2p = 1+q2; /* location of quadrature point in reference coordinates */
        const dReal qdmm = qdm*q2m,qdmp = qdm*q2p,qdpm = qdp*q2m,qdpp = qdp*q2p;
        const dReal qmdm = qmd*q2m,qmdp = qmd*q2p,qpdm = qpd*q2m,qpdp = qpd*q2p;
        const dReal qmmm = qmm*q2m,qmmp = qmm*q2p,qmpm = qmp*q2m,qmpp = qmp*q2p;
        const dReal qpmm = qpm*q2m,qpmp = qpm*q2p,qppm = qpp*q2m,qppp = qpp*q2p;
        dReal J[3][3],Jdet;
        for (dInt l=0; l<3; l++) {
          qg[p][l] = (+ x[0][l]*qmmm
                      + x[1][l]*qpmm
                      + x[2][l]*qppm
                      + x[3][l]*qmpm
                      + x[4][l]*qmmp
                      + x[5][l]*qpmp
                      + x[6][l]*qppp
                      + x[7][l]*qmpp);
          J[l][0] =  (- x[0][l]*qdmm
                      + x[1][l]*qdmm
                      + x[2][l]*qdpm
                      - x[3][l]*qdpm
                      - x[4][l]*qdmp
                      + x[5][l]*qdmp
                      + x[6][l]*qdpp
                      - x[7][l]*qdpp);
          J[l][1] =  (- x[0][l]*qmdm
                      - x[1][l]*qpdm
                      + x[2][l]*qpdm
                      + x[3][l]*qmdm
                      - x[4][l]*qmdp
                      - x[5][l]*qpdp
                      + x[6][l]*qpdp
                      + x[7][l]*qmdp);
          J[l][2] =  (- x[0][l]*qmm
                      - x[1][l]*qpm
                      - x[2][l]*qpp
                      - x[3][l]*qmp
                      + x[4][l]*qmm
                      + x[5][l]*qpm
                      + x[6][l]*qpp
                      + x[7][l]*qmp);
        }
        err = dGeomInvert3(&J[0][0],&jinv[p][0][0],&Jdet);dCHK(err);
        if (Jdet <= 0.0) dERROR(PETSC_COMM_SELF,1,"Negative Jacobian at %d,%d,%d",i,j,k);
        jw[p] = Jdet * qw[0][i] * qw[1][j] * qw[2][k]; /* Weight the determinant */
      }
    }
  }
  err = PetscLogFlops(Q[0]*(4 + Q[1]*8) + QQ*(18 /* prep */ + 3*4*15 /* qg,J */ + 42 /* invert */ + 3 /* weight */));dCHK(err);
  dFunctionReturn(0);
}

static dErr dRuleGetPatches_Tensor_All(dRule grule,dInt *npatches,dInt *patchsize,const dInt **ind,const dReal **weight)
{
  dRule_Tensor *rule = (dRule_Tensor*)grule;
  dFunctionBegin;
  *npatches  = rule->npatches;
  *patchsize = rule->patchsize;
  *ind       = rule->patchind;
  *weight    = rule->patchweight;
  dFunctionReturn(0);
}

dErr dQuadratureRuleOpsSetUp_Tensor(dQuadrature quad)
{
  static const struct _dRuleOps ruleOpsLine = {
    .view                = dRuleView_Tensor_Line,
    .getSize             = dRuleGetSize_Tensor_Line,
    .getNodeWeight       = NULL, /* dRuleGetNodeWeight_Tensor_Line, */
    .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Line,
    .getPatches          = dRuleGetPatches_Tensor_All,
    .computeGeometry     = NULL, /* Not implemented */
  };
  static const struct _dRuleOps ruleOpsQuad = {
    .view                = dRuleView_Tensor_Quad,
    .getSize             = dRuleGetSize_Tensor_Quad,
    .getNodeWeight       = NULL, /* dRuleGetNodeWeight_Tensor_Quad, */
    .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Quad,
    .getPatches          = dRuleGetPatches_Tensor_All,
    .computeGeometry     = NULL, /* Not implemented */
  };
  static const struct _dRuleOps ruleOpsHex  = {
    .view                = dRuleView_Tensor_Hex,
    .getSize             = dRuleGetSize_Tensor_Hex,
    .getNodeWeight       = NULL, /* dRuleGetNodeWeight_Tensor_Hex, */
    .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Hex,
    .getPatches          = dRuleGetPatches_Tensor_All,
    .computeGeometry     = dRuleComputeGeometry_Tensor_Hex,
  };
  dQuadrature_Tensor *tnsr = quad->data;
  dErr err;

  dFunctionBegin;
  err = dMalloc(sizeof(struct _dRuleOps),&tnsr->ruleOpsLine);dCHK(err);
  err = dMemcpy(tnsr->ruleOpsLine,&ruleOpsLine,sizeof(struct _dRuleOps));
  err = dMalloc(sizeof(struct _dRuleOps),&tnsr->ruleOpsQuad);dCHK(err);
  err = dMemcpy(tnsr->ruleOpsQuad,&ruleOpsQuad,sizeof(struct _dRuleOps));
  err = dMalloc(sizeof(struct _dRuleOps),&tnsr->ruleOpsHex);dCHK(err);
  err = dMemcpy(tnsr->ruleOpsHex,&ruleOpsHex,sizeof(struct _dRuleOps));
  dFunctionReturn(0);
}
