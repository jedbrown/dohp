#include "dohpgeom.h"
#include "tensor.h"

static dErr dRuleView_Tensor_Private(const char*,dInt,TensorRule*,PetscViewer);
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
#define _F(f) static dErr f(dRule*,dReal*,dReal*)
/* _F(dRuleGetNodeWeight_Tensor_Line); */
/* _F(dRuleGetNodeWeight_Tensor_Quad); */
/* _F(dRuleGetNodeWeight_Tensor_Line); */
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

/**
* Set up the rule ops table for each topology.
*
* @param jac
*
* @return
*/
dErr dJacobiRuleOpsSetUp_Tensor(dJacobi jac)
{
  static const struct _dRuleOps ruleOpsLine = {
    .view                = dRuleView_Tensor_Line,
    .getSize             = dRuleGetSize_Tensor_Line,
    .getNodeWeight       = NULL, /* dRuleGetNodeWeight_Tensor_Line, */
    .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Line,
    .computeGeometry     = NULL, /* Not implemented */
  };
  static const struct _dRuleOps ruleOpsQuad = {
    .view                = dRuleView_Tensor_Quad,
    .getSize             = dRuleGetSize_Tensor_Quad,
    .getNodeWeight       = NULL, /* dRuleGetNodeWeight_Tensor_Quad, */
    .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Quad,
    .computeGeometry     = NULL, /* Not implemented */
  };
  static const struct _dRuleOps ruleOpsHex  = {
    .view                = dRuleView_Tensor_Hex,
    .getSize             = dRuleGetSize_Tensor_Hex,
    .getNodeWeight       = NULL, /* dRuleGetNodeWeight_Tensor_Hex, */
    .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Hex,
    .computeGeometry     = dRuleComputeGeometry_Tensor_Hex,
  };
  Tensor this = (Tensor)jac->impl;
  dErr err;

  dFunctionBegin;
  if (!this->ruleOpsLine) {
    err = dMalloc(sizeof(struct _dRuleOps),&this->ruleOpsLine);dCHK(err);
    err = dMemcpy(this->ruleOpsLine,&ruleOpsLine,sizeof(struct _dRuleOps));
  }
  if (!this->ruleOpsQuad) {
    err = dMalloc(sizeof(struct _dRuleOps),&this->ruleOpsQuad);dCHK(err);
    err = dMemcpy(this->ruleOpsQuad,&ruleOpsQuad,sizeof(struct _dRuleOps));
  }
  if (!this->ruleOpsHex) {
    err = dMalloc(sizeof(struct _dRuleOps),&this->ruleOpsHex);dCHK(err);
    err = dMemcpy(this->ruleOpsHex,&ruleOpsHex,sizeof(struct _dRuleOps));
  }
  dFunctionReturn(0);
}

dErr dJacobiRuleOpsDestroy_Tensor(dJacobi jac)
{
  Tensor this = (Tensor)jac->impl;
  dErr err;

  dFunctionBegin;
  if (this->ruleOpsLine) { err = dFree(this->ruleOpsLine);dCHK(err); }
  if (this->ruleOpsQuad) { err = dFree(this->ruleOpsQuad);dCHK(err); }
  if (this->ruleOpsHex)  { err = dFree(this->ruleOpsHex);dCHK(err); }
  dFunctionReturn(0);
}

static dErr dRuleView_Tensor_Private(const char *name,dInt n,TensorRule *tr,PetscViewer viewer)
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
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

static dErr dRuleComputeGeometry_Tensor_Hex(dRule rule,const dReal x[restrict][3],dReal qg[restrict][3],dReal jinv[restrict][3][3],dReal jdet[restrict])
{
  const TensorRule *r = ((dRule_Tensor*)rule)->trule;
  const dInt Q[3] = {r[0]->size,r[1]->size,r[2]->size},QQ = Q[0]*Q[1]*Q[2];
  const dReal *qx[3] = {r[0]->coord,r[1]->coord,r[2]->coord};
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<Q[0]; i++) {
    for (dInt j=0; j<Q[1]; j++) {
      for (dInt k=0; k<Q[2]; k++) {
        const dInt p = (i*Q[1]+j)*Q[2]+k;                /* Index of quadrature point */
        const dReal q[3] = {qx[0][i],qx[1][j],qx[2][k]}; /* location of quadrature point in reference coordinates */
        dReal J[3][3],f[6][3];
        err = dGeomConvexComb_2_4(q[0],q[2],x,dMeshConnectHexQuad[0],f[0]);dCHK(err);
        err = dGeomConvexComb_2_4(q[1],q[2],x,dMeshConnectHexQuad[1],f[1]);dCHK(err);
        err = dGeomConvexComb_2_4(-q[0],q[2],x,dMeshConnectHexQuad[2],f[2]);dCHK(err);
        err = dGeomConvexComb_2_4(-q[1],q[2],x,dMeshConnectHexQuad[3],f[3]);dCHK(err);
        err = dGeomConvexComb_2_4(q[0],q[1],x,dMeshConnectHexQuad[4],f[4]);dCHK(err);
        err = dGeomConvexComb_2_4(q[0],q[1],x,dMeshConnectHexQuad[5],f[5]);dCHK(err);
        for (dInt l; l<3; l++) {
          J[l][0] = 0.5 * (f[1][l] - f[3][l]);
          J[l][1] = 0.5 * (f[2][l] - f[0][l]);
          J[l][2] = 0.5 * (f[5][l] - f[4][l]);
          qg[p][l] = 0.5 * ((1.0-q[1])*f[0][l] + (1.0+q[1])*f[2][l]);
        }
        err = dGeomInvert3(&J[0][0],&jinv[p][0][0],&jdet[p]);dCHK(err);
      }
    }
  }
  /* Assume the optimizer has eliminated common subexpressions, this estimate is about right, some of the convex
  * combinations could be hoisted out of the inner loop */
  PetscLogFlops(QQ*(6+6*30+6+42));
  dFunctionReturn(0);
}
