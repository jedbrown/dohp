#include "tensor.h"

static dErr dRuleView_Tensor_Private(const char*,dInt,TensorRule*,PetscViewer);
static dErr dRuleView_Tensor_Line(dRule,PetscViewer);
static dErr dRuleView_Tensor_Quad(dRule,PetscViewer);
static dErr dRuleView_Tensor_Hex(dRule,PetscViewer);
static dErr dRuleGetSize_Tensor_Line(dRule,dInt*,dInt*);
static dErr dRuleGetSize_Tensor_Quad(dRule,dInt*,dInt*);
static dErr dRuleGetSize_Tensor_Hex(dRule,dInt*,dInt*);
#if 0
static dErr dRuleGetNodeWeight_Tensor_Line(dRule*,dReal*,dReal*);
static dErr dRuleGetNodeWeight_Tensor_Quad(dRule*,dReal*,dReal*);
static dErr dRuleGetNodeWeight_Tensor_Hex(dRule*,dReal*,dReal*);
#endif
static dErr dRuleGetTensorNodeWeight_Tensor_Line(dRule,dInt*,dInt[],const dReal**,const dReal**);
static dErr dRuleGetTensorNodeWeight_Tensor_Quad(dRule,dInt*,dInt[],const dReal**,const dReal**);
static dErr dRuleGetTensorNodeWeight_Tensor_Hex(dRule,dInt*,dInt[],const dReal**,const dReal**);

/**
* Set up the rule ops table for each topology.
*
* @param jac
*
* @return
*/
dErr dJacobiRuleOpsSetUp_Tensor(dJacobi jac)
{
  static const struct _dRuleOps ruleOpsLine = { .view = dRuleView_Tensor_Line,
                                                .getSize = dRuleGetSize_Tensor_Line,
                                                .getNodeWeight = NULL, /* dRuleGetNodeWeight_Tensor_Line, */
                                                .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Line };
  static const struct _dRuleOps ruleOpsQuad = { .view = dRuleView_Tensor_Quad,
                                                .getSize = dRuleGetSize_Tensor_Quad,
                                                .getNodeWeight = NULL, /* dRuleGetNodeWeight_Tensor_Quad, */
                                                .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Quad };
  static const struct _dRuleOps ruleOpsHex  = { .view = dRuleView_Tensor_Hex,
                                                .getSize = dRuleGetSize_Tensor_Hex,
                                                .getNodeWeight = NULL, /* dRuleGetNodeWeight_Tensor_Hex, */
                                                .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Hex };
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
