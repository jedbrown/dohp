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
  static const struct v_dRuleOps ruleOpsLine = { .view = dRuleView_Tensor_Line,
                                                 .getSize = dRuleGetSize_Tensor_Line,
                                                 .getNodeWeight = NULL, /* dRuleGetNodeWeight_Tensor_Line, */
                                                 .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Line };
  static const struct v_dRuleOps ruleOpsQuad = { .view = dRuleView_Tensor_Quad,
                                                 .getSize = dRuleGetSize_Tensor_Quad,
                                                 .getNodeWeight = NULL, /* dRuleGetNodeWeight_Tensor_Quad, */
                                                 .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Quad };
  static const struct v_dRuleOps ruleOpsHex  = { .view = dRuleView_Tensor_Hex,
                                                 .getSize = dRuleGetSize_Tensor_Hex,
                                                 .getNodeWeight = NULL, /* dRuleGetNodeWeight_Tensor_Hex, */
                                                 .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Hex };
  Tensor this = (Tensor)jac->impl;
  dErr err;

  dFunctionBegin;
  if (!this->ruleOpsLine) {
    err = dMalloc(sizeof(struct v_dRuleOps),&this->ruleOpsLine);dCHK(err);
    err = dMemcpy(this->ruleOpsLine,&ruleOpsLine,sizeof(struct v_dRuleOps));
  }
  if (!this->ruleOpsQuad) {
    err = dMalloc(sizeof(struct v_dRuleOps),&this->ruleOpsQuad);dCHK(err);
    err = dMemcpy(this->ruleOpsQuad,&ruleOpsQuad,sizeof(struct v_dRuleOps));
  }
  if (!this->ruleOpsHex) {
    err = dMalloc(sizeof(struct v_dRuleOps),&this->ruleOpsHex);dCHK(err);
    err = dMemcpy(this->ruleOpsHex,&ruleOpsHex,sizeof(struct v_dRuleOps));
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
  err = dRuleView_Tensor_Private("Tensor_Line",1,((struct s_dRule_Tensor_Line*)rule)->trule,viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dRuleView_Tensor_Quad(dRule rule,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dRuleView_Tensor_Private("Tensor_Quad",2,((struct s_dRule_Tensor_Quad*)rule)->trule,viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dRuleView_Tensor_Hex(dRule rule,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dRuleView_Tensor_Private("Tensor_Hex",3,((struct s_dRule_Tensor_Hex*)rule)->trule,viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dRuleGetSize_Tensor_Line(dRule rule,dInt *dim,dInt *size)
{
  dFunctionBegin;
  if (dim) *dim = 1;
  if (size) *size = ((struct s_dRule_Tensor_Line*)rule)->trule[0]->size;
  dFunctionReturn(0);
}

static dErr dRuleGetSize_Tensor_Quad(dRule rule,dInt *dim,dInt *size)
{
  TensorRule *r = ((struct s_dRule_Tensor_Quad*)rule)->trule;
  dFunctionBegin;
  if (dim) *dim = 2;
  if (size) *size = r[0]->size * r[1]->size;
  dFunctionReturn(0);
}

static dErr dRuleGetSize_Tensor_Hex(dRule rule,dInt *dim,dInt *size)
{
  TensorRule *r = ((struct s_dRule_Tensor_Hex*)rule)->trule;
  dFunctionBegin;
  if (dim) *dim = 3;
  if (size) *size = r[0]->size * r[1]->size * r[2]->size;
  dFunctionReturn(0);
}

static dErr dRuleGetTensorNodeWeight_Tensor_Line(dRule rule,dInt *dim,dInt nnodes[],const dReal *coord[],const dReal *weight[])
{
  TensorRule *r = ((struct s_dRule_Tensor_Line*)rule)->trule;
  dFunctionBegin;
  if (dim) *dim = 1;
  if (nnodes) nnodes[0] = r[0]->size;
  if (coord) coord[0] = r[0]->coord;
  if (weight) weight[0] = r[0]->weight;
  dFunctionReturn(0);
}

static dErr dRuleGetTensorNodeWeight_Tensor_Quad(dRule rule,dInt *dim,dInt nnodes[],const dReal *coord[],const dReal *weight[])
{
  TensorRule *r = ((struct s_dRule_Tensor_Quad*)rule)->trule;
  dFunctionBegin;
  if (dim) *dim = 2;
  if (nnodes) {
    nnodes[0] = r[0]->size;
    nnodes[1] = r[1]->size;
  }
  if (coord) {
    coord[0] = r[0]->coord;
    coord[1] = r[1]->coord;
  }
  if (weight) {
    weight[0] = r[0]->weight;
    weight[1] = r[1]->weight;
  }
  dFunctionReturn(0);
}

static dErr dRuleGetTensorNodeWeight_Tensor_Hex(dRule rule,dInt *dim,dInt nnodes[],const dReal *coord[],const dReal *weight[])
{
  TensorRule *r = ((struct s_dRule_Tensor_Hex*)rule)->trule;
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
