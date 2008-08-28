#include "tensor.h"

static dErr dRuleView_Tensor_Line(dRule*,PetscViewer);
static dErr dRuleView_Tensor_Quad(dRule*,PetscViewer);
static dErr dRuleView_Tensor_Hex(dRule*,PetscViewer);
static dErr dRuleGetSize_Tensor_Line(dRule*,dInt*,dInt*);
static dErr dRuleGetSize_Tensor_Quad(dRule*,dInt*,dInt*);
static dErr dRuleGetSize_Tensor_Hex(dRule*,dInt*,dInt*);
static dErr dRuleGetNodeWeight_Tensor_Line(dRule*,dReal*,dReal*);
static dErr dRuleGetNodeWeight_Tensor_Quad(dRule*,dReal*,dReal*);
static dErr dRuleGetNodeWeight_Tensor_Hex(dRule*,dReal*,dReal*);
static dErr dRuleGetTensorNodeWeight_Tensor_Line(dRule*,dInt*,dInt[],const dReal**,const dReal**);
static dErr dRuleGetTensorNodeWeight_Tensor_Quad(dRule*,dInt*,dInt[],const dReal**,const dReal**);
static dErr dRuleGetTensorNodeWeight_Tensor_Hex(dRule*,dInt*,dInt[],const dReal**,const dReal**);

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
                                                 .getNodeWeight = dRuleGetNodeWeight_Tensor_Line,
                                                 .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Line };
  static const struct v_dRuleOps ruleOpsQuad = { .view = dRuleView_Tensor_Quad,
                                                 .getSize = dRuleGetSize_Tensor_Quad,
                                                 .getNodeWeight = dRuleGetNodeWeight_Tensor_Quad,
                                                 .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Quad };
  static const struct v_dRuleOps ruleOpsHex  = { .view = dRuleView_Tensor_Hex,
                                                 .getSize = dRuleGetSize_Tensor_Hex,
                                                 .getNodeWeight = dRuleGetNodeWeight_Tensor_Hex,
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

