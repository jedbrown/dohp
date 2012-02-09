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

static dErr dRuleGetPatches_Tensor_All(dRule grule,dInt *npatches,dInt *patchsize,const dInt **ind,const dReal **weight)
{
  dRule_Tensor *rule = (dRule_Tensor*)grule;
  dFunctionBegin;
  if (npatches)  *npatches  = rule->npatches;
  if (patchsize) *patchsize = rule->patchsize;
  if (ind)       *ind       = rule->patchind;
  if (weight)    *weight    = rule->patchweight;
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
  };
  static const struct _dRuleOps ruleOpsQuad = {
    .view                = dRuleView_Tensor_Quad,
    .getSize             = dRuleGetSize_Tensor_Quad,
    .getNodeWeight       = NULL, /* dRuleGetNodeWeight_Tensor_Quad, */
    .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Quad,
    .getPatches          = dRuleGetPatches_Tensor_All,
  };
  static const struct _dRuleOps ruleOpsHex  = {
    .view                = dRuleView_Tensor_Hex,
    .getSize             = dRuleGetSize_Tensor_Hex,
    .getNodeWeight       = NULL, /* dRuleGetNodeWeight_Tensor_Hex, */
    .getTensorNodeWeight = dRuleGetTensorNodeWeight_Tensor_Hex,
    .getPatches          = dRuleGetPatches_Tensor_All,
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
