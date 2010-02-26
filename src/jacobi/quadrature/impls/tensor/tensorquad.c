#include <dohp.h>
#include "tensorquad.h"
#include "../../../impls/tensor/polylib.h"

dErr TensorRuleView(TensorRule rule,PetscViewer viewer)
{
  dTruth ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  err = PetscViewerASCIIPrintf(viewer,"TensorRule with %d nodes.\n",rule->size);dCHK(err);
  err = dRealTableView(1,rule->size,rule->coord,"q",viewer);dCHK(err);
  err = dRealTableView(1,rule->size,rule->weight,"w",viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr TensorRuleDestroy(TensorRule rule)
{
  dErr err;

  dFunctionBegin;
  if (!rule) dFunctionReturn(0);
  err = dFree2(rule->weight,rule->coord);dCHK(err);
  err = dFree(rule);dCHK(err);
  dFunctionReturn(0);
}

static dErr TensorGetRule(dQuadrature_Tensor *tnsr,dInt size,TensorRule *rule)
{
  dErr err;
  int key,new;
  khiter_t kiter;

  dFunctionBegin;
  *rule = 0;
  if (!(0 < size && size < 50)) dERROR(1,"rule size out of bounds.");
  if (tnsr->family != GAUSS) dERROR(1,"GaussFamily %d not supported",tnsr->family);
  key = ((int)tnsr->family << 16) | size;
  kiter = kh_put_rule(tnsr->rules,key,&new);
  if (new) {
    TensorRule r;

    err = dNew(struct s_TensorRule,&r);dCHK(err);
    err = dMallocA2(size,&r->weight,size,&r->coord);dCHK(err);
    r->size = size;
    if (size == 1) {              /* Polylib function fails for this size. */
      r->weight[0] = 2.0;
      r->coord[0] = 0.0;
    } else {
      zwgj(r->coord,r->weight,size,tnsr->alpha,tnsr->beta); /* polylib function */
    }
    kh_val(tnsr->rules,kiter) = r;
  }
  *rule = kh_val(tnsr->rules,kiter);
  dFunctionReturn(0);
}

static dErr dQuadratureGetRule_Tensor(dQuadrature quad,dInt n,const dEntTopology topo[],const dInt rsize[],dRule firstrule)
{
  dQuadrature_Tensor *tnsr = quad->data;
  dRule_Tensor *rule = (dRule_Tensor*)firstrule;
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<n; i++) {
    switch (topo[i]) {
      case dTOPO_LINE:
        rule[i].ops = tnsr->ruleOpsLine;
        err = TensorGetRule(tnsr,rsize[3*i+0],&rule[i].trule[0]);dCHK(err);
        if (rsize[3*i+1] != 1 || rsize[3*i+2] != 1)
          dERROR(1,"Invalid rule size for Line");
        rule[i].trule[1] = rule[i].trule[2] = NULL;
        break;
      case dTOPO_QUAD:
        rule[i].ops = tnsr->ruleOpsQuad;
        err = TensorGetRule(tnsr,rsize[3*i+0],&rule[i].trule[0]);dCHK(err);
        err = TensorGetRule(tnsr,rsize[3*i+1],&rule[i].trule[1]);dCHK(err);
        if (rsize[3*i+2] != 1) dERROR(1,"Invalid rule size for Line");
        rule[i].trule[2] = NULL;
        break;
      case dTOPO_HEX:
        rule[i].ops = tnsr->ruleOpsHex;
        err = TensorGetRule(tnsr,rsize[3*i+0],&rule[i].trule[0]);dCHK(err);
        err = TensorGetRule(tnsr,rsize[3*i+1],&rule[i].trule[1]);dCHK(err);
        err = TensorGetRule(tnsr,rsize[3*i+2],&rule[i].trule[2]);dCHK(err);
        break;
      default: dERROR(1,"no rule available for given topology");
    }
  }
  dFunctionReturn(0);
}

static dErr dQuadratureDestroy_Tensor(dQuadrature quad)
{
  dQuadrature_Tensor *tnsr = (dQuadrature_Tensor*)quad->data;
  dErr err;
  khiter_t k;

  dFunctionBegin;
  for (k=kh_begin(tnsr->rules); k!= kh_end(tnsr->rules); k++) {
    if (!kh_exist(tnsr->rules,k)) continue;
    err = TensorRuleDestroy(kh_val(tnsr->rules,k));dCHK(err);
  }
  kh_destroy_rule(tnsr->rules);
  if (tnsr->ruleOpsLine) { err = dFree(tnsr->ruleOpsLine);dCHK(err); }
  if (tnsr->ruleOpsQuad) { err = dFree(tnsr->ruleOpsQuad);dCHK(err); }
  if (tnsr->ruleOpsHex)  { err = dFree(tnsr->ruleOpsHex);dCHK(err); }
  err = dFree(tnsr);dCHK(err);
  dFunctionReturn(0);dCHK(err);
}

static dErr dQuadratureView_Tensor(dQuadrature quad,PetscViewer viewer)
{
  dQuadrature_Tensor *tnsr = (dQuadrature_Tensor*)quad->data;
  dErr err;
  dTruth iascii;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);dCHK(err);
  if (!iascii) SETERRQ(PETSC_ERR_SUP,"only ASCII");
  err = PetscViewerASCIIPrintf(viewer,"Tensor Quadrature: %s\n",GaussFamilies[tnsr->family]);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"alpha %g  beta %g\n",tnsr->alpha,tnsr->beta);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Tensor rules:\n");dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  for (khiter_t k=kh_begin(tnsr->rules); k!= kh_end(tnsr->rules); k++) {
    err = TensorRuleView(kh_val(tnsr->rules,k),viewer);dCHK(err);
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dQuadratureSetFromOptions_Tensor(dQuadrature quad)
{
  dQuadrature_Tensor *tnsr = quad->data;

  dFunctionBegin;
  if (!tnsr) SETERRQ(1,"void");
  dFunctionReturn(0);
}

dErr dQuadratureCreate_Tensor(dQuadrature quad)
{
  static const struct _dQuadratureOps myops = {
    .View = dQuadratureView_Tensor,
    .Destroy = dQuadratureDestroy_Tensor,
    .GetRule = dQuadratureGetRule_Tensor,
    .SetFromOptions = dQuadratureSetFromOptions_Tensor,
  };
  dQuadrature_Tensor *tnsr;
  dErr err;

  dFunctionBegin;
  err = dMemcpy(quad->ops,&myops,sizeof(myops));dCHK(err);
  err = dNewLog(quad,dQuadrature_Tensor,&tnsr);dCHK(err);
  quad->data = tnsr;

  tnsr->rules = kh_init_rule();
  tnsr->family = GAUSS;
  tnsr->alpha = 0.;
  tnsr->beta = 0.;

  err = dQuadratureRuleOpsSetUp_Tensor(quad);dCHK(err);
  dFunctionReturn(0);
}
