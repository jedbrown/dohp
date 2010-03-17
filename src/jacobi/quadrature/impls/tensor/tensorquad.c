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

/* This helper is really a hack. */
static void two_point_private(double nodes[],double weights[],int npoints,double dUNUSED alpha,double dUNUSED beta,double x0,double x1)
{
  int nelems = npoints/2;
  double intervals[nelems+1],unused[nelems+1];
  zwglj(intervals,unused,nelems+1,0,0); /* Get Gauss-Lobatto points (this part is the horrible hack) */
  for (int i=0; i<nelems; i++) {
    double h = intervals[i+1] - intervals[i];
    nodes[2*i+0] = intervals[i] + (x0 - (-1))*h/2.;
    nodes[2*i+1] = intervals[i] + (x1 - (-1))*h/2.;
    weights[2*i+0] = h/2.;
    weights[2*i+1] = h/2.;
  }
}

static void two_point_gauss(double nodes[],double weights[],int npoints,double alpha,double beta)
{two_point_private(nodes,weights,npoints,alpha,beta,-0.57735026918962573,0.57735026918962573);}
static void two_point_lobatto(double nodes[],double weights[],int npoints,double alpha,double beta)
{two_point_private(nodes,weights,npoints,alpha,beta,-1,1);}

static dErr TensorGetRule(dQuadrature_Tensor *tnsr,dQuadratureMethod method,dInt order,TensorRule *rule)
{
  dErr err;
  int key,new;
  khiter_t kiter;

  dFunctionBegin;
  *rule = 0;
  if (!(0 <= order && order < 100)) dERROR(1,"rule order %d out of bounds",order);
  key = ((uint32_t)method) << 8 | order;
  kiter = kh_put_tensor(tnsr->tensor,key,&new);
  if (new) {
    TensorRule r;
    void (*nodes_and_weights)(double[],double[],int,double,double);
    dInt size;

    switch (method) {
      case dQUADRATURE_METHOD_SPARSE:
        /* The meaning of this option is somewhat unfortunate, but I don't see a clearly better API at the moment.  Here
         * we just assume a Legendre-Gauss-Lobatto tensor product basis, and construct a quadrature (either Gauss or
         * Lobatto) on the subelements.  A better system would somehow have the subelement details available explicitly.
         */
        size = order;
        if (tnsr->alpha != 0 || tnsr->beta != 0) dERROR(PETSC_ERR_SUP,"only alpha=0, beta=0 (Legendre)");
        switch (tnsr->family) {
          case dGAUSS_GAUSS: nodes_and_weights = two_point_gauss; break;
          case dGAUSS_LOBATTO: nodes_and_weights = two_point_lobatto; break;
          default: dERROR(PETSC_ERR_SUP,"GaussFamily %d",tnsr->family);
        }
        break;
      case dQUADRATURE_METHOD_FAST:
        switch (tnsr->family) {
          case dGAUSS_GAUSS:
            size = 1 + order/2;       /* Gauss integrates a polynomial of order 2*size-1 exactly */
            nodes_and_weights = zwgj;
            break;
          case dGAUSS_LOBATTO:
            size = 2 + order/2;     /* Lobatto integrates a polynomial of order 2*size-3 exactly */
            nodes_and_weights = zwglj;
            break;
          default:
            dERROR(1,"GaussFamily %d not supported",tnsr->family);
        }
        break;
      default: dERROR(PETSC_ERR_SUP,"quadrature method %d",method);
    }
    err = dNew(struct s_TensorRule,&r);dCHK(err);
    err = dMallocA2(size,&r->weight,size,&r->coord);dCHK(err);
    r->size = size;
    if (size == 1) {              /* Polylib function fails for this size. */
      r->weight[0] = 2.0;
      r->coord[0] = 0.0;
    } else {
      nodes_and_weights(r->coord,r->weight,size,tnsr->alpha,tnsr->beta); /* polylib function */
    }
    kh_val(tnsr->tensor,kiter) = r;
  }
  *rule = kh_val(tnsr->tensor,kiter);
  dFunctionReturn(0);
}

static dErr dQuadratureGetRule_Tensor_Private(dQuadrature quad,dInt n,const dEntTopology topo[],const dPolynomialOrder order[],dRule rules[],dQuadratureMethod method)
{
  dQuadrature_Tensor *tnsr = quad->data;
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<n; i++) {
    int       new;
    khint64_t key   = ((uint64_t)topo[i]) << 32 | ((uint64_t)order[i]);
    khiter_t   kiter = kh_put_rule(tnsr->rules,key,&new);
    if (new) {
      dRule_Tensor *newrule;
      err = dNewLog(quad,dRule_Tensor,&newrule);dCHK(err);
      newrule->topo = topo[i];
      switch (topo[i]) {
        case dTOPO_LINE:
          err = dMemcpy(&newrule->ops,tnsr->ruleOpsLine,sizeof(struct _dRuleOps));dCHK(err);
          err = TensorGetRule(tnsr,method,dPolynomialOrder1D(order[i],0),&newrule->trule[0]);dCHK(err);
          break;
        case dTOPO_QUAD:
          err = dMemcpy(&newrule->ops,tnsr->ruleOpsQuad,sizeof(struct _dRuleOps));dCHK(err);
          err = TensorGetRule(tnsr,method,dPolynomialOrder1D(order[i],0),&newrule->trule[0]);dCHK(err);
          err = TensorGetRule(tnsr,method,dPolynomialOrder1D(order[i],1),&newrule->trule[1]);dCHK(err);
          break;
        case dTOPO_HEX:
          err = dMemcpy(&newrule->ops,tnsr->ruleOpsHex,sizeof(struct _dRuleOps));dCHK(err);
          err = TensorGetRule(tnsr,method,dPolynomialOrder1D(order[i],0),&newrule->trule[0]);dCHK(err);
          err = TensorGetRule(tnsr,method,dPolynomialOrder1D(order[i],1),&newrule->trule[1]);dCHK(err);
          err = TensorGetRule(tnsr,method,dPolynomialOrder1D(order[i],2),&newrule->trule[2]);dCHK(err);
          break;
        default: dERROR(1,"no rule available for given topology");
      }
      kh_val(tnsr->rules,kiter) = newrule;
    }
    rules[i] = (dRule)kh_val(tnsr->rules,kiter);
  }
  dFunctionReturn(0);
}

static dErr dQuadratureGetRule_Tensor_FAST(dQuadrature quad,dInt n,const dEntTopology topo[],const dPolynomialOrder order[],dRule rules[])
{return dQuadratureGetRule_Tensor_Private(quad,n,topo,order,rules,dQUADRATURE_METHOD_FAST);}
static dErr dQuadratureGetRule_Tensor_SPARSE(dQuadrature quad,dInt n,const dEntTopology topo[],const dPolynomialOrder order[],dRule rules[])
{return dQuadratureGetRule_Tensor_Private(quad,n,topo,order,rules,dQUADRATURE_METHOD_SPARSE);}

static dErr dQuadratureDestroy_Tensor(dQuadrature quad)
{
  dQuadrature_Tensor *tnsr = (dQuadrature_Tensor*)quad->data;
  dErr err;
  khiter_t k;

  dFunctionBegin;
  /* Destroy all the cached (n-dimensional) rules */
  for (k=kh_begin(tnsr->rules); k!=kh_end(tnsr->rules); k++) {
    if (!kh_exist(tnsr->rules,k)) continue;
    err = dFree(kh_val(tnsr->rules,k));dCHK(err);
  }
  kh_destroy_rule(tnsr->rules);
  /* Destroy all cached 1D tensor rules */
  for (k=kh_begin(tnsr->tensor); k!= kh_end(tnsr->tensor); k++) {
    if (!kh_exist(tnsr->tensor,k)) continue;
    err = TensorRuleDestroy(kh_val(tnsr->tensor,k));dCHK(err);
  }
  kh_destroy_tensor(tnsr->tensor);
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
  err = PetscViewerASCIIPrintf(viewer,"Tensor Quadrature: %s\n",dGaussFamilies[tnsr->family]);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"alpha %g  beta %g\n",tnsr->alpha,tnsr->beta);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Tensor rules:\n");dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  for (khiter_t k=kh_begin(tnsr->tensor); k!= kh_end(tnsr->tensor); k++) {
    err = TensorRuleView(kh_val(tnsr->tensor,k),viewer);dCHK(err);
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dQuadratureSetMethod_Tensor(dQuadrature quad,dQuadratureMethod method)
{
  dQuadrature_Tensor *tnsr = quad->data;

  dFunctionBegin;
  tnsr->method = method;
  switch (method) {
    case dQUADRATURE_METHOD_FAST:
      quad->ops->GetRule = dQuadratureGetRule_Tensor_FAST; break;
    case dQUADRATURE_METHOD_SPARSE:
      quad->ops->GetRule = dQuadratureGetRule_Tensor_SPARSE; break;
    default: dERROR(PETSC_ERR_SUP,"Quadrature method '%s'",dQuadratureMethods[method]);
  }
  dFunctionReturn(0);
}

static dErr dQuadratureSetFromOptions_Tensor(dQuadrature quad)
{
  dQuadrature_Tensor *tnsr  = quad->data;
  dQuadratureMethod  method = dQUADRATURE_METHOD_FAST;
  dTruth             flg;
  dErr               err;

  dFunctionBegin;
  err = PetscOptionsHead("Quadrature Tensor Options");dCHK(err);
  err = PetscOptionsEnum("-dquad_tensor_method","Quadrature method","dQuadratureSetMethod",dQuadratureMethods,tnsr->method,(PetscEnum*)&method,&flg);dCHK(err);
  if (flg || !quad->ops->GetRule) {
    err = dQuadratureSetMethod(quad,method);dCHK(err);
  }
  err = PetscOptionsEnum("-dquad_tensor_gauss_family","Gauss type","None",dGaussFamilies,tnsr->family,(PetscEnum*)&tnsr->family,NULL);dCHK(err);
  err = PetscOptionsTail();dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureCreate_Tensor(dQuadrature quad)
{
  static const struct _dQuadratureOps myops = {
    .View           = dQuadratureView_Tensor,
    .Destroy        = dQuadratureDestroy_Tensor,
    .GetRule        = 0,        /* Does not exist until method is set */
    .SetFromOptions = dQuadratureSetFromOptions_Tensor,
    .SetMethod      = dQuadratureSetMethod_Tensor,
  };
  dQuadrature_Tensor *tnsr;
  dErr err;

  dFunctionBegin;
  err = dMemcpy(quad->ops,&myops,sizeof(myops));dCHK(err);
  err = dNewLog(quad,dQuadrature_Tensor,&tnsr);dCHK(err);
  quad->data = tnsr;

  tnsr->tensor = kh_init_tensor();
  tnsr->rules  = kh_init_rule();
  tnsr->family = dGAUSS_GAUSS;
  tnsr->alpha = 0.;
  tnsr->beta = 0.;

  err = dQuadratureRuleOpsSetUp_Tensor(quad);dCHK(err);
  dFunctionReturn(0);
}
