#if !defined _TRULE_H
#define _TRULE_H

#include <dohpjacimpl.h>
#include <dohpkhash.h>

typedef struct s_TensorRule *TensorRule;
/**
* We make no particular attempt to align the beginning of the structure, however we would like to align the arrays,
* especially the large derivative arrays, to 16-byte boundaries (to enable SSE).
**/
struct s_TensorRule {
  dInt  size;                      /**< number of quadrature points */
  dReal *weight,*coord;            /**< weights and nodal coordinates */
};

typedef struct {
  dRuleHEADER;
  dEntTopology topo;
  TensorRule trule[3];
} dRule_Tensor;

KHASH_MAP_INIT_INT(tensor,TensorRule)  /* Holds the 1D pieces */
KHASH_MAP_INIT_INT64(rule,dRule_Tensor*) /* Holds n-dimensional rules (basically just aggregates of 1D rules) */

typedef struct {
  dReal alpha,beta;
  dGaussFamily family;
  dQuadratureMethod method;
  khash_t(tensor) *tensor;
  khash_t(rule)   *rules;
  struct _dRuleOps *ruleOpsLine,*ruleOpsQuad,*ruleOpsHex;
} dQuadrature_Tensor;

extern dErr TensorRuleView(TensorRule,PetscViewer); /* exported so that topology implementation can use it */
extern dErr dQuadratureRuleOpsSetUp_Tensor(dQuadrature quad);

#endif
