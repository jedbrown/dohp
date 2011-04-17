static const char help[] = "Solve nonlinear elasticity using dual order hp elements.\n"
  "The model problem is\n"
  "  -div(lambda tr(Du)I + 2 mu (Du - gamma*dot(grad(U),grad(U))/2)) = f \n"
  "where\n"
  "  D is the symmetric gradient operator\n"
  "  mu and lambda are the Lame parameters\n"
  "  gamma controls the nonlinearity, gamma=1 is the standard finite strain tensor and gamma=0 is linear elasticity\n\n";

#include <petscsnes.h>
#include <dohpfs.h>
#include <dohpviewer.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>

static PetscLogEvent LOG_ElastShellMult;

struct ElastParam {
  dReal lambda,mu,gamma;
  dBool bdy100;
};

struct ElastExactCtx {
  dReal a,b,c;
};
struct ElastExact {
  void (*solution)(const struct ElastExactCtx*,const struct ElastParam*,const dReal x[3],dScalar u[3],dScalar du[9]);
  void (*forcing)(const struct ElastExactCtx*,const struct ElastParam*,const dReal x[3],dScalar f[3]);
};

/* Exact solutions are generated using sympy:
from sympy import *
from sympy.abc import *
xyz = [x,y,z]
u0 = cos(a*x) * exp(b*y) * sin(c*z); u1 = sin(a*x) * tanh(b*y) * cosh(c*z); u2 = exp(a*x) * sinh(b*y) * log(1+(c*z)**2); u = [u0,u1,u2]
du = [[diff(v,x),diff(v,y),diff(v,z)] for v in u]
s = [ [l*(du[0][0]+du[1][1]+du[2][2])*(i==j) + mu*(du[i][j] + du[j][i]) + gamma*np.sum([du[k][i]*du[k][j] for k in range(3)]) for j in range(3)] for i in range(3) ]
[ccode(-np.sum([diff(s[i][j],xyz[j]) for j in range(3)])) for i in range(3)]
*/
static void ElastExact_0_Solution(const struct ElastExactCtx *ctx,const struct ElastParam dUNUSED *prm,const dReal xyz[3],dScalar u[3],dScalar du_flat[9])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2];
  dScalar (*du)[3] = (dScalar(*)[3])du_flat;
  u[0] = cos(a*x) * exp(b*y) * sin(c*z);
  u[1] = sin(a*x) * tanh(b*y) * cosh(c*z);
  u[2] = exp(a*x) * sinh(b*y) * log(1+dSqr(c*z));
  du[0][0] = -a*exp(b*y)*sin(a*x)*sin(c*z);
  du[0][1] = b*cos(a*x)*exp(b*y)*sin(c*z);
  du[0][2] = c*cos(a*x)*cos(c*z)*exp(b*y);
  du[1][0] = a*cos(a*x)*cosh(c*z)*tanh(b*y);
  du[1][1] = b*(1 - pow(tanh(b*y),2))*cosh(c*z)*sin(a*x);
  du[1][2] = c*sin(a*x)*sinh(c*z)*tanh(b*y);
  du[2][0] = a*exp(a*x)*log(1 + pow(c,2)*pow(z,2))*sinh(b*y);
  du[2][1] = b*cosh(b*y)*exp(a*x)*log(1 + pow(c,2)*pow(z,2));
  du[2][2] = 2*z*pow(c,2)*exp(a*x)*sinh(b*y)/(1 + pow(c,2)*pow(z,2));
}
static void ElastExact_0_Forcing(const struct ElastExactCtx *ctx,const struct ElastParam *prm,const dReal xyz[3],dScalar f[3])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2],l = prm->lambda,mu = prm->mu,gamma = prm->gamma;
  f[0] = -gamma*(-2*pow(a,3)*pow(cosh(c*z),2)*pow(tanh(b*y),2)*cos(a*x)*sin(a*x) + 2*pow(a,3)*pow(sin(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) + 2*pow(a,3)*pow(log(1 + pow(c,2)*pow(z,2)),2)*pow(sinh(b*y),2)*exp(2*a*x)) - gamma*(a*pow(b,2)*pow((1 - pow(tanh(b*y),2)),2)*pow(cosh(c*z),2)*cos(a*x)*sin(a*x) - 2*a*pow(b,2)*pow(sin(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) - 2*a*pow(b,2)*pow(cosh(c*z),2)*pow(tanh(b*y),2)*(1 - pow(tanh(b*y),2))*cos(a*x)*sin(a*x) + a*pow(b,2)*pow(cosh(b*y),2)*pow(log(1 + pow(c,2)*pow(z,2)),2)*exp(2*a*x) + a*pow(b,2)*pow(log(1 + pow(c,2)*pow(z,2)),2)*pow(sinh(b*y),2)*exp(2*a*x)) - gamma*(a*pow(c,2)*pow(cosh(c*z),2)*pow(tanh(b*y),2)*cos(a*x)*sin(a*x) + a*pow(c,2)*pow(sin(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) + a*pow(c,2)*pow(sinh(c*z),2)*pow(tanh(b*y),2)*cos(a*x)*sin(a*x) - a*pow(c,2)*pow(cos(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) + 2*a*pow(c,2)*pow(sinh(b*y),2)*exp(2*a*x)*log(1 + pow(c,2)*pow(z,2))/(1 + pow(c,2)*pow(z,2)) - 4*a*pow(c,4)*pow(z,2)*pow(sinh(b*y),2)*exp(2*a*x)*log(1 + pow(c,2)*pow(z,2))/pow((1 + pow(c,2)*pow(z,2)),2) + 4*a*pow(c,4)*pow(z,2)*pow(sinh(b*y),2)*exp(2*a*x)/pow((1 + pow(c,2)*pow(z,2)),2)) - l*(-pow(a,2)*cos(a*x)*exp(b*y)*sin(c*z) + a*b*(1 - pow(tanh(b*y),2))*cos(a*x)*cosh(c*z) + 2*a*z*pow(c,2)*exp(a*x)*sinh(b*y)/(1 + pow(c,2)*pow(z,2))) - mu*(pow(b,2)*cos(a*x)*exp(b*y)*sin(c*z) + a*b*(1 - pow(tanh(b*y),2))*cos(a*x)*cosh(c*z)) - mu*(-pow(c,2)*cos(a*x)*exp(b*y)*sin(c*z) + 2*a*z*pow(c,2)*exp(a*x)*sinh(b*y)/(1 + pow(c,2)*pow(z,2))) + 2*mu*pow(a,2)*cos(a*x)*exp(b*y)*sin(c*z);
  f[1] = -gamma*(-4*pow(b,3)*pow((1 - pow(tanh(b*y),2)),2)*pow(cosh(c*z),2)*pow(sin(a*x),2)*tanh(b*y) + 2*pow(b,3)*pow(log(1 + pow(c,2)*pow(z,2)),2)*cosh(b*y)*exp(2*a*x)*sinh(b*y) + 2*pow(b,3)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)) - gamma*(b*pow(a,2)*pow(cos(a*x),2)*pow(cosh(c*z),2)*(1 - pow(tanh(b*y),2))*tanh(b*y) - b*pow(a,2)*pow(cosh(c*z),2)*pow(sin(a*x),2)*(1 - pow(tanh(b*y),2))*tanh(b*y) + 2*b*pow(a,2)*pow(log(1 + pow(c,2)*pow(z,2)),2)*cosh(b*y)*exp(2*a*x)*sinh(b*y) + b*pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y) - b*pow(a,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)) - gamma*(b*pow(c,2)*pow(cosh(c*z),2)*pow(sin(a*x),2)*(1 - pow(tanh(b*y),2))*tanh(b*y) + b*pow(c,2)*pow(sin(a*x),2)*pow(sinh(c*z),2)*(1 - pow(tanh(b*y),2))*tanh(b*y) + 2*b*pow(c,2)*cosh(b*y)*exp(2*a*x)*log(1 + pow(c,2)*pow(z,2))*sinh(b*y)/(1 + pow(c,2)*pow(z,2)) + 4*b*pow(c,4)*pow(z,2)*cosh(b*y)*exp(2*a*x)*sinh(b*y)/pow((1 + pow(c,2)*pow(z,2)),2) - 4*b*pow(c,4)*pow(z,2)*cosh(b*y)*exp(2*a*x)*log(1 + pow(c,2)*pow(z,2))*sinh(b*y)/pow((1 + pow(c,2)*pow(z,2)),2) + b*pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y) - b*pow(c,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)) - l*(-a*b*exp(b*y)*sin(a*x)*sin(c*z) - 2*pow(b,2)*(1 - pow(tanh(b*y),2))*cosh(c*z)*sin(a*x)*tanh(b*y) + 2*b*z*pow(c,2)*cosh(b*y)*exp(a*x)/(1 + pow(c,2)*pow(z,2))) - mu*(pow(c,2)*cosh(c*z)*sin(a*x)*tanh(b*y) + 2*b*z*pow(c,2)*cosh(b*y)*exp(a*x)/(1 + pow(c,2)*pow(z,2))) - mu*(-pow(a,2)*cosh(c*z)*sin(a*x)*tanh(b*y) - a*b*exp(b*y)*sin(a*x)*sin(c*z)) + 4*mu*pow(b,2)*(1 - pow(tanh(b*y),2))*cosh(c*z)*sin(a*x)*tanh(b*y);
  f[2] = -gamma*(-2*pow(c,3)*pow(cos(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) + 2*pow(c,3)*pow(sin(a*x),2)*pow(tanh(b*y),2)*cosh(c*z)*sinh(c*z) - 16*pow(c,6)*pow(z,3)*pow(sinh(b*y),2)*exp(2*a*x)/pow((1 + pow(c,2)*pow(z,2)),3) + 8*z*pow(c,4)*pow(sinh(b*y),2)*exp(2*a*x)/pow((1 + pow(c,2)*pow(z,2)),2)) - gamma*(c*pow(a,2)*pow(cos(a*x),2)*pow(tanh(b*y),2)*cosh(c*z)*sinh(c*z) + c*pow(a,2)*pow(sin(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) - c*pow(a,2)*pow(cos(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) - c*pow(a,2)*pow(sin(a*x),2)*pow(tanh(b*y),2)*cosh(c*z)*sinh(c*z) + 4*z*pow(a,2)*pow(c,2)*pow(sinh(b*y),2)*exp(2*a*x)*log(1 + pow(c,2)*pow(z,2))/(1 + pow(c,2)*pow(z,2))) - gamma*(c*pow(b,2)*pow((1 - pow(tanh(b*y),2)),2)*pow(sin(a*x),2)*cosh(c*z)*sinh(c*z) + 2*c*pow(b,2)*pow(cos(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) - 2*c*pow(b,2)*pow(sin(a*x),2)*pow(tanh(b*y),2)*(1 - pow(tanh(b*y),2))*cosh(c*z)*sinh(c*z) + 2*z*pow(b,2)*pow(c,2)*pow(cosh(b*y),2)*exp(2*a*x)*log(1 + pow(c,2)*pow(z,2))/(1 + pow(c,2)*pow(z,2)) + 2*z*pow(b,2)*pow(c,2)*pow(sinh(b*y),2)*exp(2*a*x)*log(1 + pow(c,2)*pow(z,2))/(1 + pow(c,2)*pow(z,2))) - l*(2*pow(c,2)*exp(a*x)*sinh(b*y)/(1 + pow(c,2)*pow(z,2)) + b*c*(1 - pow(tanh(b*y),2))*sin(a*x)*sinh(c*z) - a*c*cos(c*z)*exp(b*y)*sin(a*x) - 4*pow(c,4)*pow(z,2)*exp(a*x)*sinh(b*y)/pow((1 + pow(c,2)*pow(z,2)),2)) - mu*(pow(a,2)*exp(a*x)*log(1 + pow(c,2)*pow(z,2))*sinh(b*y) - a*c*cos(c*z)*exp(b*y)*sin(a*x)) - mu*(pow(b,2)*exp(a*x)*log(1 + pow(c,2)*pow(z,2))*sinh(b*y) + b*c*(1 - pow(tanh(b*y),2))*sin(a*x)*sinh(c*z)) - 4*mu*pow(c,2)*exp(a*x)*sinh(b*y)/(1 + pow(c,2)*pow(z,2)) + 8*mu*pow(c,4)*pow(z,2)*exp(a*x)*sinh(b*y)/pow((1 + pow(c,2)*pow(z,2)),2);
}

static void ElastExact_1_Solution(const struct ElastExactCtx *ctx,const struct ElastParam dUNUSED *prm,const dReal xyz[3],dScalar u[3],dScalar du_flat[9])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2];
  dScalar (*du)[3] = (dScalar(*)[3])du_flat;
  u[0] = a*x;
  u[1] = b*y;
  u[2] = c*z;
  du[0][0] = a;
  du[0][1] = 0;
  du[0][2] = 0;
  du[1][0] = 0;
  du[1][1] = b;
  du[1][2] = 0;
  du[2][0] = 0;
  du[2][1] = 0;
  du[2][2] = c;
}
static void ElastExact_1_Forcing(const struct ElastExactCtx dUNUSED *ctx,const struct ElastParam dUNUSED *prm,const dReal dUNUSED xyz[3],dScalar f[3])
{
  //const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2],l = prm->lambda,mu = prm->mu,gamma = prm->gamma;
  f[0] = 0;
  f[1] = 0;
  f[2] = 0;
}

static void ElastExact_2_Solution(const struct ElastExactCtx dUNUSED *ctx,const struct ElastParam dUNUSED *prm,const dReal dUNUSED xyz[3],dScalar u[3],dScalar du_flat[9])
{
  dScalar (*du)[3] = (dScalar(*)[3])du_flat;
  u[0] = 0;
  u[1] = 0;
  u[2] = 0;
  du[0][0] = 0;
  du[0][1] = 0;
  du[0][2] = 0;
  du[1][0] = 0;
  du[1][1] = 0;
  du[1][2] = 0;
  du[2][0] = 0;
  du[2][1] = 0;
  du[2][2] = 0;
}
static void ElastExact_2_Forcing(const struct ElastExactCtx dUNUSED *ctx,const struct ElastParam dUNUSED *prm,const dReal dUNUSED xyz[3],dScalar f[3])
{
  f[0] = 0;
  f[1] = 1;
  f[2] = 0;
}


struct ElastStore {
  dReal Du[3][3];
};

typedef enum {EVAL_FUNCTION,EVAL_JACOBIAN, EVAL_UB} ElastEvaluation;

typedef struct ElastCtx *Elast;
struct ElastCtx {
  MPI_Comm              comm;
  struct ElastParam     param;
  struct ElastExact     exact;
  struct ElastExactCtx  exactctx;
  dJacobi               jac;
  dMesh                 mesh;
  dFS                   fs;
  Vec                   x,y;
  dInt                  constBDeg;
  dBool                 errorview;
  dQuadratureMethod     function_qmethod,jacobian_qmethod;
  dRulesetIterator      regioniter[EVAL_UB];
};

static dErr ElastCreate(MPI_Comm comm,Elast *elast)
{
  Elast elt;
  dErr err;

  dFunctionBegin;
  *elast = 0;
  err = dNew(struct ElastCtx,&elt);dCHK(err);
  elt->comm = comm;

  elt->constBDeg = 3;
  elt->param.lambda = 10;
  elt->param.mu     = 100;
  elt->param.gamma  = 0;
  elt->param.bdy100 = dFALSE;
  elt->function_qmethod = dQUADRATURE_METHOD_FAST;
  elt->jacobian_qmethod = dQUADRATURE_METHOD_SPARSE;

  *elast = elt;
  dFunctionReturn(0);
}

static dErr ElastSetFromOptions(Elast elt)
{
  struct ElastParam *prm = &elt->param;
  struct ElastExactCtx *exc = &elt->exactctx;
  dMesh mesh;
  dFS fs;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag dtag;
  dInt exact;
  dErr err;

  dFunctionBegin;
  exact = 0; exc->a = exc->b = exc->c = 1;
  err = PetscOptionsBegin(elt->comm,NULL,"Elasticity options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-const_bdeg","Use constant isotropic degree on all elements","",elt->constBDeg,&elt->constBDeg,NULL);dCHK(err);
    err = PetscOptionsBool("-error_view","View errors","",elt->errorview,&elt->errorview,NULL);dCHK(err);
    err = PetscOptionsReal("-elast_lambda","first Lame parameter","",prm->lambda,&prm->lambda,NULL);dCHK(err);
    err = PetscOptionsReal("-elast_mu","Second Lame parameter","",prm->mu,&prm->mu,NULL);dCHK(err);
    err = PetscOptionsReal("-elast_gamma","Strength of nonlinearity [0,1]","",prm->gamma,&prm->gamma,NULL);dCHK(err);
    err = PetscOptionsBool("-bdy100","Only use boundary 100","",prm->bdy100,&prm->bdy100,NULL);dCHK(err);
    err = PetscOptionsEnum("-elast_f_qmethod","Quadrature method for residual evaluation/matrix-free","",dQuadratureMethods,(PetscEnum)elt->function_qmethod,(PetscEnum*)&elt->function_qmethod,NULL);dCHK(err);
    err = PetscOptionsEnum("-elast_jac_qmethod","Quadrature to use for Jacobian assembly","",dQuadratureMethods,(PetscEnum)elt->jacobian_qmethod,(PetscEnum*)&elt->jacobian_qmethod,NULL);dCHK(err);
    err = PetscOptionsInt("-exact","Exact solution choice","",exact,&exact,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_a","First scale parameter","",exc->a,&exc->a,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_b","Second scale parameter","",exc->b,&exc->b,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_c","Third scale parameter","",exc->c,&exc->c,NULL);dCHK(err);

  } err = PetscOptionsEnd();dCHK(err);

  switch (exact) {
    case 0:
      elt->exact.solution = ElastExact_0_Solution;
      elt->exact.forcing = ElastExact_0_Forcing;
      break;
    case 1:
      elt->exact.solution = ElastExact_1_Solution;
      elt->exact.forcing = ElastExact_1_Forcing;
      break;
    case 2:
      elt->exact.solution = ElastExact_2_Solution;
      elt->exact.forcing = ElastExact_2_Forcing;
      break;
    default: dERROR(PETSC_COMM_SELF,1,"Exact solution %d not implemented");
  }

  err = dMeshCreate(elt->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);dCHK(err);
  elt->mesh = mesh;
  err = dMeshGetRoot(mesh,&domain);dCHK(err); /* Need a taggable set */
  err = dMeshSetDuplicateEntsOnly(mesh,domain,&domain);dCHK(err);

  err = dJacobiCreate(elt->comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  elt->jac = jac;

  err = dMeshCreateRuleTagIsotropic(mesh,domain,"elast_efs_degree",elt->constBDeg,&dtag);dCHK(err);

  err = dFSCreate(elt->comm,&fs);dCHK(err);
  err = dFSSetBlockSize(fs,3);dCHK(err);
  err = dFSSetMesh(fs,mesh,domain);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  err = dFSRegisterBoundary(fs,100,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  if (!elt->param.bdy100) {
    err = dFSRegisterBoundary(fs,200,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
    err = dFSRegisterBoundary(fs,300,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  }
  err = dFSSetFromOptions(fs);dCHK(err);
  err = PetscObjectSetName((PetscObject)fs,"dFS_0");dCHK(err);
  err = PetscObjectSetName((PetscObject)mesh,"dMesh_0");dCHK(err);
  elt->fs = fs;

  err = dFSCreateExpandedVector(fs,&elt->x);dCHK(err);
  err = VecDuplicate(elt->x,&elt->y);dCHK(err);
  dFunctionReturn(0);
}

static dErr ElastGetRegionIterator(Elast elt,ElastEvaluation eval,dRulesetIterator *riter)
{
  dErr err;

  dFunctionBegin;
  if (!elt->regioniter[eval]) {
    dRulesetIterator iter;
    dRuleset ruleset;
    dFS cfs;
    dMeshESH domain;
    dQuadratureMethod qmethod;
    switch (eval) {
    case EVAL_FUNCTION: qmethod = elt->function_qmethod; break;
    case EVAL_JACOBIAN: qmethod = elt->jacobian_qmethod; break;
    default: dERROR(elt->comm,PETSC_ERR_ARG_OUTOFRANGE,"Unknown evaluation context");
    }
    err = dFSGetDomain(elt->fs,&domain);dCHK(err);
    err = dFSGetPreferredQuadratureRuleSet(elt->fs,domain,dTYPE_REGION,dTOPO_ALL,qmethod,&ruleset);dCHK(err);
    err = dFSGetCoordinateFS(elt->fs,&cfs);dCHK(err);
    err = dRulesetCreateIterator(ruleset,cfs,&iter);dCHK(err);
    err = dRulesetDestroy(ruleset);dCHK(err); /* Give ownership to iterator */
    err = dRulesetIteratorAddFS(iter,elt->fs);dCHK(err);
    if (eval == EVAL_FUNCTION) {err = dRulesetIteratorAddStash(iter,0,sizeof(struct ElastStore));dCHK(err);}
    elt->regioniter[eval] = iter;
  }
  *riter = elt->regioniter[eval];
  dFunctionReturn(0);
}

static dErr ElastDestroy(Elast elt)
{
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(elt->fs);dCHK(err);
  err = dJacobiDestroy(elt->jac);dCHK(err);
  err = dMeshDestroy(elt->mesh);dCHK(err);
  if (elt->x) {err = VecDestroy(elt->x);dCHK(err);}
  if (elt->y) {err = VecDestroy(elt->y);dCHK(err);}
  for (dInt i=0; i<EVAL_UB; i++) {
    if (elt->regioniter[i]) {dRulesetIteratorDestroy(elt->regioniter[i]);dCHK(err);}
  }
  err = dFree(elt);dCHK(err);
  dFunctionReturn(0);
}

static inline void ElastPointwiseComputeStore(struct ElastParam dUNUSED *prm,const dReal dUNUSED x[3],const dScalar dUNUSED u[3],const dScalar Du[9],struct ElastStore *st)
{
  memcpy(st->Du,Du,9*sizeof Du[0]);
}

static inline dScalar DotColumn(const dScalar a[3][3],const dScalar b[3][3],dInt i,dInt j)
{return a[0][i]*b[0][j] + a[1][i]*b[1][j] + a[2][i]*b[2][j];}

static inline void ElastPointwiseFunction(struct ElastParam *prm,struct ElastExact *exact,struct ElastExactCtx *exactctx,
                                          const dReal x[3],dReal weight,const dScalar u[3],const dScalar Du_flat[9],
                                          struct ElastStore *st,dScalar v[3],dScalar Dv_flat[9])
{
  const dScalar (*restrict Du)[3] = (const dScalar(*)[3])Du_flat;
  dScalar       (*restrict Dv)[3] = (dScalar(*)[3])Dv_flat;
  dScalar f[3],E0,E1,E2,E3,E4,E5,volume;
  ElastPointwiseComputeStore(prm,x,u,&Du[0][0],st);
  exact->forcing(exactctx,prm,x,f);
  v[0] = -weight * f[0];       /* Coefficient of \a v in weak form */
  v[1] = -weight * f[1];
  v[2] = -weight * f[2];
                                /* Twice the Green-Lagrangian strain tensor */
  E0 = 2*Du[0][0] + prm->gamma*DotColumn(Du,Du,0,0);
  E1 = 2*Du[1][1] + prm->gamma*DotColumn(Du,Du,1,1);
  E2 = 2*Du[2][2] + prm->gamma*DotColumn(Du,Du,2,2);
  E3 = Du[0][1] + Du[1][0] + prm->gamma*DotColumn(Du,Du,0,1);
  E4 = Du[0][2] + Du[2][0] + prm->gamma*DotColumn(Du,Du,0,2);
  E5 = Du[1][2] + Du[2][1] + prm->gamma*DotColumn(Du,Du,1,2);
                                /* volume strain, should use \f$ e = 3\det{Du}^{1/3} \f$ */
  volume = Du[0][0] + Du[1][1] + Du[2][2];
                                /* Coefficients in weak form of deformation gradient of test functions */
  Dv[0][0] = weight*(prm->lambda*volume + prm->mu*E0);
  Dv[1][1] = weight*(prm->lambda*volume + prm->mu*E1);
  Dv[2][2] = weight*(prm->lambda*volume + prm->mu*E2);
  Dv[0][1] = Dv[1][0] = weight*prm->mu*E3;
  Dv[0][2] = Dv[2][0] = weight*prm->mu*E4;
  Dv[1][2] = Dv[2][1] = weight*prm->mu*E5;
}

static inline void ElastPointwiseJacobian(struct ElastParam dUNUSED *prm,const struct ElastStore *restrict st,dReal weight,
                                          const dScalar dUNUSED u[restrict static 3],const dScalar Du_flat[restrict static 9],
                                          dScalar v[restrict static 3],dScalar Dv_flat[restrict static 9])
{
  const dScalar (*Du)[3] = (const dScalar(*restrict)[3])Du_flat;
  dScalar (*Dv)[3] = (dScalar (*)[3])Dv_flat;
  const dScalar volume = Du[0][0] + Du[1][1] + Du[2][2]; /* volume strain, should use \f$ e = 3\det{Du}^{1/3} \f$ */
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;
                                /* Coefficients in weak form of Jacobian */
  for (dInt i=0; i<3; i++) {
    for (dInt j=0; j<3; j++) {
      Dv[i][j] = weight*((i==j)*prm->lambda*volume + prm->mu*(Du[i][j] + Du[j][i]
                                                              + prm->gamma*(+ st->Du[0][i]*Du[0][j] + st->Du[0][j]*Du[0][i]
                                                                            + st->Du[1][i]*Du[1][j] + st->Du[1][j]*Du[1][i]
                                                                            + st->Du[2][i]*Du[2][j] + st->Du[2][j]*Du[2][i])));
    }
  }
}

static dErr ElastFunction(SNES dUNUSED snes,Vec gx,Vec gy,void *ctx)
{
  Elast elt = ctx;
  Vec Coords;
  dRulesetIterator iter;
  dErr err;

  dFunctionBegin;
  err = VecZeroEntries(gy);dCHK(err);
  err = dFSGetGeometryVectorExpanded(elt->fs,&Coords);dCHK(err);
  err = ElastGetRegionIterator(elt,EVAL_FUNCTION,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,gy,dFS_INHOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar *x,*dx,*u,*du,*v,*dv;
    dInt Q;
    struct ElastStore *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,&v,&dv);dCHK(err);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      ElastPointwiseFunction(&elt->param,&elt->exact,&elt->exactctx,&x[i*3],jw[i],&u[i*3],&du[i*9],&stash[i],&v[i*3],&dv[i*9]);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES,NULL,NULL,v,dv);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  dFunctionReturn(0);
}

static dErr ElastShellMatMult(Mat J,Vec gx,Vec gy)
{
  Elast elt;
  Vec Coords;
  dRulesetIterator iter;
  dErr err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_ElastShellMult,J,gx,gy,0);dCHK(err);
  err = MatShellGetContext(J,(void**)&elt);dCHK(err);
  err = VecZeroEntries(gy);dCHK(err);
  err = dFSGetGeometryVectorExpanded(elt->fs,&Coords);dCHK(err);
  err = ElastGetRegionIterator(elt,EVAL_FUNCTION,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_HOMOGENEOUS,gy,dFS_HOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar *x,*dx,*u,*du,*v,*dv;
    dInt Q;
    const struct ElastStore *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,&v,&dv);dCHK(err);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      ElastPointwiseJacobian(&elt->param,&stash[i],jw[i],&u[i*3],&du[i*9],&v[i*3],&dv[i*9]);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES,NULL,NULL,v,dv);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = PetscLogEventEnd(LOG_ElastShellMult,J,gx,gy,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr ElastJacobian(SNES dUNUSED snes,Vec gx,Mat *J,Mat *Jp,MatStructure *structure,void *ctx)
{
  Elast            elt = ctx;
  Vec              Coords;
  dRulesetIterator iter;
  dScalar          *Kflat;
  dErr             err;

  dFunctionBegin;
  err = MatZeroEntries(*Jp);dCHK(err);
  err = dFSGetGeometryVectorExpanded(elt->fs,&Coords);dCHK(err);
  err = ElastGetRegionIterator(elt,EVAL_JACOBIAN,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSplit(iter,NULL,NULL,NULL,&Kflat);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw,*interp_flat,*deriv_flat;
    const dInt *rowcol;
    dScalar (*x)[3],(*dx)[3][3],(*u)[3],(*du)[3][3];
    dInt Q,P;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, (dScalar**)&u,(dScalar**)&du,NULL,NULL);dCHK(err);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat);dCHK(err);
    {                           /* Scope so that we can declare new VLA pointers for convenient assembly */
      const dReal (*interp)[P] = (const dReal(*)[P])interp_flat;
      const dReal (*deriv)[P][3] = (const dReal(*)[P][3])deriv_flat;
      dScalar (*K)[3][P][3] = (dScalar(*)[3][P][3])Kflat;
      err = PetscMemzero(K,P*3*P*3*sizeof(K[0][0][0][0]));dCHK(err);
      for (dInt q=0; q<Q; q++) {
        struct ElastStore store;
        ElastPointwiseComputeStore(&elt->param,x[q],u[q],&du[q][0][0],&store);
        for (dInt j=0; j<P; j++) {
          for (dInt fj=0; fj<3; fj++) {
            dScalar uu[3] = {0},duu[3][3] = {{0},{0},{0}},v[3],dv[3][3];
            uu[fj] = interp[q][j];
            duu[fj][0] = deriv[q][j][0];
            duu[fj][1] = deriv[q][j][1];
            duu[fj][2] = deriv[q][j][2];
            ElastPointwiseJacobian(&elt->param,&store,jw[q],uu,&duu[0][0],v,&dv[0][0]);
            for (dInt i=0; i<P; i++) {
              for (dInt fi=0; fi<3; fi++) {
                K[i][fi][j][fj] += (interp[q][i] * v[fi]
                                    + deriv[q][i][0] * dv[fi][0]
                                    + deriv[q][i][1] * dv[fi][1]
                                    + deriv[q][i][2] * dv[fi][2]);
              }
            }
          }
        }
      }
      err = dFSMatSetValuesBlockedExpanded(elt->fs,*Jp,P,rowcol,P,rowcol,&K[0][0][0][0],ADD_VALUES);dCHK(err);
    }
    err = dRulesetIteratorRestorePatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);

  err = MatAssemblyBegin(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  if (J != Jp) {
    err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  }
  *structure = SAME_NONZERO_PATTERN;
  dFunctionReturn(0);
}

static dErr ElastErrorNorms(Elast elt,Vec gx,dReal errorNorms[static 3],dReal gerrorNorms[static 3])
{
  dErr err;
  dInt patchcnt;
  Vec Coords;
  dRulesetIterator iter;

  dFunctionBegin;
  err = dMemzero(errorNorms,3*sizeof(errorNorms));dCHK(err);
  err = dMemzero(gerrorNorms,3*sizeof(gerrorNorms));dCHK(err);
  err = dFSGetGeometryVectorExpanded(elt->fs,&Coords);dCHK(err);
  err = ElastGetRegionIterator(elt,EVAL_FUNCTION,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  patchcnt = 0;
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[3][3],(*u)[3],(*du)[3][3];
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,NULL);dCHK(err);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[3],duu[3][3],r[3],gr[3],rsum=0,grsum=0;
      elt->exact.solution(&elt->exactctx,&elt->param,x[i],uu,&duu[0][0]);
      for (dInt j=0; j<3; j++) {
        r[j] = u[i][j] - uu[j]; /* Function error at point */
        rsum += dSqr(r[j]);
        gr[j] = dSqrt(dSqr(du[i][j][0]-duu[j][0]) + dSqr(du[i][j][1]-duu[j][1]) + dSqr(du[i][j][2]-duu[j][2])); /* Gradient error at point */
        grsum += dSqr(gr[j]);
      }
      if (elt->errorview) {
        printf("e,q = %3d %3d (% 5f,% 5f,% 5f) dohp %10.2e %10.2e %10.2e   exact %10.2e %10.2e %10.2e   error %10.e\n",
               patchcnt,i,x[i][0],x[i][1],x[i][2],u[i][0],u[i][1],u[i][2],uu[0],uu[1],uu[2],rsum);
      }
      errorNorms[0] += (dAbs(r[0]) + dAbs(r[1]) + dAbs(r[2])) * jw[i];                   /* 1-norm */
      errorNorms[1] += grsum * jw[i];                                                    /* 2-norm */
      errorNorms[2] = dMax(errorNorms[2],dMax(dAbs(r[0]),dMax(dAbs(r[1]),dAbs(r[2])))); /* Sup-norm */
      gerrorNorms[0] += (dAbs(gr[0]) + dAbs(gr[1]) + dAbs(gr[2])) * jw[i];
      gerrorNorms[1] += grsum * jw[i];
      gerrorNorms[2] = dMax(gerrorNorms[2],dMax(dAbs(gr[0]),dMax(dAbs(gr[1]),dAbs(gr[2]))));
#if 0
      printf("pointwise stats %8g %8g %8g %8g\n",jw[i],r[0],dSqr(r[0]),errorNorms[1]);
      printf("pointwise grads %8g %8g %8g (%8g)\n",gr[0],gr[1],gr[2],grsum);
# if 0
      printf("jinv[%2d][%3d]   %+3.1f %+3.1f %+3.1f    %+3.1f %+3.1f %+3.1f    %+3.1f %+3.1f %+3.1f\n",e,i,
             jinv[i][0][0],jinv[i][0][1],jinv[i][0][2],
             jinv[i][1][0],jinv[i][1][1],jinv[i][1][2],
             jinv[i][2][0],jinv[i][2][1],jinv[i][2][2]);
# endif
#endif
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
    patchcnt++;
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  errorNorms[1] = dSqrt(errorNorms[1]);
  gerrorNorms[1] = dSqrt(gerrorNorms[1]);
  dFunctionReturn(0);
}

static dErr ElastGetSolutionVector(Elast elt,Vec *insoln)
{
  dErr err;
  Vec soln,xc,cvecg,cvec;
  dScalar *x;
  const dScalar *coords;
  dInt n,bs;

  dFunctionBegin;
  *insoln = 0;
  err = dFSCreateGlobalVector(elt->fs,&soln);dCHK(err);
  err = VecDohpGetClosure(soln,&xc);dCHK(err);
  err = dFSGetNodalCoordinatesGlobal(elt->fs,&cvecg);dCHK(err);
  err = VecDohpGetClosure(cvecg,&cvec);dCHK(err);
  err = VecGetLocalSize(xc,&n);dCHK(err);
  err = VecGetBlockSize(xc,&bs);dCHK(err);
  {
    dInt nc;
    err = VecGetLocalSize(cvec,&nc);dCHK(err);
    if (nc*bs != n*3) dERROR(PETSC_COMM_SELF,1,"Coordinate vector has inconsistent size");
  }
  err = VecGetArray(xc,&x);dCHK(err);
  err = VecGetArrayRead(cvec,&coords);dCHK(err);
  for (dInt i=0; i<n/bs; i++) {
    dScalar du_unused[3*bs];
    elt->exact.solution(&elt->exactctx,&elt->param,&coords[3*i],&x[i*bs],du_unused);
    /* printf("Node %3d: coords %+8f %+8f %+8f   exact %+8f %+8f %+8f\n",i,coords[3*i],coords[3*i+1],coords[3*i+2],x[3*i],x[3*i+1],x[3*i+2]); */
  }
  err = VecRestoreArray(xc,&x);dCHK(err);
  err = VecRestoreArrayRead(cvec,&coords);dCHK(err);
  err = VecDohpRestoreClosure(cvecg,&cvec);dCHK(err);
  err = dFSInhomogeneousDirichletCommit(elt->fs,xc);dCHK(err);
  err = VecDohpRestoreClosure(soln,&xc);dCHK(err);
  *insoln = soln;
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  char mtype[256] = MATBAIJ;
  Elast elt;
  dFS fs;
  MPI_Comm comm;
  PetscViewer viewer;
  Mat J,Jp;
  Vec r,x,soln;
  SNES snes;
  dBool  nojshell = dFALSE,nocheck = dFALSE,view_soln = dFALSE;
  dErr err;

  err = dInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = PetscLogEventRegister("ElastShellMult",MAT_CLASSID,&LOG_ElastShellMult);dCHK(err);

  err = ElastCreate(comm,&elt);dCHK(err);
  err = ElastSetFromOptions(elt);dCHK(err);
  fs = elt->fs;

  err = dFSCreateGlobalVector(fs,&r);dCHK(err);
  err = PetscOptionsGetString(NULL,"-q1mat_type",mtype,sizeof(mtype),NULL);dCHK(err);
  err = dFSGetMatrix(fs,mtype,&Jp);dCHK(err);
  err = MatSetOptionsPrefix(Jp,"q1");dCHK(err);

  err = PetscOptionsBegin(elt->comm,NULL,"Elasticity solver options",__FILE__);dCHK(err); {
    err = PetscOptionsBool("-nojshell","Do not use shell Jacobian","",nojshell,&nojshell,NULL);dCHK(err);
    err = PetscOptionsBool("-nocheck_error","Do not compute errors","",nocheck,&nocheck,NULL);dCHK(err);
    err = PetscOptionsBool("-view_soln","View the solution","",view_soln,&view_soln,NULL);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);
  if (nojshell) {
    /* Use the preconditioning matrix in place of the Jacobin.  This will NOT converge unless the elements are actually
    * Q1 (bdeg=2).  This option is nullified by -snes_mf_operator which will still only use the assembled Jacobian for
    * preconditioning. */
    J = Jp;
  } else {
    dInt m,n,M,N;
    err = MatGetLocalSize(Jp,&m,&n);dCHK(err);
    err = MatGetSize(Jp,&M,&N);dCHK(err);
    err = MatCreateShell(comm,m,n,M,N,elt,&J);dCHK(err);
    err = MatSetOptionsPrefix(J,"j");dCHK(err);
    err = MatShellSetOperation(J,MATOP_MULT,(void(*)(void))ElastShellMatMult);dCHK(err);
  }
  err = SNESCreate(comm,&snes);dCHK(err);
  err = SNESSetFunction(snes,r,ElastFunction,elt);dCHK(err);
  err = SNESSetJacobian(snes,J,Jp,ElastJacobian,elt);dCHK(err);
  err = SNESSetTolerances(snes,PETSC_DEFAULT,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);dCHK(err);
  err = SNESSetFromOptions(snes);dCHK(err);
  err = VecZeroEntries(r);dCHK(err);
  err = VecDuplicate(r,&x);dCHK(err);
  err = PetscObjectSetName((PetscObject)x,"Displacement");dCHK(err);
  err = ElastGetSolutionVector(elt,&soln);dCHK(err);
  {
    Vec sc;
    err = VecDohpGetClosure(soln,&sc);dCHK(err);
    err = dFSInhomogeneousDirichletCommit(elt->fs,sc);dCHK(err);
    err = VecDohpRestoreClosure(soln,&sc);dCHK(err);
  }
  err = VecZeroEntries(x);dCHK(err);
  err = SNESSolve(snes,NULL,x);dCHK(err);
  if (!nocheck) {
    dReal anorm[2],anorminf,inorm[3],enorm[3],gnorm[3];
    err = ElastErrorNorms(elt,x,enorm,gnorm);dCHK(err);
    err = VecNorm(r,NORM_1_AND_2,anorm);dCHK(err);
    err = VecNorm(r,NORM_INFINITY,&anorminf);dCHK(err);
    err = VecWAXPY(r,-1,soln,x);dCHK(err);
    err = VecNorm(r,NORM_1_AND_2,inorm);dCHK(err);
    err = VecNorm(r,NORM_INFINITY,&inorm[2]);dCHK(err);
    err = dPrintf(comm,"Algebraic residual        |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",anorm[0],anorm[1],anorminf);dCHK(err);
    err = dPrintf(comm,"Interpolation residual    |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",inorm[0],inorm[1],inorm[2]);dCHK(err);
    err = dPrintf(comm,"Pointwise solution error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",enorm[0],enorm[1],enorm[2]);dCHK(err);
    err = dPrintf(comm,"Pointwise gradient error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",gnorm[0],gnorm[1],gnorm[2]);dCHK(err);
  }
  if (view_soln) {
    dViewer viewdhm;
    err = PetscViewerCreate(comm,&viewdhm);dCHK(err);
    err = PetscViewerSetType(viewdhm,PETSCVIEWERDHM);dCHK(err);
    err = PetscViewerFileSetName(viewdhm,"elast.dhm");dCHK(err);
    err = PetscViewerFileSetMode(viewdhm,FILE_MODE_WRITE);dCHK(err);
    err = dViewerDHMSetTimeUnits(viewdhm,"hour",PETSC_PI*1e7/3600);dCHK(err);
    err = dViewerDHMSetTime(viewdhm,0.1);dCHK(err);
    err = VecView(x,viewdhm);dCHK(err);
    err = PetscViewerDestroy(viewdhm);dCHK(err);
  }

  err = VecDestroy(r);dCHK(err);
  err = VecDestroy(x);dCHK(err);
  err = VecDestroy(soln);dCHK(err);
  err = SNESDestroy(snes);dCHK(err);
  if (J != Jp) {err = MatDestroy(J);dCHK(err);}
  err = MatDestroy(Jp);dCHK(err);
  err = ElastDestroy(elt);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
