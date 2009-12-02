static const char help[] = "Solve nonlinear elasticity using dual order hp elements.\n"
  "The model problem is\n"
  "  -div(lambda tr(Du)I + 2 mu (Du - gamma*dot(grad(U),grad(U))/2)) = f \n"
  "where\n"
  "  D is the symmetric gradient operator\n"
  "  mu and lambda are the Lame parameters\n"
  "  gamma controls the nonlinearity, gamma=1 is the standard finite strain tensor and gamma=0 is linear elasticity\n\n";

#include <petscsnes.h>
#include <dohpfs.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>

static PetscLogEvent LOG_ElastShellMult;

struct ElastParam {
  dReal lambda,mu,gamma;
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


struct ElastStore {
  dReal Du[3][3];
};

typedef struct ElastCtx *Elast;
struct ElastCtx {
  MPI_Comm              comm;
  struct ElastParam     param;
  struct ElastExact     exact;
  struct ElastExactCtx  exactctx;
  struct ElastStore    *store;
  dInt                 *storeoff;
  dJacobi               jac;
  dMesh                 mesh;
  dFS                   fs;
  Vec                   x,y;
  dInt                  constBDeg,nominalRDeg;
  dTruth                errorview;
};

static dErr ElastCreate(MPI_Comm comm,Elast *elast)
{
  Elast elt;
  dErr err;

  dFunctionBegin;
  *elast = 0;
  err = dNew(struct ElastCtx,&elt);dCHK(err);
  elt->comm = comm;

  elt->constBDeg = 4;
  elt->nominalRDeg = 0;
  elt->param.lambda = 10;
  elt->param.mu     = 100;
  elt->param.gamma  = 0;
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
  dMeshTag rtag,dtag;
  dInt exact;
  dErr err;

  dFunctionBegin;
  exact = 0; exc->a = exc->b = exc->c = 1;
  err = PetscOptionsBegin(elt->comm,NULL,"Elasticity options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-const_bdeg","Use constant isotropic degree on all elements","",elt->constBDeg,&elt->constBDeg,NULL);dCHK(err);
    err = PetscOptionsInt("-nominal_rdeg","Nominal rule degree (will be larger if basis requires it)","",elt->nominalRDeg,&elt->nominalRDeg,NULL);dCHK(err);
    err = PetscOptionsTruth("-error_view","View errors","",elt->errorview,&elt->errorview,NULL);dCHK(err);
    err = PetscOptionsReal("-elast_lambda","first Lame parameter","",prm->lambda,&prm->lambda,NULL);dCHK(err);
    err = PetscOptionsReal("-elast_mu","Second Lame parameter","",prm->mu,&prm->mu,NULL);dCHK(err);
    err = PetscOptionsReal("-elast_gamma","Strength of nonlinearity [0,1]","",prm->gamma,&prm->gamma,NULL);dCHK(err);
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
    default: dERROR(1,"Exact solution %d not implemented");
  }

  err = dMeshCreate(elt->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);dCHK(err);
  elt->mesh = mesh;
  err = dMeshGetRoot(mesh,&domain);dCHK(err); /* Need a taggable set */
  err = dMeshSetDuplicateEntsOnly(mesh,domain,&domain);dCHK(err);

  err = dJacobiCreate(elt->comm,&jac);dCHK(err);
  err = dJacobiSetDegrees(jac,9,2);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);
  elt->jac = jac;

  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"elast_rule_degree",elt->nominalRDeg,&rtag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"elast_efs_degree",elt->constBDeg,&dtag);dCHK(err);

  err = dFSCreate(elt->comm,&fs);dCHK(err);
  err = dFSSetBlockSize(fs,3);dCHK(err);
  err = dFSSetMesh(fs,mesh,domain);dCHK(err);
  err = dFSSetRuleTag(fs,jac,rtag);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  err = dFSRegisterBoundary(fs,100,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSRegisterBoundary(fs,200,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSRegisterBoundary(fs,300,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSSetFromOptions(fs);dCHK(err);
  elt->fs = fs;

  err = dFSCreateExpandedVector(fs,&elt->x);dCHK(err);
  err = VecDuplicate(elt->x,&elt->y);dCHK(err);

  {                             /* Allocate space for stored values */
    dInt n;
    s_dRule *rule;
    err = dFSGetElements(fs,&n,NULL,&rule,NULL,NULL,NULL);dCHK(err);
    err = dMallocA(n+1,&elt->storeoff);dCHK(err);
    elt->storeoff[0] = 0;
    for (dInt i=0; i<n; i++) {
      dInt q;
      err = dRuleGetSize(&rule[i],NULL,&q);dCHK(err);
      elt->storeoff[i+1] = elt->storeoff[i] + q;
    }
    err = dMallocA(elt->storeoff[n],&elt->store);dCHK(err);
    err = dMemzero(elt->store,elt->storeoff[n]*sizeof(elt->store[0]));dCHK(err);
    err = dFSRestoreElements(fs,&n,NULL,&rule,NULL,NULL,NULL);dCHK(err);
  }

  dFunctionReturn(0);
}

static dErr ElastDestroy(Elast elt)
{
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(elt->fs);dCHK(err);
  err = dJacobiDestroy(elt->jac);dCHK(err);
  err = dMeshDestroy(elt->mesh);dCHK(err);
  err = dFree(elt->storeoff);dCHK(err);
  err = dFree(elt->store);dCHK(err);
  if (elt->x) {err = VecDestroy(elt->x);dCHK(err);}
  if (elt->y) {err = VecDestroy(elt->y);dCHK(err);}
  err = dFree(elt);dCHK(err);
  dFunctionReturn(0);
}

static inline void ElastPointwiseComputeStore(struct ElastParam dUNUSED *prm,const dReal dUNUSED x[3],const dScalar dUNUSED u[3],const dScalar Du[],struct ElastStore *st)
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
  dFS fs = elt->fs;
  dInt n,*off,*geomoff;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*restrict geom)[3],(*restrict q)[3],(*restrict jinv)[3][3],*restrict jw;
  dScalar *x,*y,*restrict u,*restrict v,*restrict du,*restrict dv;
  dErr err;

  dFunctionBegin;
  err = dFSGlobalToExpanded(fs,gx,elt->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(elt->x,&x);dCHK(err);
  err = VecGetArray(elt->y,&y);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err); /* \note \a off is in terms of \e nodes, not \e dofs */
  err = dFSGetWorkspace(fs,__func__,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  for (dInt e=0; e<n; e++) {
    dInt Q;
    err = dRuleComputeGeometry(&rule[e],(const dReal(*)[3])(geom+geomoff[e]),q,jinv,jw);dCHK(err);
    err = dRuleGetSize(&rule[e],0,&Q);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,x+3*off[e],u,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,x+3*off[e],du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      struct ElastStore *restrict st = &elt->store[elt->storeoff[e]+i];
      ElastPointwiseFunction(&elt->param,&elt->exact,&elt->exactctx,q[i],jw[i],&u[i*3],&du[i*9],st,&v[i*3],&dv[i*9]);
    }
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,v,y+3*off[e],dAPPLY_INTERP_TRANSPOSE,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,dv,y+3*off[e],dAPPLY_GRAD_TRANSPOSE,ADD_VALUES);dCHK(err);
  }
  err = dFSRestoreWorkspace(fs,__func__,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(elt->x,&x);dCHK(err);
  err = VecRestoreArray(elt->y,&y);dCHK(err);
  err = dFSExpandedToGlobal(fs,elt->y,gy,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  dFunctionReturn(0);
}

static dErr ElastShellMatMult(Mat J,Vec gx,Vec gy)
{
  Elast elt;
  dFS fs;
  dInt n,*off,*geomoff;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*restrict geom)[3],(*restrict q)[3],(*restrict jinv)[3][3],*restrict jw;
  dScalar *x,*y,*restrict u,*restrict v,*restrict du,*restrict dv;
  dErr err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_ElastShellMult,J,gx,gy,0);dCHK(err);
  err = MatShellGetContext(J,(void**)&elt);dCHK(err);
  fs = elt->fs;
  err = dFSGlobalToExpanded(fs,gx,elt->x,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(elt->x,&x);dCHK(err);
  err = VecGetArray(elt->y,&y);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  for (dInt e=0; e<n; e++) {
    dInt Q;
    err = dRuleComputeGeometry(&rule[e],(const dReal(*)[3])(geom+geomoff[e]),q,jinv,jw);dCHK(err);
    err = dRuleGetSize(&rule[e],0,&Q);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,x+3*off[e],u,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,x+3*off[e],du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      struct ElastStore *restrict st = &elt->store[elt->storeoff[e]+i];
      ElastPointwiseJacobian(&elt->param,st,jw[i],&u[i*3],&du[i*9],&v[i*3],&dv[i*9]);
    }
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,v,y+3*off[e],dAPPLY_INTERP_TRANSPOSE,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,dv,y+3*off[e],dAPPLY_GRAD_TRANSPOSE,ADD_VALUES);dCHK(err);
  }
  err = dFSRestoreWorkspace(fs,__func__,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(elt->x,&x);dCHK(err);
  err = VecRestoreArray(elt->y,&y);dCHK(err);
  err = dFSExpandedToGlobal(fs,elt->y,gy,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = PetscLogEventEnd(LOG_ElastShellMult,J,gx,gy,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr ElastJacobian(SNES dUNUSED snes,Vec gx,Mat *J,Mat *Jp,MatStructure *structure,void *ctx)
{
  Elast elt = ctx;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*nx)[3];
  dScalar *x;
  dFS fs = elt->fs;
  dInt n,*off,*geomoff;
  dReal (*geom)[3];
  dErr err;

  dFunctionBegin;
  err = MatZeroEntries(*Jp);dCHK(err);
  err = dFSGlobalToExpanded(fs,gx,elt->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(elt->x,&x);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&nx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err); /* We only need space for nodal coordinates */
  for (dInt e=0; e<n; e++) {
    dInt three,P[3];
    err = dEFSGetGlobalCoordinates(&efs[e],(const dReal(*)[3])(geom+geomoff[e]),&three,P,nx);dCHK(err);
    if (three != 3) dERROR(1,"Dimension not equal to 3");
    for (dInt i=0; i<P[0]-1; i++) { /* P-1 = number of sub-elements in each direction */
      for (dInt j=0; j<P[1]-1; j++) {
        for (dInt k=0; k<P[2]-1; k++) {
          dQ1CORNER_CONST_DECLARE(c,rowcol,corners,off[e],nx,P,i,j,k);
          const dScalar (*uc)[3] = (const dScalar(*)[3])x+off[e]; /* function values, indexed at subelement corners \c uc[c[#]][0] */
          const dReal (*qx)[3],*jw,(*basis)[8],(*deriv)[8][3];
          dInt qn;
          dScalar K[8*3][8*3];
          err = dMemzero(K,sizeof(K));dCHK(err);
          err = dQ1HexComputeQuadrature(corners,&qn,&qx,&jw,(const dReal**)&basis,(const dReal**)&deriv);dCHK(err);
          for (dInt lq=0; lq<qn; lq++) { /* loop over quadrature points */
            struct ElastStore st;
            { /* Set up store */
              dReal st_u[3] = {0,0,0},st_Du[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
              for (dInt lp=0; lp<8; lp++) { /* Evaluate function values and gradients at this quadrature point */
                for (dInt f=0; f<3; f++) {  /* for each field */
                  st_u[f] += basis[lq][lp] * uc[c[lp]][f];
                  st_Du[f][0] += deriv[lq][lp][0] * uc[c[lp]][f];
                  st_Du[f][1] += deriv[lq][lp][1] * uc[c[lp]][f];
                  st_Du[f][2] += deriv[lq][lp][2] * uc[c[lp]][f];
                }
              }
              ElastPointwiseComputeStore(&elt->param,qx[lq],st_u,&st_Du[0][0],&st);
            }
            for (dInt ltest=0; ltest<8; ltest++) {              /* Loop over test basis functions (corners) */
              for (dInt lp=0; lp<8; lp++) {                     /* Loop over trial basis functions (corners) */
                for (dInt fp=0; fp<3; fp++) {                   /* Each field component of trial function */
                  dReal u[3] = {0,0,0},Du[3][3] = {{0,0,0},{0,0,0},{0,0,0}},v[3],Dv[3][3];
                  u[fp] = basis[lq][lp]; /* Trial function for only this field component */
                  Du[fp][0] = deriv[lq][lp][0];
                  Du[fp][1] = deriv[lq][lp][1];
                  Du[fp][2] = deriv[lq][lp][2];
                  /* Get the coefficients of test functions for each field component */
                  ElastPointwiseJacobian(&elt->param,&st,jw[lq],u,&Du[0][0],v,&Dv[0][0]);
                  for (dInt ftest=0; ftest<3; ftest++) { /* Insert contribution from each test function field component */
                    K[ltest*3+ftest][lp*3+fp] += //basis[lq][ltest] * v[0]
                      + deriv[lq][ltest][0] * Dv[ftest][0]
                      + deriv[lq][ltest][1] * Dv[ftest][1]
                      + deriv[lq][ltest][2] * Dv[ftest][2];
                  }
                }
              }
            }
          }
          err = dFSMatSetValuesBlockedExpanded(fs,*Jp,8,rowcol,8,rowcol,&K[0][0],ADD_VALUES);dCHK(err);
        }
      }
    }
  }
  err = dFSRestoreWorkspace(fs,__func__,&nx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(elt->x,&x);dCHK(err);

  err = MatAssemblyBegin(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  *structure = SAME_NONZERO_PATTERN;
  dFunctionReturn(0);
}

static dErr ElastErrorNorms(Elast elt,Vec gx,dReal errorNorms[static 3],dReal gerrorNorms[static 3])
{
  dFS fs = elt->fs;
  dInt n,*off,*geomoff;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*geom)[3],(*q)[3],(*jinv)[3][3],*jw;
  dScalar *x,(*u)[3],(*du)[9];
  dErr err;

  dFunctionBegin;
  err = dMemzero(errorNorms,3*sizeof(errorNorms));dCHK(err);
  err = dMemzero(gerrorNorms,3*sizeof(gerrorNorms));dCHK(err);
  err = dFSGlobalToExpanded(fs,gx,elt->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(elt->x,&x);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&q,&jinv,&jw,(dReal**)&u,NULL,(dReal**)&du,NULL);dCHK(err);
  for (dInt e=0; e<n; e++) {
    dInt Q;
    err = dRuleComputeGeometry(&rule[e],(const dReal(*)[3])(geom+geomoff[e]),q,jinv,jw);dCHK(err);
    err = dRuleGetSize(&rule[e],0,&Q);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,x+3*off[e],&u[0][0],dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,x+3*off[e],&du[0][0],dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[3],duu[9],r[3],gr[3],rsum=0,grsum=0;
      elt->exact.solution(&elt->exactctx,&elt->param,q[i],uu,duu);
      for (dInt j=0; j<3; j++) {
        r[j] = u[i][j] - uu[j]; /* Function error at point */
        rsum += dSqr(r[j]);
        gr[j] = dSqrt(dSqr(du[i][j*3+0]-duu[j*3+0]) + dSqr(du[i][j*3+1]-duu[j*3+1]) + dSqr(du[i][j*3+2]-duu[j*3+2])); /* Gradient error at point */
        grsum += dSqr(gr[j]);
      }
      if (elt->errorview) {
        printf("e,q = %3d %3d (% 5f,% 5f,% 5f) dohp %10.2e %10.2e %10.2e   exact %10.2e %10.2e %10.2e   error %10.e\n",
               e,i,q[i][0],q[i][1],q[i][2],u[i][0],u[i][1],u[i][2],uu[0],uu[1],uu[2],rsum);
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
  }
  err = dFSRestoreWorkspace(fs,__func__,&q,&jinv,&jw,(dReal**)&u,NULL,NULL,NULL);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(elt->x,&x);dCHK(err);
  errorNorms[1] = dSqrt(errorNorms[1]);
  gerrorNorms[1] = dSqrt(gerrorNorms[1]);
  dFunctionReturn(0);
}

static dErr ElastGetSolutionVector(Elast elt,Vec *insoln)
{
  dErr err;
  Vec soln,xc,cvec;
  dScalar *x,*coords;
  dInt n,bs;

  dFunctionBegin;
  *insoln = 0;
  err = dFSCreateGlobalVector(elt->fs,&soln);dCHK(err);
  err = VecDohpGetClosure(soln,&xc);dCHK(err);
  err = dFSGetCoordinates(elt->fs,&cvec);dCHK(err);
  err = VecGetLocalSize(xc,&n);dCHK(err);
  err = VecGetBlockSize(xc,&bs);dCHK(err);
  {
    dInt nc;
    err = VecGetLocalSize(cvec,&nc);dCHK(err);
    if (nc*bs != n*3) dERROR(1,"Coordinate vector has inconsistent size");
  }
  err = VecGetArray(xc,&x);dCHK(err);
  err = VecGetArray(cvec,&coords);dCHK(err);
  for (dInt i=0; i<n/bs; i++) {
    dScalar du_unused[3*bs];
    elt->exact.solution(&elt->exactctx,&elt->param,&coords[3*i],&x[i*bs],du_unused);
    /* printf("Node %3d: coords %+8f %+8f %+8f   exact %+8f %+8f %+8f\n",i,coords[3*i],coords[3*i+1],coords[3*i+2],x[3*i],x[3*i+1],x[3*i+2]); */
  }
  err = VecRestoreArray(xc,&x);dCHK(err);
  err = VecRestoreArray(cvec,&coords);dCHK(err);
  err = VecDestroy(cvec);dCHK(err);
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
  dTruth nojshell,nocheck;
  dErr err;

  err = dInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = PetscLogEventRegister("ElastShellMult",MAT_COOKIE,&LOG_ElastShellMult);dCHK(err);

  err = ElastCreate(comm,&elt);dCHK(err);
  err = ElastSetFromOptions(elt);dCHK(err);
  fs = elt->fs;

  err = dFSCreateGlobalVector(fs,&r);dCHK(err);
  err = PetscOptionsGetString(NULL,"-q1mat_type",mtype,sizeof(mtype),NULL);dCHK(err);
  err = dFSGetMatrix(fs,mtype,&Jp);dCHK(err);
  err = MatSetOptionsPrefix(Jp,"q1");dCHK(err);

  err = PetscOptionsBegin(elt->comm,NULL,"Elasticity solver options",__FILE__);dCHK(err); {
    err = PetscOptionsName("-nojshell","Do not use shell Jacobian","",&nojshell);dCHK(err);
    err = PetscOptionsName("-nocheck_error","Do not compute errors","",&nocheck);dCHK(err);
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
