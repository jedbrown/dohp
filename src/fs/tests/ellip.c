static const char help[] = "Solve a scalar elliptic problem, a regularized p-Bratu using dual order hp elements.\n"
  "The model problem is\n"
  "  \\int_\\Omega (\\eta Dv \\cdot Du - \\lambda e^u - f v) - \\int_\\Gamma v (\\eta Du \\cdot n) = 0\n"
  "where\n"
  "  \\eta(u) = (\\epsilon + 1/2 Du . Du)^((p-2)/2)\n"
  "  (\\eta Du \\cdot n) = known OR function of u OR self (\"No condition\" outflow)\n\n";

#include <dohp.h>
#include <dohpfs.h>
#include <dohpvec.h>
#include <dohpviewer.h>
#include <dohpsys.h>
#include <petscsnes.h>

PetscLogEvent LOG_EllipShellMatMult;

typedef struct EllipCtx *Ellip;

typedef struct {
  Ellip elp;
  Vec Mdiag,work0,work1;
  Mat Mq1;
  KSP ksp;
} PC_Ellip;

struct EllipParam {
  dReal epsilon;
  dReal p;
  dReal lambda;
  dBool  bdy100;
};

struct EllipExactCtx {
  dReal a,b,c;
};
struct EllipExact {
  void (*solution)(const struct EllipExactCtx*,const struct EllipParam*,const dReal x[3],dScalar u[1],dScalar du[3]);
  void (*forcing)(const struct EllipExactCtx*,const struct EllipParam*,const dReal x[3],dScalar f[1]);
};

/* Exact solutions are generated using sympy:
from sympy import *
from sympy.abc import *
u = cos(a*x)*exp(b*y)*sin(c*z)                                                                          # exact solution
ux,uy,uz=diff(u,x),diff(u,y),diff(u,z);                                                                 # convenience
gam = (ux**2+uy**2+uz**2)/2; eta = (e**2+gam)**((p-2)/2); f = -diff(eta*ux,x)-diff(eta*uy,y)-diff(eta*uz,z) # Physics
ccode(f)
*/
static void EllipExact_0_Solution(const struct EllipExactCtx *ctx,const struct EllipParam dUNUSED *prm,const dReal xyz[3],dScalar u[1],dScalar du[3])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2];
  u[0] = cos(a*x) * exp(b*y) * sin(c*z);
  du[0] = -a*sin(a*x) * exp(b*y) * sin(c*z);
  du[1] = cos(a*x) * b*exp(b*y) * sin(c*z);
  du[2] = cos(a*x) * exp(b*y) * cos(c*z);
}
static void EllipExact_0_Forcing(const struct EllipExactCtx *ctx,const struct EllipParam *prm,const dReal xyz[3],dScalar f[1])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2],e = prm->epsilon,p = prm->p;
  f[0] = pow(a,2)*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-1 + p/2))*cos(a*x)*exp(b*y)*sin(c*z) + pow(c,2)*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-1 + p/2))*cos(a*x)*exp(b*y)*sin(c*z) - pow(b,2)*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-1 + p/2))*cos(a*x)*exp(b*y)*sin(c*z) + b*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-2 + p/2))*(1 - p/2)*(pow(b,3)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y) + b*pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y) + b*pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y))*cos(a*x)*exp(b*y)*sin(c*z) - a*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-2 + p/2))*(1 - p/2)*(pow(a,3)*pow(sin(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) - a*pow(b,2)*pow(sin(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) - a*pow(c,2)*pow(cos(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x))*exp(b*y)*sin(a*x)*sin(c*z) + c*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-2 + p/2))*(1 - p/2)*(-pow(c,3)*pow(cos(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) + c*pow(a,2)*pow(sin(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) + c*pow(b,2)*pow(cos(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z))*cos(a*x)*cos(c*z)*exp(b*y);
}

static void EllipExact_1_Solution(const struct EllipExactCtx *ctx,const struct EllipParam dUNUSED *prm,const dReal xyz[3],dScalar u[1],dScalar du[3])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2];
  u[0] = a*dSqr(x) + b*dSqr(y) + c*(1-dSqr(z));
  du[0] = 2*a*x;
  du[1] = 2*b*y;
  du[2] = -2*c*z;
}
static void EllipExact_1_Forcing(const struct EllipExactCtx *ctx,const struct EllipParam *prm,const dReal xyz[3],dScalar f[1])
{
  const dUNUSED dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2],e = prm->epsilon,p = prm->p;
  f[0] = -8*pow(c,3)*pow(z,2)*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-2 + p/2))*(1 - p/2) + 8*pow(a,3)*pow(x,2)*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-2 + p/2))*(1 - p/2) + 8*pow(b,3)*pow(y,2)*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-2 + p/2))*(1 - p/2) - 2*a*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-1 + p/2)) - 2*b*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-1 + p/2)) + 2*c*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-1 + p/2));
}

static void EllipExact_2_Solution(const struct EllipExactCtx *ctx,const struct EllipParam dUNUSED *prm,const dReal xyz[3],dScalar u[1],dScalar du[3])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2];
  u[0] = a*x + b*y + c*z;
  du[0] = a;
  du[1] = b;
  du[2] = c;
}
static void EllipExact_2_Forcing(const struct EllipExactCtx dUNUSED *ctx,const struct EllipParam dUNUSED *prm,const dReal dUNUSED xyz[3],dScalar f[1])
{
  f[0] = 0;
}

static void EllipExact_3_Solution(const struct EllipExactCtx *ctx,const struct EllipParam dUNUSED *prm,const dReal xyz[3],dScalar u[1],dScalar du[3])
{
  const dReal x = xyz[0],y = xyz[1],z = xyz[2],a = ctx->a,b = ctx->b,c = ctx->c;
  u[0] = a*x*x*x + b*y*y*z + c*(1-z*z)*x;
  du[0] = 3*a*x*x + c*(1-z*z);
  du[1] = 2*b*y*z;
  du[2] = b*y*y - 2*c*z*x;
}
static void EllipExact_3_Forcing(const struct EllipExactCtx *ctx,const struct EllipParam *prm,const dReal xyz[3],dScalar f[1])
{
  const dReal x = xyz[0],y = xyz[1],z = xyz[2],a = ctx->a,b = ctx->b,c = ctx->c,e = prm->epsilon,p = prm->p;
  f[0] = pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-2 + p/2))*(1 - p/2)*(-2*c*x*z + b*pow(y,2))*(-2*c*x*(-2*c*x*z + b*pow(y,2)) - 2*c*z*(c*(1 - pow(z,2)) + 3*a*pow(x,2)) + 4*z*pow(b,2)*pow(y,2)) + pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-2 + p/2))*(1 - p/2)*(c*(1 - pow(z,2)) + 3*a*pow(x,2))*(-2*c*z*(-2*c*x*z + b*pow(y,2)) + 6*a*x*(c*(1 - pow(z,2)) + 3*a*pow(x,2))) + 2*b*y*z*pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-2 + p/2))*(1 - p/2)*(2*b*y*(-2*c*x*z + b*pow(y,2)) + 4*y*pow(b,2)*pow(z,2)) - 6*a*x*pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-1 + p/2)) - 2*b*z*pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-1 + p/2)) + 2*c*x*pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-1 + p/2));
}

struct EllipStore {
  dReal eta;
  dReal sqrt_mdeta_Du[3];
  dReal lambda_exp_u;
};

typedef enum {EVAL_FUNCTION,EVAL_JACOBIAN, EVAL_UB} EllipEvaluation;

struct EllipCtx {
  MPI_Comm              comm;
  struct EllipParam     param;
  struct EllipExact     exact;
  struct EllipExactCtx  exactctx;
  dJacobi               jac;
  dMesh                 mesh;
  dFS                   fs;
  dInt                  constBDeg;
  dBool                 errorview;
  dBool                 eta_monitor;
  dQuadratureMethod     function_qmethod,jacobian_qmethod;
  dRulesetIterator      regioniter[EVAL_UB];
};

static dErr EllipCreate(MPI_Comm comm,Ellip *ellip)
{
  struct EllipParam *prm;
  Ellip elp;
  dErr err;

  dFunctionBegin;
  *ellip = 0;
  err = dNew(struct EllipCtx,&elp);dCHK(err);
  elp->comm = comm;

  elp->constBDeg = 3;

  prm = &elp->param;
  prm->p           = 2.0;       /* p in p-Laplacian */
  prm->epsilon     = 1.0;
  prm->lambda      = 0.0;       /* Bratu nonlinearity */
  elp->function_qmethod = dQUADRATURE_METHOD_FAST;
  elp->jacobian_qmethod = dQUADRATURE_METHOD_SPARSE;

  *ellip = elp;
  dFunctionReturn(0);
}

typedef struct {
  dReal morph,twist,stretch;
} MorphCtx;

static void Morph(void *vctx,double *coords)
{
  MorphCtx *ctx = vctx;
  dReal morph = ctx->morph,twist = ctx->twist,stretch = ctx->stretch;
  double x = coords[0],y = coords[1],z = coords[2];
  coords[0] =  cos(twist*z)*x + sin(twist*z)*y;
  coords[1] = -sin(twist*z)*x + cos(twist*z)*y;
  coords[2] = -1 + (z + 1) * (1 + morph * (dSqr(x + 1) + dSqr(y + 1))) + stretch*z;
}

static dErr EllipSetFromOptions(Ellip elp)
{
  char mesh_out_name[256] = "ellip.h5m";
  struct EllipParam *prm = &elp->param;
  struct EllipExactCtx *exc = &elp->exactctx;
  dMesh mesh;
  dFS fs;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag dtag;
  dBool  mesh_out;
  dReal morph,twist,stretch;
  dInt exact;
  dErr err;

  dFunctionBegin;
  exact = 0; morph = twist = stretch = 0.0; mesh_out = dFALSE; exc->a = exc->b = exc->c = 1;
  err = PetscOptionsBegin(elp->comm,NULL,"Elliptic (p-Laplacian) options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-const_bdeg","Use constant isotropic degree on all elements","",elp->constBDeg,&elp->constBDeg,NULL);dCHK(err);
    err = PetscOptionsBool("-error_view","View errors","",elp->errorview,&elp->errorview,NULL);dCHK(err);
    err = PetscOptionsBool("-eta_monitor","Monitor nonlinearity","",elp->eta_monitor,&elp->eta_monitor,NULL);dCHK(err);
    err = PetscOptionsReal("-ellip_p","p in p-Laplacian","",prm->p,&prm->p,NULL);dCHK(err);
    err = PetscOptionsReal("-ellip_eps","Regularization in p-Laplacian","",prm->epsilon,&prm->epsilon,NULL);dCHK(err);
    err = PetscOptionsReal("-ellip_lam","Strength of Bratu nonlinearity","",prm->lambda,&prm->lambda,NULL);dCHK(err);
    err = PetscOptionsEnum("-ellip_f_qmethod","Quadrature method for residual evaluation/matrix-free","",dQuadratureMethods,(PetscEnum)elp->function_qmethod,(PetscEnum*)&elp->function_qmethod,NULL);dCHK(err);
    err = PetscOptionsEnum("-ellip_jac_qmethod","Quadrature to use for Jacobian assembly","",dQuadratureMethods,(PetscEnum)elp->jacobian_qmethod,(PetscEnum*)&elp->jacobian_qmethod,NULL);dCHK(err);
    err = PetscOptionsBool("-bdy100","Only use boundary 100","",prm->bdy100,&prm->bdy100,NULL);dCHK(err);
    err = PetscOptionsInt("-exact","Exact solution choice","",exact,&exact,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_a","First scale parameter","",exc->a,&exc->a,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_b","Second scale parameter","",exc->b,&exc->b,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_c","Third scale parameter","",exc->c,&exc->c,NULL);dCHK(err);
    err = PetscOptionsReal("-morph","Deform the mesh by this factor","",morph,&morph,NULL);dCHK(err);
    err = PetscOptionsReal("-twist","Twist the mesh by this factor","",twist,&twist,NULL);dCHK(err);
    err = PetscOptionsReal("-stretch","Stretch the mesh by this factor","",stretch,&stretch,NULL);dCHK(err);
    err = PetscOptionsString("-mesh_out","Write the (morphed) mesh with this name","",mesh_out_name,mesh_out_name,sizeof mesh_out_name,&mesh_out);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);

  switch (exact) {
    case 0:
      elp->exact.solution = EllipExact_0_Solution;
      elp->exact.forcing = EllipExact_0_Forcing;
      break;
    case 1:
      elp->exact.solution = EllipExact_1_Solution;
      elp->exact.forcing = EllipExact_1_Forcing;
      break;
    case 2:
      elp->exact.solution = EllipExact_2_Solution;
      elp->exact.forcing = EllipExact_2_Forcing;
      break;
    case 3:
      elp->exact.solution = EllipExact_3_Solution;
      elp->exact.forcing = EllipExact_3_Forcing;
      break;
    default: dERROR(PETSC_COMM_SELF,1,"Exact solution %d not implemented");
  }

  err = dMeshCreate(elp->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);dCHK(err);
  {
    MorphCtx ctx = {.morph   = morph,
                    .twist   = twist,
                    .stretch = stretch};
    err = dMeshMorph(mesh,Morph,&ctx);dCHK(err);
  }
  elp->mesh = mesh;
  err = dMeshGetRoot(mesh,&domain);dCHK(err); /* Need a taggable set */
  err = dMeshSetDuplicateEntsOnly(mesh,domain,&domain);dCHK(err);

  err = dJacobiCreate(elp->comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  elp->jac = jac;

  err = dMeshCreateRuleTagIsotropic(mesh,domain,"ellip_efs_degree",elp->constBDeg,&dtag);dCHK(err);

  err = dFSCreate(elp->comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,domain);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  err = dFSRegisterBoundary(fs,100,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  if (!elp->param.bdy100) {
    err = dFSRegisterBoundary(fs,200,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
    err = dFSRegisterBoundary(fs,300,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  }
  err = dFSSetFromOptions(fs);dCHK(err);
  err = PetscObjectSetName((PetscObject)fs,"dFS_0");dCHK(err);
  err = PetscObjectSetName((PetscObject)mesh,"dMesh_0");dCHK(err);
  elp->fs = fs;

  if (mesh_out) {
    iMesh_Instance mi;
    dIInt          ierr;
    err = dMeshGetInstance(mesh,&mi);dCHK(err);
    iMesh_save(mi,domain,mesh_out_name,"",&ierr,(int)strlen(mesh_out_name),0);dICHK(mi,ierr);
  }

  dFunctionReturn(0);
}

static dErr EllipGetRegionIterator(Ellip elp,EllipEvaluation eval,dRulesetIterator *riter)
{
  dErr err;

  dFunctionBegin;
  if (!elp->regioniter[eval]) {
    dRulesetIterator iter;
    dRuleset ruleset;
    dFS cfs;
    dMeshESH domain;
    dQuadratureMethod qmethod;
    switch (eval) {
    case EVAL_FUNCTION: qmethod = elp->function_qmethod; break;
    case EVAL_JACOBIAN: qmethod = elp->jacobian_qmethod; break;
    default: dERROR(elp->comm,PETSC_ERR_ARG_OUTOFRANGE,"Unknown evaluation context");
    }
    err = dFSGetDomain(elp->fs,&domain);dCHK(err);
    err = dFSGetPreferredQuadratureRuleSet(elp->fs,domain,dTYPE_REGION,dTOPO_ALL,qmethod,&ruleset);dCHK(err);
    err = dFSGetCoordinateFS(elp->fs,&cfs);dCHK(err);
    err = dRulesetCreateIterator(ruleset,cfs,&iter);dCHK(err);
    err = dRulesetDestroy(&ruleset);dCHK(err); /* Give ownership to iterator */
    err = dRulesetIteratorAddFS(iter,elp->fs);dCHK(err);
    if (eval == EVAL_FUNCTION) {err = dRulesetIteratorAddStash(iter,0,sizeof(struct EllipStore));dCHK(err);}
    elp->regioniter[eval] = iter;
  }
  *riter = elp->regioniter[eval];
  dFunctionReturn(0);
}

static dErr EllipDestroy(Ellip elp)
{
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(&elp->fs);dCHK(err);
  err = dJacobiDestroy(&elp->jac);dCHK(err);
  err = dMeshDestroy(&elp->mesh);dCHK(err);
  for (dInt i=0; i<EVAL_UB; i++) {
    err = dRulesetIteratorDestroy(&elp->regioniter[i]);dCHK(err);
  }
  err = dFree(elp);dCHK(err);
  dFunctionReturn(0);
}

static inline void EllipPointwiseComputeStore(struct EllipParam *prm,const dReal dUNUSED x[3],const dScalar dUNUSED u[1],const dScalar Du[3],struct EllipStore *st)
{
  dReal gamma,espg,p=prm->p,sqrt_mdeta;
  gamma = 0.5 * (dSqr(Du[0]) + dSqr(Du[1]) + dSqr(Du[2]));
  espg = dSqr(prm->epsilon) + gamma;
  st->eta = pow(espg,(p-2)/2);
  sqrt_mdeta = sqrt(-((p-2)/2) * st->eta / espg);
  for (dInt i=0; i<3; i++) st->sqrt_mdeta_Du[i] = sqrt_mdeta*Du[i];
  st->lambda_exp_u = prm->lambda * exp(u[0]);
}

static inline void EllipPointwiseFunction(struct EllipParam *prm,struct EllipExact *exact,struct EllipExactCtx *exactctx,
                                          const dReal x[3],dReal weight,const dScalar u[1],const dScalar Du[3],
                                          struct EllipStore *st,dScalar v[1],dScalar Dv[3])
{
  dScalar f[1];
  EllipPointwiseComputeStore(prm,x,u,Du,st);
  exact->forcing(exactctx,prm,x,f);
  v[0] = - weight * f[0];       /* Coefficient of \a v in weak form */
  v[0] += -weight * st->lambda_exp_u;
  for (dInt i=0; i<3; i++) {
    Dv[i] = weight * st->eta * Du[i]; /* Coefficient of Dv in weak form */
  }
}

static inline void EllipPointwiseJacobian(struct EllipParam dUNUSED *prm,const struct EllipStore *restrict st,dReal weight,
                                          const dScalar u[restrict static 1],const dScalar Du[restrict static 3],
                                          dScalar v[restrict static 1],dScalar Dv[restrict static 3])
{
  const dScalar dotw = dDotScalar3(st->sqrt_mdeta_Du,Du)*weight;
  const dReal etaw = st->eta*weight;
  v[0] = -weight * st->lambda_exp_u * u[0];
  for (dInt i=0; i<3; i++) Dv[i] = etaw*Du[i] - dotw*st->sqrt_mdeta_Du[i];
}

static dErr EllipFunction(SNES dUNUSED snes,Vec gx,Vec gy,void *ctx)
{
  Ellip elp = ctx;
  Vec Coords;
  dRulesetIterator iter;
  dReal mineta=1e10,maxeta=1e-10;
  dErr err;

  dFunctionBegin;
  err = VecZeroEntries(gy);dCHK(err);
  err = dFSGetGeometryVectorExpanded(elp->fs,&Coords);dCHK(err);
  err = EllipGetRegionIterator(elp,EVAL_FUNCTION,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,gy,dFS_INHOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar *x,*dx,*u,*du,*v,*dv;
    dInt Q;
    struct EllipStore *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,&v,&dv);dCHK(err);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      EllipPointwiseFunction(&elp->param,&elp->exact,&elp->exactctx,&x[i*3],jw[i],&u[i],&du[i*3],&stash[i],&v[i],&dv[i*3]);
      maxeta = dMax(maxeta,stash[i].eta);
      mineta = dMin(mineta,stash[i].eta);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES,NULL,NULL,v,dv);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  if (elp->eta_monitor) {
    err = dPrintf(elp->comm,"## Eta min %f, max %f, ratio %f\n",mineta,maxeta,maxeta/mineta);
  }
  dFunctionReturn(0);
}

static dErr EllipShellMatMult(Mat J,Vec gx,Vec gy)
{
  Ellip elp;
  Vec Coords;
  dRulesetIterator iter;
  dErr err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_EllipShellMatMult,J,0,0,0);dCHK(err);
  err = MatShellGetContext(J,(void**)&elp);dCHK(err);
  err = VecZeroEntries(gy);dCHK(err);
  err = dFSGetGeometryVectorExpanded(elp->fs,&Coords);dCHK(err);
  err = EllipGetRegionIterator(elp,EVAL_FUNCTION,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_HOMOGENEOUS,gy,dFS_HOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar *x,*dx,*u,*du,*v,*dv;
    dInt Q;
    const struct EllipStore *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,&v,&dv);dCHK(err);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      EllipPointwiseJacobian(&elp->param,&stash[i],jw[i],&u[i],&du[i*3],&v[i],&dv[i*3]);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES,NULL,NULL,v,dv);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = PetscLogEventEnd(LOG_EllipShellMatMult,J,0,0,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr EllipJacobian(SNES snes,Vec gx,Mat *J,Mat *Jp,MatStructure *structure,void *ctx)
{
  Ellip            elp = ctx;
  KSP              ksp;
  PC               pc;
  PC_Ellip         *pce;
  Vec              Coords;
  dRulesetIterator iter;
  dErr             err;
  dScalar          *Kflat;

  dFunctionBegin;
  err = SNESGetKSP(snes,&ksp);dCHK(err);
  err = KSPGetPC(ksp,&pc);dCHK(err);
  err = PCShellGetContext(pc,(void**)&pce);dCHK(err);
  err = MatZeroEntries(*Jp);dCHK(err);
  err = dFSGetGeometryVectorExpanded(elp->fs,&Coords);dCHK(err);
  err = EllipGetRegionIterator(elp,EVAL_JACOBIAN,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSplit(iter,NULL,NULL,NULL,&Kflat);dCHK(err);

  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw,*interp_flat,*deriv_flat;
    const dInt *rowcol;
    dScalar *x,*dx,*u,*du;
    dInt Q,P;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,NULL,NULL);dCHK(err);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat);dCHK(err);
    {                           /* Scope so that we can declare new VLA pointers for convenient assembly */
      const dReal (*interp)[P] = (const dReal(*)[P])interp_flat;
      const dReal (*deriv)[P][3] = (const dReal(*)[P][3])deriv_flat;
      dScalar (*K)[P] = (dScalar(*)[P])Kflat;
      err = PetscMemzero(K,P*P*sizeof(K[0][0]));dCHK(err);
      for (dInt q=0; q<Q; q++) {
        struct EllipStore store;
        EllipPointwiseComputeStore(&elp->param,&x[3*q],&u[q],&du[q*3],&store);
        for (dInt j=0; j<P; j++) {
          dScalar v[1],dv[3];
          EllipPointwiseJacobian(&elp->param,&store,jw[q],&interp[q][j],deriv[q][j],v,dv);
          for (dInt i=0; i<P; i++) {
            K[i][j] += (interp[q][i] * v[0]
                        + deriv[q][i][0] * dv[0]
                        + deriv[q][i][1] * dv[1]
                        + deriv[q][i][2] * dv[2]);
          }
        }
      }
      err = dFSMatSetValuesBlockedExpanded(elp->fs,*Jp,P,rowcol,P,rowcol,&K[0][0],ADD_VALUES);dCHK(err);
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

static dErr EllipErrorNorms(Ellip elp,Vec gx,dReal errorNorms[static 3],dReal gerrorNorms[static 3])
{
  dErr err;
  dInt patchcnt;
  Vec Coords;
  dRulesetIterator iter;

  dFunctionBegin;
  err = dMemzero(errorNorms,3*sizeof(errorNorms));dCHK(err);
  err = dMemzero(gerrorNorms,3*sizeof(gerrorNorms));dCHK(err);
  err = dFSGetGeometryVectorExpanded(elp->fs,&Coords);dCHK(err);
  err = EllipGetRegionIterator(elp,EVAL_FUNCTION,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  patchcnt = 0;
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[3][3],(*u)[1],(*du)[1][3];
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,NULL);dCHK(err);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[1],duu[1][3],r[1],gr[3];             /* Scalar problem */
      dReal grsum;
      elp->exact.solution(&elp->exactctx,&elp->param,x[i],uu,(dScalar*)duu);
      r[0] = u[i][0] - uu[0];   /* Function error at point */
      gr[0] = du[i][0][0] - duu[0][0]; /* Gradient error at point */
      gr[1] = du[i][0][1] - duu[0][1];
      gr[2] = du[i][0][2] - duu[0][2];
      if (elp->errorview) {
        printf("e,q = %3d %3d (% 5f,% 5f,% 5f) dohp %10.2e   exact %10.2e   error %10.e\n",patchcnt,i,x[i][0],x[i][1],x[i][2],u[i][0],uu[0],r[0]);
      }
      grsum = dAbs(dDotScalar3(gr,gr));
      errorNorms[0] += dAbs(r[0]) * jw[i];               /* 1-norm */
      errorNorms[1] += dSqr(r[0]) * jw[i];               /* 2-norm */
      errorNorms[2] = dMax(errorNorms[2],dAbs(r[0])); /* Sup-norm */
      gerrorNorms[0] += grsum * jw[i];
      gerrorNorms[1] += dSqr(grsum) * jw[i];
      gerrorNorms[2] = dMax(gerrorNorms[2],grsum);
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

static dErr EllipGetNodalSolutionVector(Ellip elp,Vec *insoln)
{
  dErr err;
  Vec soln,xc,cvec;
  dScalar *x;
  const dScalar *coords;
  dInt n,bs;

  dFunctionBegin;
  *insoln = 0;
  err = dFSCreateGlobalVector(elp->fs,&soln);dCHK(err);
  err = VecDohpGetClosure(soln,&xc);dCHK(err);
  err = dFSGetNodalCoordinatesGlobal(elp->fs,&cvec);dCHK(err);
  err = VecGetLocalSize(xc,&n);dCHK(err);
  err = VecGetBlockSize(xc,&bs);dCHK(err);
  err = VecGetArray(xc,&x);dCHK(err);
  err = VecGetArrayRead(cvec,&coords);dCHK(err);
  for (dInt i=0; i<n/bs; i++) {
    dScalar du_unused[3*bs];
    elp->exact.solution(&elp->exactctx,&elp->param,&coords[3*i],&x[i*bs],du_unused);
  }
  err = VecRestoreArray(xc,&x);dCHK(err);
  err = VecRestoreArrayRead(cvec,&coords);dCHK(err);
  err = VecDohpRestoreClosure(soln,&xc);dCHK(err);
  *insoln = soln;
  dFunctionReturn(0);
}

static dErr PCApply_Ellip(PC pc,Vec x,Vec y)
{
  dErr err;
  PC_Ellip *pce;

  dFunctionBegin;
  err = PCShellGetContext(pc,(void**)&pce);dCHK(err);
  if (1) {
    err = VecPointwiseDivide(pce->work0,x,pce->Mdiag);dCHK(err);
  } else {
    err = VecCopy(x,pce->work0);dCHK(err);
  }
  if (1) {
    err = MatMult(pce->Mq1,pce->work0,pce->work1);dCHK(err);
  } else {
    err = VecCopy(pce->work0,pce->work1);dCHK(err);
  }
  err = KSPSolve(pce->ksp,pce->work1,y);dCHK(err);
  dFunctionReturn(0);
}

static dErr PCSetUp_Ellip(PC pc)
{
  dErr      err;
  PC_Ellip *pce;
  Mat       pmat;
  MatStructure mstruct;

  dFunctionBegin;
  err = PCShellGetContext(pc,(void**)&pce);dCHK(err);
  if (!pce->Mdiag) dERROR(PETSC_COMM_SELF,1,"Mdiag has not been set");
  if (!pce->Mq1) dERROR(PETSC_COMM_SELF,1,"Mq1 has not been set");
  if (!pce->work0) {err = VecDuplicate(pce->Mdiag,&pce->work0);dCHK(err);}
  if (!pce->work1) {err = VecDuplicate(pce->Mdiag,&pce->work1);dCHK(err);}
  if (!pce->ksp) {
    err = KSPCreate(pce->elp->comm,&pce->ksp);dCHK(err);
    err = KSPSetOptionsPrefix(pce->ksp,"ellip_");dCHK(err);
    err = KSPSetFromOptions(pce->ksp);dCHK(err);
  }
  err = PCGetOperators(pc,NULL,&pmat,&mstruct);dCHK(err);
  err = KSPSetOperators(pce->ksp,pmat,pmat,mstruct);dCHK(err);
  dFunctionReturn(0);
}

static dErr PCDestroy_Ellip(PC pc)
{
  dErr err;
  PC_Ellip *pce;

  dFunctionBegin;
  err = PCShellGetContext(pc,(void**)&pce);dCHK(err);
  err = VecDestroy(&pce->Mdiag);dCHK(err);
  err = VecDestroy(&pce->work0);dCHK(err);
  err = VecDestroy(&pce->work1);dCHK(err);
  err = MatDestroy(&pce->Mq1);dCHK(err);
  err = KSPDestroy(&pce->ksp);dCHK(err);
  err = dFree(pce);dCHK(err);
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  char mtype[256] = MATAIJ;
  Ellip elp;
  dFS fs;
  MPI_Comm comm;
  Mat J,Jp;
  Vec r,x,soln;
  SNES snes;
  dBool  nojshell,nocheck,viewdhm;
  dErr err;

  err = dInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;

  err = PetscLogEventRegister("EllipShellMult",MAT_CLASSID,&LOG_EllipShellMatMult);dCHK(err);

  err = EllipCreate(comm,&elp);dCHK(err);
  err = EllipSetFromOptions(elp);dCHK(err);
  fs = elp->fs;

  err = dFSCreateGlobalVector(fs,&r);dCHK(err);
  err = PetscOptionsGetString(NULL,"-q1mat_type",mtype,sizeof(mtype),NULL);dCHK(err);
  err = dFSGetMatrix(fs,mtype,&Jp);dCHK(err);
  err = MatSetOptionsPrefix(Jp,"q1");dCHK(err);
  err = MatSetFromOptions(Jp);dCHK(err);

  nojshell = dFALSE; nocheck = dFALSE; viewdhm = dFALSE;
  err = PetscOptionsBegin(elp->comm,NULL,"Elliptic solver options",__FILE__);dCHK(err); {
    err = PetscOptionsBool("-nojshell","Do not use shell Jacobian","",nojshell,&nojshell,NULL);dCHK(err);
    err = PetscOptionsBool("-nocheck_error","Do not compute errors","",nocheck,&nocheck,NULL);dCHK(err);
    err = PetscOptionsBool("-viewdhm","View to a file using DHM","",viewdhm,&viewdhm,NULL);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);
  if (nojshell) {
    /* Use the preconditioning matrix in place of the Jacobian.  This will NOT converge unless the elements are actually
    * Q1 (bdeg=1).  This option is nullified by -snes_mf_operator which will still only use the assembled Jacobian for
    * preconditioning. */
    J = Jp;
  } else {
    dInt m,n,M,N;
    err = MatGetLocalSize(Jp,&m,&n);dCHK(err);
    err = MatGetSize(Jp,&M,&N);dCHK(err);
    err = MatCreateShell(comm,m,n,M,N,elp,&J);dCHK(err);
    err = MatSetOptionsPrefix(J,"j");dCHK(err);
    err = MatShellSetOperation(J,MATOP_MULT,(void(*)(void))EllipShellMatMult);dCHK(err);
  }
  err = SNESCreate(comm,&snes);dCHK(err);
  err = SNESSetFunction(snes,r,EllipFunction,elp);dCHK(err);
  err = SNESSetJacobian(snes,J,Jp,EllipJacobian,elp);dCHK(err);
  err = SNESSetTolerances(snes,PETSC_DEFAULT,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);dCHK(err);
  {                             /* Set default PC */
    KSP ksp;
    PC pc;
    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPGetPC(ksp,&pc);dCHK(err);
    err = PCSetType(pc,PCSHELL);dCHK(err);
  }
  err = SNESSetFromOptions(snes);dCHK(err);
  err = VecZeroEntries(r);dCHK(err);
  err = VecDuplicate(r,&x);dCHK(err);
  err = PetscObjectSetName((PetscObject)x,"U");dCHK(err);
  err = EllipGetNodalSolutionVector(elp,&soln);dCHK(err);
  {
    Vec sc;
    err = VecDohpGetClosure(soln,&sc);dCHK(err);
    err = dFSInhomogeneousDirichletCommit(elp->fs,sc);dCHK(err);
    err = VecDohpRestoreClosure(soln,&sc);dCHK(err);
  }
  {
    KSP ksp;
    PC pc;
    dBool  isshell;
    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPGetPC(ksp,&pc);dCHK(err);
    err = PetscTypeCompare((dObject)pc,PCSHELL,&isshell);dCHK(err);
    if (isshell) {
      PC_Ellip *pce;
      err = dNew(PC_Ellip,&pce);dCHK(err);
      pce->elp = elp;
      err = PetscOptionsGetString(NULL,"-mq1mat_type",mtype,sizeof(mtype),NULL);dCHK(err);
      err = dFSGetMatrix(fs,mtype,&pce->Mq1);dCHK(err);
      err = MatSetOptionsPrefix(pce->Mq1,"mq1");dCHK(err);
      err = MatSetFromOptions(pce->Mq1);dCHK(err);
      err = VecDuplicate(x,&pce->Mdiag);dCHK(err);
      err = PCShellSetContext(pc,pce);dCHK(err);
      err = PCShellSetApply(pc,PCApply_Ellip);dCHK(err);
      err = PCShellSetSetUp(pc,PCSetUp_Ellip);dCHK(err);
      err = PCShellSetDestroy(pc,PCDestroy_Ellip);dCHK(err);
    }
  }
  err = VecZeroEntries(x);dCHK(err);
  err = SNESSolve(snes,NULL,x);dCHK(err);
  if (!nocheck) {
    dInt n;
    dReal anorm[2],anorminf,inorm[3],enorm[3],gnorm[3];
    err = EllipErrorNorms(elp,x,enorm,gnorm);dCHK(err);
    err = VecNorm(r,NORM_1_AND_2,anorm);dCHK(err);
    err = VecNorm(r,NORM_INFINITY,&anorminf);dCHK(err);
    err = VecWAXPY(r,-1,soln,x);dCHK(err);
    err = VecNorm(r,NORM_1_AND_2,inorm);dCHK(err);
    err = VecNorm(r,NORM_INFINITY,&inorm[2]);dCHK(err);
    /* Scale 1-norms for constant domain size (to mimic L^p instead of l^p) */
    err = VecGetSize(r,&n);dCHK(err);
    anorm[0] /= 1.*n;
    inorm[0] /= 1.*n;
    enorm[0] /= 1.*n;
    gnorm[0] /= 1.*n;
    /* Correct 2-norms for domain size */
    anorm[1] /= sqrt(1.*n);
    inorm[1] /= sqrt(1.*n);
    enorm[1] /= sqrt(1.*n);
    gnorm[1] /= sqrt(1.*n);
    /* Limit anorm so it is not reported as so small that different rounding modes change the result */
    anorm[0] *= (anorm[0] > 1e-14);
    anorm[1] *= (anorm[1] > 1e-14);
    anorminf *= (anorminf > 1e-14);
    err = dPrintf(comm,"Algebraic residual        |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",anorm[0],anorm[1],anorminf);dCHK(err);
    err = dPrintf(comm,"Interpolation residual    |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",inorm[0],inorm[1],inorm[2]);dCHK(err);
    err = dPrintf(comm,"Pointwise solution error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",enorm[0],enorm[1],enorm[2]);dCHK(err);
    err = dPrintf(comm,"Pointwise gradient error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",gnorm[0],gnorm[1],gnorm[2]);dCHK(err);
  }

  if (viewdhm) {
    dViewer viewer;
    err = PetscViewerCreate(comm,&viewer);dCHK(err);
    err = PetscViewerSetType(viewer,PETSCVIEWERDHM);dCHK(err);
    err = PetscViewerFileSetName(viewer,"ellip.dhm");dCHK(err);
    err = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);dCHK(err);
    err = dViewerDHMSetTimeUnits(viewer,"hour",PETSC_PI*1e7/3600);dCHK(err);
    err = dViewerDHMSetTime(viewer,0.1);dCHK(err);
    err = VecView(x,viewer);dCHK(err);
    err = PetscViewerDestroy(&viewer);dCHK(err);
  }

  err = VecDestroy(&r);dCHK(err);
  err = VecDestroy(&x);dCHK(err);
  err = VecDestroy(&soln);dCHK(err);
  err = SNESDestroy(&snes);dCHK(err);
  if (J != Jp) {err = MatDestroy(&J);dCHK(err);}
  err = MatDestroy(&Jp);dCHK(err);
  err = EllipDestroy(elp);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
