static const char help[] = "Solve the for the first eigen-pair of the p-Laplacian.\n"
  "The model problem is\n"
  "  \\int_\\Omega (\\eta Dv \\cdot Du - f v) = 0\n"
  "where\n"
  "  \\eta(u) = (\\epsilon + 1/2 Du . Du)^((p-2)/2)\n\n";

#include <dohp.h>
#include <dohpfs.h>
#include <dohpvec.h>
#include <dohpviewer.h>
#include <dohpsys.h>
#include <petscsnes.h>
#include <petscpf.h>

#define ALEN(a) ((dInt)(sizeof(a)/sizeof(a)[0]))
PetscLogEvent LOG_EllipShellMatMult;

typedef struct EllipCtx *Ellip;

struct EllipParam {
  dReal epsilon;
  dReal p;
  dReal lambda;
  dInt  dirichlet[16];           /* Dirichlet sets, 0 means unused */
};

struct EllipExactCtx {
  dReal a,b,c;
};
struct EllipExact {
  void (*solution)(const struct EllipExactCtx*,const struct EllipParam*,const dReal x[3],dScalar u[1],dScalar du[3]);
  void (*forcing)(const struct EllipExactCtx*,const struct EllipParam*,const dReal x[3],dScalar f[1]);
};

static void EllipExact_0_Solution(const struct EllipExactCtx dUNUSED *ctx,const struct EllipParam dUNUSED *prm,const dReal dUNUSED xyz[3],dScalar u[1],dScalar du[3])
{ // This is not actually an exact solution
  u[0] = 0;
  du[0] = 0;
  du[1] = 0;
  du[2] = 0;
}
static void EllipExact_0_Forcing(const struct EllipExactCtx *ctx,const struct EllipParam dUNUSED *prm,const dReal dUNUSED xyz[3],dScalar f[1])
{
  f[0] = ctx->a;
}

/* Exact solutions are generated using sympy:
from sympy import *
from sympy.abc import *
u = cos(a*x)*exp(b*y)*sin(c*z)                                                                          # exact solution
ux,uy,uz=diff(u,x),diff(u,y),diff(u,z);                                                                 # convenience
gam = (ux**2+uy**2+uz**2)/2; eta = (e**2+gam)**((p-2)/2); f = -diff(eta*ux,x)-diff(eta*uy,y)-diff(eta*uz,z) # Physics
ccode(f)
*/
static void EllipExact_1_Solution(const struct EllipExactCtx *ctx,const struct EllipParam dUNUSED *prm,const dReal xyz[3],dScalar u[1],dScalar du[3])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2];
  u[0] = cos(a*x) * exp(b*y) * sin(c*z);
  du[0] = -a*sin(a*x) * exp(b*y) * sin(c*z);
  du[1] = cos(a*x) * b*exp(b*y) * sin(c*z);
  du[2] = cos(a*x) * exp(b*y) * c*cos(c*z);
}
static void EllipExact_1_Forcing(const struct EllipExactCtx *ctx,const struct EllipParam *prm,const dReal xyz[3],dScalar f[1])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2],e = prm->epsilon,p = prm->p;
  f[0] = pow(a,2)*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-1 + p/2))*cos(a*x)*exp(b*y)*sin(c*z) + pow(c,2)*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-1 + p/2))*cos(a*x)*exp(b*y)*sin(c*z) - pow(b,2)*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-1 + p/2))*cos(a*x)*exp(b*y)*sin(c*z) + b*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-2 + p/2))*(1 - p/2)*(pow(b,3)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y) + b*pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y) + b*pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y))*cos(a*x)*exp(b*y)*sin(c*z) - a*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-2 + p/2))*(1 - p/2)*(pow(a,3)*pow(sin(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) - a*pow(b,2)*pow(sin(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) - a*pow(c,2)*pow(cos(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x))*exp(b*y)*sin(a*x)*sin(c*z) + c*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-2 + p/2))*(1 - p/2)*(-pow(c,3)*pow(cos(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) + c*pow(a,2)*pow(sin(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) + c*pow(b,2)*pow(cos(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z))*cos(a*x)*cos(c*z)*exp(b*y);
}

struct EllipStore {
  dReal eta;
  dReal sqrt_mdeta_Du[3];
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
  Vec                   rhs;
  dQuadratureMethod     function_qmethod,jacobian_qmethod;
  dRulesetIterator      regioniter[EVAL_UB];
};

struct PLapEig {
  dReal mu,q;
};
static dErr PFApply_PLapEig(void *ctx,PetscInt n,const PetscScalar *x,PetscScalar *y)
{
  struct PLapEig *pleig = ctx;
  dReal mu = pleig->mu,q = pleig->q;
  for (dInt i=0; i<n; i++) y[i] = mu * pow(x[i],q-1);
  return 0;
}

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
  prm->dirichlet[0] = 100;
  prm->dirichlet[1] = 200;
  prm->dirichlet[2] = 300;
  elp->function_qmethod = dQUADRATURE_METHOD_FAST;
  elp->jacobian_qmethod = dQUADRATURE_METHOD_SPARSE;

  *ellip = elp;
  dFunctionReturn(0);
}

static dErr EllipSetFromOptions(Ellip elp)
{
  struct EllipParam *prm = &elp->param;
  struct EllipExactCtx *exc = &elp->exactctx;
  dMesh mesh;
  dFS fs;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag dtag;
  dInt exact;
  dErr err;

  dFunctionBegin;
  exact = 0; exc->a = exc->b = exc->c = 1;
  err = PetscOptionsBegin(elp->comm,NULL,"Elliptic (p-Laplacian) options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-const_bdeg","Use constant isotropic degree on all elements","",elp->constBDeg,&elp->constBDeg,NULL);dCHK(err);
    err = PetscOptionsBool("-error_view","View errors","",elp->errorview,&elp->errorview,NULL);dCHK(err);
    err = PetscOptionsBool("-eta_monitor","Monitor nonlinearity","",elp->eta_monitor,&elp->eta_monitor,NULL);dCHK(err);
    err = PetscOptionsReal("-ellip_p","p in p-Laplacian","",prm->p,&prm->p,NULL);dCHK(err);
    err = PetscOptionsReal("-ellip_eps","Regularization in p-Laplacian","",prm->epsilon,&prm->epsilon,NULL);dCHK(err);
    err = PetscOptionsEnum("-ellip_f_qmethod","Quadrature method for residual evaluation/matrix-free","",dQuadratureMethods,(PetscEnum)elp->function_qmethod,(PetscEnum*)&elp->function_qmethod,NULL);dCHK(err);
    err = PetscOptionsEnum("-ellip_jac_qmethod","Quadrature to use for Jacobian assembly","",dQuadratureMethods,(PetscEnum)elp->jacobian_qmethod,(PetscEnum*)&elp->jacobian_qmethod,NULL);dCHK(err);
    {
      dBool flg; dInt n = ALEN(prm->dirichlet);
      err = PetscOptionsIntArray("-dirichlet","List of boundary sets on which to impose Dirichlet conditions","",prm->dirichlet,&n,&flg);dCHK(err);
      if (flg) for (dInt i=n; i<ALEN(prm->dirichlet); i++) prm->dirichlet[i] = 0; /* Clear out any leftover values */
    }
    err = PetscOptionsInt("-exact","Exact solution choice","",exact,&exact,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_a","First scale parameter","",exc->a,&exc->a,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_b","Second scale parameter","",exc->b,&exc->b,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_c","Third scale parameter","",exc->c,&exc->c,NULL);dCHK(err);
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
  default: dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Exact solution %d not implemented",exact);
  }

  err = dMeshCreate(elp->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);
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
  for (dInt i=0; i<ALEN(elp->param.dirichlet) && elp->param.dirichlet[i]>0; i++) {
    err = dFSRegisterBoundary(fs,elp->param.dirichlet[i],dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  }
  err = dFSSetFromOptions(fs);dCHK(err);
  err = PetscObjectSetName((PetscObject)fs,"dFS_0");dCHK(err);
  err = PetscObjectSetName((PetscObject)mesh,"dMesh_0");dCHK(err);
  elp->fs = fs;
  err = dFSCreateGlobalVector(elp->fs,&elp->rhs);dCHK(err);
  err = VecZeroEntries(elp->rhs);dCHK(err);
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
    if (eval == EVAL_FUNCTION) {err = dRulesetIteratorAddFS(iter,elp->fs);dCHK(err);}
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
  err = VecDestroy(&elp->rhs);dCHK(err);
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
}

static inline void EllipPointwiseFunction(struct EllipParam *prm,struct EllipExact *exact,struct EllipExactCtx *exactctx,
                                          const dReal x[3],dReal weight,const dScalar u[1],const dScalar Du[3],const dScalar w[1],
                                          struct EllipStore *st,dScalar v[1],dScalar Dv[3])
{
  dScalar f[1];
  EllipPointwiseComputeStore(prm,x,u,Du,st);
  exact->forcing(exactctx,prm,x,f);
  v[0] = - weight * (w[0] + f[0]);    /* Coefficient of \a v in weak form */
  for (dInt i=0; i<3; i++) Dv[i] = weight * st->eta * Du[i]; /* Coefficient of Dv in weak form */
}

static inline void EllipPointwiseJacobian(struct EllipParam dUNUSED *prm,const struct EllipStore *restrict st,dReal weight,
                                          const dScalar u[restrict static 1],const dScalar Du[restrict static 3],
                                          dScalar v[restrict static 1],dScalar Dv[restrict static 3])
{
  const dScalar dotw = dDotScalar3(st->sqrt_mdeta_Du,Du)*weight;
  const dReal etaw = st->eta*weight;
  v[0] = weight * 0 * u[0];
  for (dInt i=0; i<3; i++) Dv[i] = etaw*Du[i] - dotw*st->sqrt_mdeta_Du[i];
}

static dErr EllipFunction(SNES dUNUSED snes,Vec gx,Vec gy,void *ctx)
{
  Ellip elp = ctx;
  Vec Coords,gz = elp->rhs;
  dRulesetIterator iter;
  dReal mineta=1e10,maxeta=1e-10;
  dErr err;

  dFunctionBegin;
  err = VecZeroEntries(gy);dCHK(err);
  err = dFSGetGeometryVectorExpanded(elp->fs,&Coords);dCHK(err);
  err = EllipGetRegionIterator(elp,EVAL_FUNCTION,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL, gx,dFS_INHOMOGENEOUS,gy,dFS_INHOMOGENEOUS, gz,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar *x,*dx,*u,*du,*v,*dv,*w;
    dInt Q;
    struct EllipStore *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,&v,&dv, &w,NULL,NULL,NULL);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      EllipPointwiseFunction(&elp->param,&elp->exact,&elp->exactctx,&x[i*3],jw[i],&u[i],&du[i*3],&w[i],&stash[i],&v[i],&dv[i*3]);
      maxeta = dMax(maxeta,stash[i].eta);
      mineta = dMin(mineta,stash[i].eta);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES,NULL,NULL,v,dv,NULL,NULL);dCHK(err);
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
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_HOMOGENEOUS,gy,dFS_HOMOGENEOUS,NULL,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar *x,*dx,*u,*du,*v,*dv;
    dInt Q;
    const struct EllipStore *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,&v,&dv, NULL,NULL,NULL,NULL);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      EllipPointwiseJacobian(&elp->param,&stash[i],jw[i],&u[i],&du[i*3],&v[i],&dv[i*3]);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES,NULL,NULL,v,dv,NULL,NULL);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = PetscLogEventEnd(LOG_EllipShellMatMult,J,0,0,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr EllipJacobian(SNES dUNUSED snes,Vec gx,Mat *J,Mat *Jp,MatStructure *structure,void *ctx)
{
  Ellip            elp = ctx;
  Vec              Coords;
  dRulesetIterator iter;
  dErr             err;
  dScalar          *Kflat;

  dFunctionBegin;
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
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,NULL,NULL);dCHK(err);
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
  Vec Coords;
  dRulesetIterator iter;

  dFunctionBegin;
  err = dNormsStart(errorNorms,gerrorNorms);dCHK(err);
  err = dFSGetGeometryVectorExpanded(elp->fs,&Coords);dCHK(err);
  err = EllipGetRegionIterator(elp,EVAL_FUNCTION,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL, gx,dFS_INHOMOGENEOUS,NULL, NULL,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[3][3],(*u)[1],(*du)[1][3];
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[1],duu[1][3];
      elp->exact.solution(&elp->exactctx,&elp->param,x[i],uu,(dScalar*)duu);
      err = dNormsUpdate(errorNorms,gerrorNorms,jw[i],1,uu,u[i],&duu[0][0],&du[i][0][0]);dCHK(err);
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = dNormsFinish(errorNorms,gerrorNorms);dCHK(err);
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

int main(int argc,char *argv[])
{
  char mtype[256] = MATAIJ;
  Ellip elp;
  dFS fs;
  MPI_Comm comm;
  Mat J,Jp;
  Vec r,x,soln;
  SNES snes;
  dBool  nojshell,nocheck,viewdhm,eigen;
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
    err = PetscOptionsBool("-eigen","Solve the eigen-problem, always use -exact 4 for this","",eigen=dFALSE,&eigen,NULL);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);
  if (nojshell) {
    /* Use the preconditioning matrix in place of the Jacobian.  This will NOT converge unless the elements are actually
    * Q1 (bdeg=1).  This option is nullified by -snes_mf_operator which will still only use the assembled Jacobian for
    * preconditioning. */
    J = Jp;
    err = PetscObjectReference((PetscObject)J);dCHK(err);
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
  if (eigen) {
    PF pfeig;
    struct PLapEig pleig;
    dReal mu,q,p,norm,kp,scale;
    dInt n,max_it;
    dBool eig_monitor;
    Vec z;
    p = elp->param.p;
    err = PetscOptionsBegin(comm,NULL,"Inverse iteration eigensolver parameters",__FILE__);dCHK(err); {
      err = PetscOptionsReal("-eig_mu","an arbitrary positive number, or 0 for auto","",mu=0,&mu,NULL);dCHK(err);
      err = PetscOptionsReal("-eig_q","q should be chosen close to p","",q=p-0.1,&q,NULL);dCHK(err);
      err = PetscOptionsInt("-eig_max_it","maximum iterations of eigen-solver","",max_it=10,&max_it,NULL);dCHK(err);
      err = PetscOptionsBool("-eig_monitor","monitor eigenvalue estimates","",eig_monitor=dFALSE,&eig_monitor,NULL);dCHK(err);
    } err = PetscOptionsEnd();dCHK(err);
    pleig.mu = mu;
    pleig.q = q;
    err = VecGetSize(x,&n);dCHK(err);
    err = PFCreate(comm,1,1,&pfeig);dCHK(err);
    err = PFSet(pfeig,PFApply_PLapEig,NULL,NULL,NULL,&pleig);dCHK(err);
    err = VecSet(elp->rhs,1);dCHK(err);
    err = VecSet(x,1);dCHK(err);
    err = SNESSolve(snes,NULL,x);dCHK(err);
    err = VecDuplicate(x,&z);dCHK(err);
    err = VecNorm(x,NORM_MAX,&norm);dCHK(err);
    kp = pow(norm,1-p);
    if (mu == 0) pleig.mu = mu = kp;
    scale = pow(mu/kp,1/(p-q));
    err = dPrintf(comm,"Torsion function norm %12.6e  kp %12.6e  scaling by %12.6e for supersolution\n",norm,kp,scale);
    err = VecScale(x,scale/norm);dCHK(err);
    for (dInt i=0; i<max_it; i++) {
      err = PFApplyVec(pfeig,x,elp->rhs);dCHK(err); // rhs = mu * u^{q-1}
      err = VecCopy(x,z);dCHK(err);                 // save
      err = SNESSolve(snes,NULL,x);dCHK(err);
      err = VecNorm(x,NORM_MAX,&norm);dCHK(err);
      //mu_q = mu * pow(norm,q-p);
      //err = VecScale(x,1./norm);dCHK(err);
      if (eig_monitor) {
        dReal xnorm2,znorm2,lambda_p;
        lambda_p = mu*pow(norm,q-p); // paper
        //lambda_p = norm / pow(mu,p-1);
        err = VecAYPX(z,-1,x);dCHK(err);
        err = VecNorm(x,NORM_2,&xnorm2);dCHK(err);
        err = VecNorm(z,NORM_2,&znorm2);dCHK(err);
        err = dPrintf(comm,"%3d EIG solution norm2 % 12.6e  difference norm2 % 12.6e  eigenvalue estimate % 18.12e\n",i,xnorm2,znorm2,lambda_p);dCHK(err);
      }
    }
    err = VecNorm(x,NORM_MAX,&norm);dCHK(err);
    err = VecScale(x,1./norm);dCHK(err);
    err = PFDestroy(&pfeig);dCHK(err);
  } else {
    err = VecZeroEntries(x);dCHK(err);
    err = SNESSolve(snes,NULL,x);dCHK(err);
  }
  if (!nocheck) {
    dReal anorm[3],inorm[3],enorm[3],gnorm[3];
    err = EllipErrorNorms(elp,x,enorm,gnorm);dCHK(err);
    err = dNormsAlgebraicScaled(anorm,r);dCHK(err); // Algebraic residual for solution, scaled by number of degrees of freedom
    err = VecWAXPY(r,-1,soln,x);dCHK(err);
    err = dNormsAlgebraicScaled(inorm,r);dCHK(err); // Algebraic difference between interpolated exact solution and computed solution
    err = dPrintf(comm,"Algebraic residual        |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",anorm[0],anorm[1],anorm[2]);dCHK(err);
    err = dPrintf(comm,"Interpolation residual    |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",inorm[0],inorm[1],inorm[2]);dCHK(err);
    err = dPrintf(comm,"Pointwise solution error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",enorm[0],enorm[1],enorm[2]);dCHK(err);
    err = dPrintf(comm,"Pointwise gradient error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",gnorm[0],gnorm[1],gnorm[2]);dCHK(err);
  }

  if (viewdhm) {
    dViewer viewer;
    err = PetscViewerCreate(comm,&viewer);dCHK(err);
    err = PetscViewerSetType(viewer,PETSCVIEWERDHM);dCHK(err);
    err = PetscViewerFileSetName(viewer,"plapeig.dhm");dCHK(err);
    err = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);dCHK(err);
    err = VecView(x,viewer);dCHK(err);
    err = PetscViewerDestroy(&viewer);dCHK(err);
  }

  err = VecDestroy(&r);dCHK(err);
  err = VecDestroy(&x);dCHK(err);
  err = VecDestroy(&soln);dCHK(err);
  err = SNESDestroy(&snes);dCHK(err);
  err = MatDestroy(&J);dCHK(err);
  err = MatDestroy(&Jp);dCHK(err);
  err = EllipDestroy(elp);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
