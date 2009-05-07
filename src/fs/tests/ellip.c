static const char help[] = "Solve a scalar elliptic problem, a regularized p-Laplacian using dual order hp elements.\n"
  "The model problem is\n"
  "  \\int_\\Omega (\\eta Dv \\cdot Du - f v) - \\int_\\Gamma v (\\eta Du \\cdot n) = 0\n"
  "where\n"
  "  \\eta(u) = (\\epsilon + 1/2 Du . Du)^(p-2)\n"
  "  (\\eta Du \\cdot n) = known OR function of u OR self (\"No condition\" outflow)\n\n";

#include "dohpfs.h"
#include "dohpvec.h"
#include "petscsnes.h"


struct EllipParam {
  dReal epsilon;
  dReal p;
  dTruth onlyproject;
};

struct EllipExactCtx {
  dReal a,b,c;
};
struct EllipExact {
  void (*solution)(const struct EllipExactCtx*,const struct EllipParam*,const dReal x[3],dScalar u[1],dScalar du[3]);
  void (*forcing)(const struct EllipExactCtx*,const struct EllipParam*,const dReal x[3],dScalar f[1]);
};

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
  f[0] = pow(a,2)*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-2 + p))*cos(a*x)*exp(b*y)*sin(c*z) + pow(c,2)*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-2 + p))*cos(a*x)*exp(b*y)*sin(c*z) - pow(b,2)*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-2 + p))*cos(a*x)*exp(b*y)*sin(c*z) + b*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-3 + p))*(2 - p)*(pow(b,3)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y) + b*pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y) + b*pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y))*cos(a*x)*exp(b*y)*sin(c*z) - a*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-3 + p))*(2 - p)*(pow(a,3)*pow(sin(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) - a*pow(b,2)*pow(sin(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x) - a*pow(c,2)*pow(cos(c*z),2)*cos(a*x)*exp(2*b*y)*sin(a*x))*exp(b*y)*sin(a*x)*sin(c*z) + c*pow((pow(e,2) + pow(a,2)*pow(sin(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(b,2)*pow(cos(a*x),2)*pow(sin(c*z),2)*exp(2*b*y)/2 + pow(c,2)*pow(cos(a*x),2)*pow(cos(c*z),2)*exp(2*b*y)/2),(-3 + p))*(2 - p)*(-pow(c,3)*pow(cos(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) + c*pow(a,2)*pow(sin(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z) + c*pow(b,2)*pow(cos(a*x),2)*cos(c*z)*exp(2*b*y)*sin(c*z))*cos(a*x)*cos(c*z)*exp(b*y);
#if 0
  /** The following expression was generated with Maxima:
  * u : cos(a*x)+exp(b*y)+sin(c*z); ux : diff(u,x); uy : diff(u,y); uz : diff(u,z);
  * gamma : 1/2 * (ux^2 + uy^2 + uz^2); eta : (eps^2 + gamma)^(p-2);
  * optimize(diff(eta*ux,x) + diff(eta*uy,y) + diff(eta*uz,z));
  **/
  const dReal t1=p-2,t2=a*x,t3=cos(t2),t4=dSqr(sin(t2)),t5=a*a,t6=b*b,t7=c*c,t8=c*z,t9=dSqr(cos(t8)),
    t10=(t7*t9+t6*exp(2*b*y)+t5*t4)/2+dSqr(eps),t11=pow(t10,p-3),t12=pow(t10,t1),t13=sin(t8);
  f[0] = -t7*t12*t13 - c*c*c*c*t1*t9*t11*t13 + t6*exp(b*y)*t12 - t5*t3*t12 + b*b*b*b*t1*exp(3*b*y)*t11 - a*a*a*a*t1*t3*t4*t11;
#endif
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
  /** Maxima optimized
  * u : x^2 + y^2 + 1-z^2; ux : diff(u,x); uy : diff(u,y); uz : diff(u,z);
  * gamma : 1/2 * (ux^2 + uy^2 + uz^2); eta : (eps^2 + gamma)^(p-2);
  * optimize(diff(eta*ux,x) + diff(eta*uy,y) + diff(eta*uz,z));
  *   (\%o38) \mathbf{block}\;\left(\left[ \mathrm{\%1},\linebreak[0]\mathrm{\%2},\linebreak[0]\mathrm{\%3},\linebreak[0]\mathrm{\%4},\linebreak[0]\mathrm{\%5},\linebreak[0]\mathrm{\%6} \right] ,\linebreak[0]\mathrm{\%1}:p-2,\linebreak[0]\mathrm{\%2}:x^{2},\linebreak[0]\mathrm{\%3}:y^{2},\linebreak[0]\mathrm{\%4}:z^{2},\linebreak[0]\mathrm{\%5}:\ifracn{4\*\mathrm{\%4}+4\*\mathrm{\%3}+4\*\mathrm{\%2}}{2}+\mathrm{eps}^{2},\linebreak[0]\mathrm{\%6}:\iexpt{\mathrm{\%5}}{p-3},\linebreak[0]2\*\mathrm{\%5}^{\mathrm{\%1}}-8\*\mathrm{\%1}\*\mathrm{\%4}\*\mathrm{\%6}+8\*\mathrm{\%1}\*\mathrm{\%3}\*\mathrm{\%6}+8\*\mathrm{\%1}\*\mathrm{\%2}\*\mathrm{\%6}\right)
  **/
#if 0
  const dReal t1=p-2,t2=x*x,t3=y*y,t4=z*z,t5=2*(t4+t3+t2)+e*e,t6=pow(t5,p-3);
  f[0] = -(2*pow(t5,t1) - 8*t1*t6*(t4+t3+t2));
#else
  f[0] = -8*pow(c,3)*pow(z,2)*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-3 + p))*(2 - p) + 8*pow(a,3)*pow(x,2)*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-3 + p))*(2 - p) + 8*pow(b,3)*pow(y,2)*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-3 + p))*(2 - p) - 2*a*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-2 + p)) - 2*b*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-2 + p)) + 2*c*pow((pow(e,2) + 2*pow(a,2)*pow(x,2) + 2*pow(b,2)*pow(y,2) + 2*pow(c,2)*pow(z,2)),(-2 + p));
#endif
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
  const dUNUSED dReal x = xyz[0],y = xyz[1],z = xyz[2],a = ctx->a,b = ctx->b,c = ctx->c,e = prm->epsilon,p = prm->p;
  f[0] = pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-3 + p))*(2 - p)*(-2*c*x*z + b*pow(y,2))*(-2*c*x*(-2*c*x*z + b*pow(y,2)) - 2*c*z*(c*(1 - pow(z,2)) + 3*a*pow(x,2)) + 4*z*pow(b,2)*pow(y,2)) + pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-3 + p))*(2 - p)*(c*(1 - pow(z,2)) + 3*a*pow(x,2))*(-2*c*z*(-2*c*x*z + b*pow(y,2)) + 6*a*x*(c*(1 - pow(z,2)) + 3*a*pow(x,2))) + 2*b*y*z*pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-3 + p))*(2 - p)*(2*b*y*(-2*c*x*z + b*pow(y,2)) + 4*y*pow(b,2)*pow(z,2)) - 6*a*x*pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-2 + p)) - 2*b*z*pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-2 + p)) + 2*c*x*pow((pow(e,2) + pow((-2*c*x*z + b*pow(y,2)),2)/2 + pow((c*(1 - pow(z,2)) + 3*a*pow(x,2)),2)/2 + 2*pow(b,2)*pow(y,2)*pow(z,2)),(-2 + p));
}

struct EllipStore {
  dReal eta,deta;
  dReal Du[3];
};

typedef struct EllipCtx *Ellip;
struct EllipCtx {
  MPI_Comm              comm;
  struct EllipParam     param;
  struct EllipExact     exact;
  struct EllipExactCtx  exactctx;
  struct EllipStore    *store;
  dInt                 *storeoff;
  dJacobi               jac;
  dMesh                 mesh;
  dFS                   fs;
  Vec                   x,y;
  dInt                  constBDeg,nominalRDeg;
  dTruth                errorview;
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

  elp->constBDeg = 4;
  elp->nominalRDeg = 0;

  prm = &elp->param;
  prm->p           = 2.0;       /* p in p-Laplacian */
  prm->epsilon     = 1.0;
  prm->onlyproject = dFALSE;

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
  dMeshTag rtag,dtag;
  dInt exact;
  dErr err;

  dFunctionBegin;
  exact = 0; exc->a = exc->b = exc->c = 1;
  err = PetscOptionsBegin(elp->comm,NULL,"Elliptic (p-Laplacian) options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-const_bdeg","Use constant isotropic degree on all elements","",elp->constBDeg,&elp->constBDeg,NULL);dCHK(err);
    err = PetscOptionsInt("-nominal_rdeg","Nominal rule degree (will be larger if basis requires it)","",elp->nominalRDeg,&elp->nominalRDeg,NULL);dCHK(err);
    err = PetscOptionsTruth("-error_view","View errors","",elp->errorview,&elp->errorview,NULL);dCHK(err);
    err = PetscOptionsReal("-ellip_p","p in p-Laplacian","",prm->p,&prm->p,NULL);dCHK(err);
    err = PetscOptionsReal("-ellip_epsilon","Regularization in p-Laplacian","",prm->epsilon,&prm->epsilon,NULL);dCHK(err);
    err = PetscOptionsTruth("-onlyproject","Actually just do a projection","",prm->onlyproject,&prm->onlyproject,NULL);dCHK(err);
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
    case 2:
      elp->exact.solution = EllipExact_2_Solution;
      elp->exact.forcing = EllipExact_2_Forcing;
      break;
    case 3:
      elp->exact.solution = EllipExact_3_Solution;
      elp->exact.forcing = EllipExact_3_Forcing;
      break;
    default: dERROR(1,"Exact solution %d not implemented");
  }

  err = dMeshCreate(elp->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);dCHK(err);
  elp->mesh = mesh;
  domain = 0;                   /* Root set in MOAB */

  err = dJacobiCreate(elp->comm,&jac);dCHK(err);
  err = dJacobiSetDegrees(jac,9,2);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);
  elp->jac = jac;

  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"ellip_rule_degree",elp->nominalRDeg,&rtag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"ellip_efs_degree",elp->constBDeg,&dtag);dCHK(err);

  err = dFSCreate(elp->comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,0);dCHK(err);
  err = dFSSetRuleTag(fs,jac,rtag);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  err = dFSRegisterBoundary(fs,100,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSRegisterBoundary(fs,200,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSRegisterBoundary(fs,300,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSSetFromOptions(fs);dCHK(err);
  elp->fs = fs;

  err = dFSCreateExpandedVector(fs,&elp->x);dCHK(err);
  err = VecDuplicate(elp->x,&elp->y);dCHK(err);

  {                             /* Allocate space for stored values */
    dInt n;
    s_dRule *rule;
    err = dFSGetElements(fs,&n,NULL,&rule,NULL,NULL,NULL);dCHK(err);
    err = dMallocA(n+1,&elp->storeoff);dCHK(err);
    elp->storeoff[0] = 0;
    for (dInt i=0; i<n; i++) {
      dInt q;
      err = dRuleGetSize(&rule[i],NULL,&q);dCHK(err);
      elp->storeoff[i+1] = elp->storeoff[i] + q;
    }
    err = dMallocA(elp->storeoff[n],&elp->store);dCHK(err);
    err = dMemzero(elp->store,elp->storeoff[n]*sizeof(elp->store[0]));dCHK(err);
    err = dFSRestoreElements(fs,&n,NULL,&rule,NULL,NULL,NULL);dCHK(err);
  }

  dFunctionReturn(0);
}

static dErr EllipDestroy(Ellip elp)
{
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(elp->fs);dCHK(err);
  err = dJacobiDestroy(elp->jac);dCHK(err);
  err = dMeshDestroy(elp->mesh);dCHK(err);
  err = dFree(elp->storeoff);dCHK(err);
  err = dFree(elp->store);dCHK(err);
  if (elp->x) {err = VecDestroy(elp->x);dCHK(err);}
  if (elp->y) {err = VecDestroy(elp->y);dCHK(err);}
  err = dFree(elp);dCHK(err);
  dFunctionReturn(0);
}

static inline void EllipPointwiseComputeStore(struct EllipParam *prm,const dReal dUNUSED x[3],const dScalar dUNUSED u[1],const dScalar Du[3],struct EllipStore *st)
{
  dReal gamma,espg,p=prm->p;
  gamma = 0.5 * (dSqr(Du[0]) + dSqr(Du[1]) + dSqr(Du[2]));
  espg = dSqr(prm->epsilon) + gamma;
  st->eta = pow(espg,p-2);
  st->deta = (p-2) * st->eta / espg;
  st->Du[0] = Du[0]; st->Du[1] = Du[1]; st->Du[2] = Du[2];
}

static inline void EllipPointwiseFunction(struct EllipParam *prm,struct EllipExact *exact,struct EllipExactCtx *exactctx,
                                          const dReal x[3],dReal weight,const dScalar u[1],const dScalar Du[3],
                                          struct EllipStore *st,dScalar v[1],dScalar Dv[3])
{
  dScalar f[1];
  EllipPointwiseComputeStore(prm,x,u,Du,st);
  exact->forcing(exactctx,prm,x,f);
  v[0] = - weight * f[0];       /* Coefficient of \a v in weak form */
  if (prm->onlyproject) {
    v[0] += weight * u[0];
    Dv[0] = Dv[1] = Dv[2] = 0;
  } else {
    for (dInt i=0; i<3; i++) {
      Dv[i] = weight * st->eta * Du[i]; /* Coefficient of Dv in weak form */
  }
  }
}

static inline void EllipPointwiseJacobian(struct EllipParam dUNUSED *prm,const struct EllipStore *restrict st,dReal weight,
                                          const dScalar dUNUSED u[restrict static 1],const dScalar Du[restrict static 3],
                                          dScalar v[restrict static 1],dScalar Dv[restrict static 3])
{
  const dScalar dot = dDotScalar3(st->Du,Du);
  const dReal etaw = st->eta*weight,dotdetaw = dot*st->deta*weight;
  if (prm->onlyproject) {
    v[0] = weight * u[0];
    Dv[0] = Dv[1] = Dv[2] = 0;
  } else {
    v[0] = 0;
    Dv[0] = etaw*Du[0] + dotdetaw*st->Du[0];
    Dv[1] = etaw*Du[1] + dotdetaw*st->Du[1];
    Dv[2] = etaw*Du[2] + dotdetaw*st->Du[2];
  }
}

static dErr EllipFunction(SNES dUNUSED snes,Vec gx,Vec gy,void *ctx)
{
  Ellip elp = ctx;
  dFS fs = elp->fs;
  dInt n,*off,*geomoff;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*restrict geom)[3],(*restrict q)[3],(*restrict jinv)[3][3],*restrict jw;
  dScalar *x,*y,*restrict u,*restrict v,*restrict du,*restrict dv;
  dErr err;

  dFunctionBegin;
  err = dFSGlobalToExpanded(fs,gx,elp->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(elp->x,&x);dCHK(err);
  err = VecZeroEntries(elp->y);dCHK(err);
  err = VecGetArray(elp->y,&y);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  for (dInt e=0; e<n; e++) {
    dInt Q;
    err = dRuleComputeGeometry(&rule[e],(const dReal(*)[3])(geom+geomoff[e]),q,jinv,jw);dCHK(err);
    err = dRuleGetSize(&rule[e],0,&Q);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,x+off[e],u,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,x+off[e],du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      struct EllipStore *restrict st = &elp->store[elp->storeoff[e]+i];
      EllipPointwiseFunction(&elp->param,&elp->exact,&elp->exactctx,q[i],jw[i],&u[i],&du[i*3],st,&v[i],&dv[i*3]);
    }
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,v,y+off[e],dAPPLY_INTERP_TRANSPOSE,ADD_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,dv,y+off[e],dAPPLY_GRAD_TRANSPOSE,ADD_VALUES);dCHK(err);
  }
  err = dFSRestoreWorkspace(fs,__func__,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
#if 0
  if (0) {
    dMeshESH *sets;
    dInt nsets;
    err = dFSGetBoundarySets(fs,100,&nsets,&sets);dCHK(err);
    for (dInt s=0; s<nsets; s++) {
      err = dFSGetElements(fs,sets[s],&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
      for (dInt e=0; e<n; e++) {
        // integrate weak forms over faces
      }
      err = dFSRestoreElements(fs,sets[s],&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
    }
  }
#endif
  err = VecRestoreArray(elp->x,&x);dCHK(err);
  err = VecRestoreArray(elp->y,&y);dCHK(err);
  err = dFSExpandedToGlobal(fs,elp->y,gy,dFS_INHOMOGENEOUS,ADD_VALUES);dCHK(err);
  dFunctionReturn(0);
}

static dErr EllipShellMatMult(Mat J,Vec gx,Vec gy)
{
  Ellip elp;
  dFS fs;
  dInt n,*off,*geomoff;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*restrict geom)[3],(*restrict q)[3],(*restrict jinv)[3][3],*restrict jw;
  dScalar *x,*y,*restrict u,*restrict v,*restrict du,*restrict dv;
  dErr err;

  dFunctionBegin;
  err = MatShellGetContext(J,(void**)&elp);dCHK(err);
  fs = elp->fs;
  err = dFSGlobalToExpanded(fs,gx,elp->x,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(elp->x,&x);dCHK(err);
  err = VecZeroEntries(elp->y);dCHK(err);
  err = VecGetArray(elp->y,&y);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  for (dInt e=0; e<n; e++) {
    dInt Q;
    err = dRuleComputeGeometry(&rule[e],(const dReal(*)[3])(geom+geomoff[e]),q,jinv,jw);dCHK(err);
    err = dRuleGetSize(&rule[e],0,&Q);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,x+off[e],u,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,x+off[e],du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      struct EllipStore *restrict st = &elp->store[elp->storeoff[e]+i];
      EllipPointwiseJacobian(&elp->param,st,jw[i],&u[i],&du[i*3],&v[i],&dv[i*3]);
    }
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,v,y+off[e],dAPPLY_INTERP_TRANSPOSE,ADD_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,dv,y+off[e],dAPPLY_GRAD_TRANSPOSE,ADD_VALUES);dCHK(err);
  }
  err = dFSRestoreWorkspace(fs,__func__,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(elp->x,&x);dCHK(err);
  err = VecRestoreArray(elp->y,&y);dCHK(err);
  err = dFSExpandedToGlobal(fs,elp->y,gy,dFS_HOMOGENEOUS,ADD_VALUES);dCHK(err);
  dFunctionReturn(0);
}

static dErr EllipJacobian(SNES dUNUSED snes,Vec gx,Mat *J,Mat *Jp,MatStructure *structure,void *ctx)
{
  Ellip elp = ctx;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*nx)[3];
  dScalar *x;
  dFS fs = elp->fs;
  dInt n,*off,*geomoff;
  dReal (*geom)[3];
  dErr err;

  dFunctionBegin;
  err = MatZeroEntries(*Jp);dCHK(err);
  err = dFSGlobalToExpanded(fs,gx,elp->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(elp->x,&x);dCHK(err);
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
          const dScalar (*uc)[1] = (const dScalar(*)[1])x+off[e]; /* function values, indexed at subelement corners \c uc[c[#]][0] */
          const dReal (*qx)[3],*jw,(*basis)[8],(*deriv)[8][3];
          dInt qn;
          dScalar K[8][8];
          err = dMemzero(K,sizeof(K));dCHK(err);
          err = dQ1HexComputeQuadrature(corners,&qn,&qx,&jw,(const dReal**)&basis,(const dReal**)&deriv);dCHK(err);
          for (dInt lq=0; lq<qn; lq++) { /* loop over quadrature points */
            struct EllipStore st;
            { /* Set up store */
              dReal st_u[1] = {0},st_Du[3] = {0,0,0};
              for (dInt lp=0; lp<8; lp++) { /* Evaluate function values and gradients at this quadrature point */
                st_u[0] += basis[lq][lp] * uc[c[lp]][0];
                st_Du[0] += deriv[lq][lp][0] * uc[c[lp]][0];
                st_Du[1] += deriv[lq][lp][1] * uc[c[lp]][0];
                st_Du[2] += deriv[lq][lp][2] * uc[c[lp]][0];
              }
              EllipPointwiseComputeStore(&elp->param,qx[lq],st_u,st_Du,&st);
            }
            for (dInt ltest=0; ltest<8; ltest++) {              /* Loop over test basis functions (corners) */
              for (dInt lp=0; lp<8; lp++) {                     /* loop over trial basis functions (corners) */
                const dReal *u = &basis[lq][lp],*Du = deriv[lq][lp];
                dReal v[1],Dv[3];
                EllipPointwiseJacobian(&elp->param,&st,jw[lq],u,Du,v,Dv);
#if 1
                K[ltest][lp] += //basis[lq][ltest] * v[0]
                  + deriv[lq][ltest][0] * Dv[0]
                  + deriv[lq][ltest][1] * Dv[1]
                  + deriv[lq][ltest][2] * Dv[2];
#else
                K[ltest][lp] += basis[lq][ltest] * jw[lq] * basis[lq][lp];
#endif
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
  err = VecRestoreArray(elp->x,&x);dCHK(err);

  err = MatAssemblyBegin(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  *structure = DIFFERENT_NONZERO_PATTERN;
  dFunctionReturn(0);
}

static dErr EllipErrorNorms(Ellip elp,Vec gx,dReal errorNorms[static 3],dReal gerrorNorms[static 3])
{
  dFS fs = elp->fs;
  dInt n,*off,*geomoff;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*geom)[3],(*q)[3],(*jinv)[3][3],*jw;
  dScalar *x,*u,(*du)[3];
  dErr err;

  dFunctionBegin;
  err = dMemzero(errorNorms,3*sizeof(errorNorms));dCHK(err);
  err = dMemzero(gerrorNorms,3*sizeof(gerrorNorms));dCHK(err);
  err = dFSGlobalToExpanded(fs,gx,elp->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(elp->x,&x);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&q,&jinv,&jw,&u,NULL,(dReal**)&du,NULL);dCHK(err);
  for (dInt e=0; e<n; e++) {
    dInt Q;
    err = dRuleComputeGeometry(&rule[e],(const dReal(*)[3])(geom+geomoff[e]),q,jinv,jw);dCHK(err);
    err = dRuleGetSize(&rule[e],0,&Q);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,x+off[e],u,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,1,x+off[e],&du[0][0],dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[1],duu[3],r[1],gr[3],grsum;             /* Scalar problem */
      elp->exact.solution(&elp->exactctx,&elp->param,q[i],uu,duu);
      r[0] = u[i] - uu[0];   /* Function error at point */
      gr[0] = du[i][0] - duu[0]; /* Gradient error at point */
      gr[1] = du[i][1] - duu[1];
      gr[2] = du[i][2] - duu[2];
      if (elp->errorview) {
        printf("e,q = %3d %3d (% 5f,% 5f,% 5f) dohp %10.2e   exact %10.2e   error %10.e\n",e,i,q[i][0],q[i][1],q[i][2],u[i],uu[0],r[0]);
      }
      grsum = dDotScalar3(gr,gr);
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
  }
  err = dFSRestoreWorkspace(fs,__func__,&q,&jinv,&jw,&u,NULL,NULL,NULL);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(elp->x,&x);dCHK(err);
  errorNorms[1] = dSqrt(errorNorms[1]);
  gerrorNorms[1] = dSqrt(gerrorNorms[1]);
  dFunctionReturn(0);
}

static dErr EllipGetSolutionVector(Ellip elp,Vec *insoln)
{
  dErr err;
  Vec soln,xc,cvec;
  dScalar *x,*coords;
  dInt n,bs;

  dFunctionBegin;
  *insoln = 0;
  err = dFSCreateGlobalVector(elp->fs,&soln);dCHK(err);
  err = VecDohpGetClosure(soln,&xc);dCHK(err);
  err = dFSGetCoordinates(elp->fs,&cvec);dCHK(err);
  err = VecGetLocalSize(xc,&n);dCHK(err);
  err = VecGetBlockSize(xc,&bs);dCHK(err);
  err = VecGetArray(xc,&x);dCHK(err);
  err = VecGetArray(cvec,&coords);dCHK(err);
  for (dInt i=0; i<n/bs; i++) {
    dScalar du_unused[3*bs];
    elp->exact.solution(&elp->exactctx,&elp->param,&coords[3*i],&x[i*bs],du_unused);
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
  Ellip elp;
  dFS fs;
  MPI_Comm comm;
  PetscViewer viewer;
  Mat J,Jp;
  Vec r,x,soln;
  SNES snes;
  dTruth nojshell,nocheck;
  dErr err;

  err = PetscInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;

  err = EllipCreate(comm,&elp);dCHK(err);
  err = EllipSetFromOptions(elp);dCHK(err);
  fs = elp->fs;

  err = dFSCreateGlobalVector(fs,&r);dCHK(err);
  err = dFSGetMatrix(fs,MATSEQAIJ,&Jp);dCHK(err);
  err = MatSetOptionsPrefix(Jp,"q1");dCHK(err);
  err = MatSeqAIJSetPreallocation(Jp,27,NULL);dCHK(err);

  err = PetscOptionsBegin(elp->comm,NULL,"Elliptic solver options",__FILE__);dCHK(err); {
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
    err = MatCreateShell(comm,m,n,M,N,elp,&J);dCHK(err);
    err = MatSetOptionsPrefix(J,"j");dCHK(err);
    err = MatShellSetOperation(J,MATOP_MULT,(void(*)(void))EllipShellMatMult);dCHK(err);
  }
  err = SNESCreate(comm,&snes);dCHK(err);
  err = SNESSetFunction(snes,r,EllipFunction,elp);dCHK(err);
  err = SNESSetJacobian(snes,J,Jp,EllipJacobian,elp);dCHK(err);
  err = SNESSetTolerances(snes,PETSC_DEFAULT,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);dCHK(err);
  err = SNESSetFromOptions(snes);dCHK(err);
  err = VecZeroEntries(r);dCHK(err);
  err = VecDuplicate(r,&x);dCHK(err);
  err = EllipGetSolutionVector(elp,&soln);dCHK(err);
  {
    Vec sc;
    err = VecDohpGetClosure(soln,&sc);dCHK(err);
    err = dFSInhomogeneousDirichletCommit(elp->fs,sc);dCHK(err);
    err = VecDohpRestoreClosure(soln,&sc);dCHK(err);
  }
  err = VecZeroEntries(x);dCHK(err);
  err = SNESSolve(snes,NULL,x);dCHK(err);
  if (!nocheck) {
    dReal anorm[2],anorminf,inorm[3],enorm[3],gnorm[3];
    err = EllipErrorNorms(elp,x,enorm,gnorm);dCHK(err);
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
  err = EllipDestroy(elp);dCHK(err);
  err = PetscFinalize();dCHK(err);
  return 0;
}
