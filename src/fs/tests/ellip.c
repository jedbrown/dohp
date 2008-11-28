static const char help[] = "Solve a scalar elliptic problem, a regularized p-Laplacian using dual order hp elements.\n"
  "The model problem is\n"
  "  \\int \\eta(u) Du . Dv - \\int f v\n = 0\n"
  "  \\eta(u) = (\\epsilon + 1/2 Du . Du)^((1-p)/2p)\n";

#include "dohpfs.h"
#include "petscsnes.h"


struct EllipParam {
  dReal epsilon;
  dReal exponent;
};

struct EllipExactCtx {
  dReal a,b,c;
};
struct EllipExact {
  void (*solution)(const struct EllipExactCtx*,const struct EllipParam*,const dReal x[3],dReal f[1]);
  void (*forcing)(const struct EllipExactCtx*,const struct EllipParam*,const dReal x[3],dReal f[1]);
};

static void EllipExact_0_Solution(const struct EllipExactCtx *ctx,const struct EllipParam dUNUSED *prm,const dReal x[3],dReal u[1])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c;
  u[0] = cos(a*x[0]) * exp(b*x[1]) * sin(c*x[2]);
}
static void EllipExact_0_Forcing(const struct EllipExactCtx *ctx,const struct EllipParam *prm,const dReal xyz[3],dReal f[1])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2],eps = prm->epsilon,n = prm->exponent;
  /** The following expression was generated with Maxima:
  * u : cos(a*x)+exp(b*y)+sin(c*z); ux : diff(u,x); uy : diff(u,y); uz : diff(u,z);
  * gamma : 1/2 * (ux^2 + uy^2 + uz^2); eta : (eps^2 + gamma)^((1-n)/n);
  * optimize(diff(eta*ux,x) + diff(eta*uy,y) + diff(eta*uz,z));
  **/
  const dReal t1 = a*a,t2 = a*x,t3=cos(t2),t4=dSqr(sin(t2)),t5=b*b,t6=c*c,t7=c*z,t8=dSqr(cos(t7)),
    t9=0.5*(t6*t8+t5*exp(2*b*y)+t1*t4)+eps*eps,t10=1-n,t11=1/n,t12=t10*t11,t13=pow(t9,t12),t14=pow(t9,t12-1),
    t15=sin(t7);
  f[0] = -dSqr(c*c)*t10*t11*t8*t14*t15 - t6*t13*t15 + dSqr(b*b)*t12*t11*exp(3*b*y)*t14 - dSqr(a*a)*t12*t11*t3*t4*t14
    + t5*exp(b*y)*t13 - t1*t3*t13;
}

static void EllipExact_1_Solution(const struct EllipExactCtx *ctx,const struct EllipParam dUNUSED *prm,const dReal xyz[3],dReal u[1])
{
  const dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2];
  u[0] = a*dSqr(x) + b*dSqr(y) + c*(1-dSqr(z));
  u[0] = 10*x+y;                /* DEBUG */
}
static void EllipExact_1_Forcing(const struct EllipExactCtx *ctx,const struct EllipParam *prm,const dReal xyz[3],dReal f[1])
{
  const dUNUSED dReal a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2],eps = prm->epsilon,n = prm->exponent;
  /** Maxima optimized
  * u : x^2 + y^2 + 1-z^2; ux : diff(u,x); uy : diff(u,y); uz : diff(u,z);
  * gamma : 1/2 * (ux^2 + uy^2 + uz^2); eta : (eps^2 + gamma)^((1-n)/n);
  * optimize(diff(eta*ux,x) + diff(eta*uy,y) + diff(eta*uz,z));
  **/
  const dReal t1=x*x,t2=y*y,t3=z*z,t4=2*(t1+t2+t3)+eps*eps,t5=1-n,t6=1/n,t7=t5*t6,t8=pow(t4,t7-1);
  f[0] = 8*t5*t6*t8*(-t3+t2+t1) + 2*pow(t4,t7);
  f[0] = 10*x+y;                /* DEBUG */
}

struct EllipStore {
  dReal eta,deta;
  dReal Du[3];
};

typedef struct EllipCtx *Ellip;
struct EllipCtx {
  MPI_Comm comm;
  struct EllipParam param;
  struct EllipExact exact;
  struct EllipExactCtx exactctx;
  struct EllipStore *store;
  dInt *storeoff;
  dJacobi jac;
  dMesh mesh;
  dFS fs;
  Vec x,y;
  dInt constBDeg,nominalRDeg;
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
  prm->exponent = 1.0;
  prm->epsilon  = 1.0;

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
    err = PetscOptionsReal("-ellip_exponent","Exponent in p-Laplacian","",prm->exponent,&prm->exponent,NULL);dCHK(err);
    err = PetscOptionsReal("-ellip_epsilon","Regularization in p-Laplacian","",prm->epsilon,&prm->epsilon,NULL);dCHK(err);
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
    default: dERROR(1,"Exact solution %d not implemented");
  }

  err = dMeshCreate(elp->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);dCHK(err);
  elp->mesh = mesh;
  domain = 0;                   /* Root set in MOAB */

  err = dJacobiCreate(elp->comm,&jac);dCHK(err);
  err = dJacobiSetDegrees(jac,8,2);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);
  elp->jac = jac;

  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"ellip_rule_degree",elp->nominalRDeg,&rtag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"ellip_efs_degree",elp->constBDeg,&dtag);dCHK(err);

  err = dFSCreate(elp->comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,0);dCHK(err);
  err = dFSSetRuleTag(fs,jac,rtag);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
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

static inline void EllipPointwiseComputeStore(struct EllipParam *prm,const dReal dUNUSED x[3],const dReal dUNUSED u[1],const dReal Du[3],struct EllipStore *st)
{
  dReal gamma,espg,power;
  gamma = 0.5 * (dSqr(Du[0]) + dSqr(Du[1]) + dSqr(Du[2]));
  espg = dSqr(prm->epsilon) + gamma;
  power = (1-prm->exponent)/(2*prm->exponent);
  st->eta = pow(espg,power);
  st->deta = power * st->eta / espg;
  st->Du[0] = Du[0]; st->Du[1] = Du[1]; st->Du[2] = Du[2];
}

static inline void EllipPointwiseFunction(struct EllipParam *prm,struct EllipExact *exact,struct EllipExactCtx *exactctx,
                                          const dReal x[3],dReal weight,const dReal u[1],const dReal Du[3],
                                          struct EllipStore *st,dReal v[1],dReal Dv[3])
{
  dScalar f[1];
  EllipPointwiseComputeStore(prm,x,u,Du,st);
  exact->forcing(exactctx,prm,x,f);
  v[0] = - weight * f[0];       /* Coefficient of \a v in weak form */
  for (dInt i=0; i<3; i++) {
    Dv[i] = weight * st->eta * Du[i]; /* Coefficient of Dv in weak form */
  }
}

static inline void EllipPointwiseJacobian(struct EllipParam dUNUSED *prm,const struct EllipStore *restrict st,dReal weight,
                                          const dReal dUNUSED u[restrict static 1],const dReal Du[restrict static 3],
                                          dScalar v[restrict static 1],dReal Dv[restrict static 3])
{
  const dScalar dot = dDotScalar3(st->Du,Du);
  const dReal etaw = st->eta*weight,dotdetaw = dot*st->deta*weight;
  v[0] = 0;
  Dv[0] = etaw*Du[0] + dotdetaw*st->Du[0];
  Dv[1] = etaw*Du[1] + dotdetaw*st->Du[1];
  Dv[2] = etaw*Du[2] + dotdetaw*st->Du[2];
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
  err = dFSGlobalToExpandedBegin(fs,gx,INSERT_VALUES,elp->x);dCHK(err);
  err = dFSGlobalToExpandedEnd(fs,gx,INSERT_VALUES,elp->x);dCHK(err);
  err = VecGetArray(elp->x,&x);dCHK(err);
  err = VecZeroEntries(elp->y);dCHK(err);
  err = VecGetArray(elp->y,&y);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
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
  err = dFSRestoreWorkspace(fs,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(elp->x,&x);dCHK(err);
  err = VecRestoreArray(elp->y,&y);dCHK(err);
  err = VecZeroEntries(gy);dCHK(err); /* Necessary? */
  err = dFSExpandedToGlobalBegin(fs,elp->y,INSERT_VALUES,gy);dCHK(err);
  err = dFSExpandedToGlobalEnd(fs,elp->y,INSERT_VALUES,gy);dCHK(err);
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
  err = dFSGlobalToExpandedBegin(fs,gx,INSERT_VALUES,elp->x);dCHK(err);
  err = dFSGlobalToExpandedEnd(fs,gx,INSERT_VALUES,elp->x);dCHK(err);
  err = VecGetArray(elp->x,&x);dCHK(err);
  err = VecZeroEntries(elp->y);dCHK(err);
  err = VecGetArray(elp->y,&y);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
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
  err = dFSRestoreWorkspace(fs,&q,&jinv,&jw,&u,&v,&du,&dv);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(elp->x,&x);dCHK(err);
  err = VecRestoreArray(elp->y,&y);dCHK(err);
  err = VecZeroEntries(gy);dCHK(err); /* Necessary? */
  err = dFSExpandedToGlobalBegin(fs,elp->y,INSERT_VALUES,gy);dCHK(err);
  err = dFSExpandedToGlobalEnd(fs,elp->y,INSERT_VALUES,gy);dCHK(err);
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
  err = dFSGlobalToExpandedBegin(fs,gx,INSERT_VALUES,elp->x);dCHK(err);
  err = dFSGlobalToExpandedEnd(fs,gx,INSERT_VALUES,elp->x);dCHK(err);
  err = VecGetArray(elp->x,&x);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,&nx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err); /* We only need space for nodal coordinates */
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
                K[ltest][lp] += basis[lq][ltest] * v[0]
                  + deriv[lq][ltest][0] * Dv[0]
                  + deriv[lq][ltest][1] * Dv[1]
                  + deriv[lq][ltest][2] * Dv[2];
              }
            }
          }
          err = dFSMatSetValuesExpanded(fs,*Jp,8,rowcol,8,rowcol,&K[0][0],ADD_VALUES);dCHK(err);
        }
      }
    }
  }
  err = dFSRestoreWorkspace(fs,&nx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(elp->x,&x);dCHK(err);

  err = MatAssemblyBegin(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  *structure = DIFFERENT_NONZERO_PATTERN;
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  Ellip elp;
  dFS fs;
  MPI_Comm comm;
  PetscViewer viewer;
  Mat J,Jp;
  Vec r,x;
  SNES snes;
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

  {
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
  err = VecDuplicate(r,&x);dCHK(err);
  err = VecZeroEntries(x);dCHK(err);
  err = VecZeroEntries(r);dCHK(err);
  err = SNESSolve(snes,NULL,x);dCHK(err);

  err = VecDestroy(r);dCHK(err);
  err = VecDestroy(x);dCHK(err);
  err = SNESDestroy(snes);dCHK(err);
  err = MatDestroy(Jp);dCHK(err);
  err = MatDestroy(J);dCHK(err);
  err = EllipDestroy(elp);dCHK(err);
  err = PetscFinalize();dCHK(err);
  return 0;
}
