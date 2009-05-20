static const char help[] = "Solve non-Newtonian Stokes problem using dual order hp elements.\n"
  "The model problem is\n"
  "  -div(eta Du) + grad(p) = f\n"
  "                  div(u) = g\n"
  "where\n"
  "  D is the symmetric gradient operator\n"
  "  eta(gamma) = A (eps^2 + gamma)^{p-2}\n"
  "  gamma = Du : Du/2\n"
  "The weak form is\n"
  "  int_Omega eta Dv:Du - p div(v) - q div(u) - f_u.v - f_p.q = 0\n"
  "with Jacobian\n"
  "  int_Omega eta Dv:Du + eta' (Dv:Dw)(Dw:Du) - p div(v) - q div(u) = 0\n"
  "The problem is linear for p=2, an incompressible for g=0\n\n";

#include "dohpfs.h"
#include "dohpvec.h"
#include "petscsnes.h"

static PetscLogEvent LOG_StokesShellMult;
typedef struct StokesCtx *Stokes;

struct StokesRheology {
  dReal A,eps,p;
};

struct StokesExactCtx {
  dReal a,b,c;
};
struct StokesExact {
  void (*solution)(const struct StokesExactCtx*,const struct StokesRheology*,const dReal x[3],dScalar u[],dScalar *p,dScalar du[],dScalar dp[]);
  void (*forcing)(const struct StokesExactCtx*,const struct StokesRheology*,const dReal x[3],dScalar fu[],dScalar *fp);
};

static void StokesExact_0_Solution(const struct StokesExactCtx *ctx,const struct StokesRheology dUNUSED *rheo,const dReal xyz[3],
                                   dScalar u[],dScalar p[],dScalar du[],dScalar dp[])
{
  const dReal dUNUSED a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2];
  u[0] = x*x*y;
  u[1] = -x*y*y;
  u[2] = 0;
  *p   = x + y - 1.0;
  /* \todo this is incorrect */
  du[0*3+0] = 0;
  du[0*3+1] = 0;
  du[0*3+2] = 0;
  du[1*3+0] = 0;
  du[1*3+1] = 0;
  du[1*3+2] = 0;
  du[2*3+0] = 0;
  du[2*3+1] = 0;
  du[2*3+2] = 0;
  dp[0]     = 0;
  dp[1]     = 0;
  dp[2]     = 0;
}
static void StokesExact_0_Forcing(const struct StokesExactCtx *ctx,const struct StokesRheology *rheo,const dReal xyz[3],dScalar fu[],dScalar *fp)
{
  const dReal dUNUSED a = ctx->a,b = ctx->b,c = ctx->c,x = xyz[0],y = xyz[1],z = xyz[2];
  fu[0] = -rheo->A*y+1;
  fu[1] = rheo->A*x+1;
  fu[2] = 0;
  *fp   = 0;
}

struct StokesStore {
  dReal eta,deta;
  dReal Du[6];
};

static dErr StokesGetNullSpace(Stokes stk,MatNullSpace *matnull);
static dErr StokesShellMatMult_All_IorA(Mat A,Vec gx,Vec gy,Vec gz,InsertMode);
static dErr StokesShellMatMult_All(Mat A,Vec gx,Vec gy)
{return StokesShellMatMult_All_IorA(A,gx,gy,NULL,INSERT_VALUES);}
static dErr StokesShellMatMultAdd_All(Mat A,Vec gx,Vec gy,Vec gz)
{return StokesShellMatMult_All_IorA(A,gx,gy,gz,ADD_VALUES);}

struct StokesCtx {
  MPI_Comm               comm;
  struct StokesRheology  rheo;
  struct StokesExact     exact;
  struct StokesExactCtx  exactctx;
  struct StokesStore    *store;
  dInt                  *storeoff;
  dJacobi                jac;
  dMesh                  mesh;
  dFS                    fsu,fsp;
  Vec                    xu,xp,yu,yp;
  Vec                    gvelocity,gvelocity_extra,gpressure,gpressure_extra,gpacked;
  VecScatter             extractVelocity,extractPressure;
  dInt                   constBDeg,nominalRDeg;
  dTruth                 errorview;
  /* Physics-based preconditioner */
  Mat Pu;                       /* preconditioner for velocity block */
  Mat Pp;                       /* preconditioner for pressure block, auxilliary system */
  Mat Au,Ap,B;                  /* shell operators */
  Mat S;                        /* Schur complement in pressure space */
  KSP kspA;                     /* Solver for Au */
  KSP kspS;                     /* Solver for S */
};

static dErr StokesCreate(MPI_Comm comm,Stokes *stokes)
{
  Stokes stk;
  dErr err;

  dFunctionBegin;
  *stokes = 0;
  err = dNew(struct StokesCtx,&stk);dCHK(err);
  stk->comm = comm;

  stk->constBDeg   = 4;
  stk->nominalRDeg = 0;
  stk->rheo.A      = 1;
  stk->rheo.eps    = 1;
  stk->rheo.p      = 2;
  *stokes = stk;
  dFunctionReturn(0);
}

static dErr MatGetVecs_Stokes(Mat A,Vec *x,Vec *y)
{
  Stokes stk;
  dInt m,n,nu,np;
  dErr err;

  dFunctionBegin;
  err = MatShellGetContext(A,(void**)&stk);dCHK(err);
  err = MatGetLocalSize(A,&m,&n);dCHK(err);
  err = VecGetLocalSize(stk->gvelocity,&nu);dCHK(err);
  err = VecGetLocalSize(stk->gpressure,&np);dCHK(err);
  if (nu==np) dERROR(1,"Degenerate case, don't know which space to copy");
  if (x) {
    if (n == nu) {
      err = VecDuplicate(stk->gvelocity,x);dCHK(err);
    } else if (n == np) {
      err = VecDuplicate(stk->gpressure,x);dCHK(err);
    } else dERROR(1,"sizes do not agree with either space");
  }
  if (y) {
    if (n == nu) {
      err = VecDuplicate(stk->gvelocity,y);dCHK(err);
    } else if (n == np) {
      err = VecDuplicate(stk->gpressure,y);dCHK(err);
    } else dERROR(1,"sizes do not agree with either space");
  }
  dFunctionReturn(0);
}

static dErr StokesSetFromOptions(Stokes stk)
{
  struct StokesRheology *rheo = &stk->rheo;
  struct StokesExactCtx *exc = &stk->exactctx;
  dMesh mesh;
  dFS fsu,fsp;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag rtag,dtag,dptag;
  dInt exact;
  dErr err;

  dFunctionBegin;
  exact = 0; exc->a = exc->b = exc->c = 1;
  err = PetscOptionsBegin(stk->comm,NULL,"Stokesicity options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-const_bdeg","Use constant isotropic degree on all elements","",stk->constBDeg,&stk->constBDeg,NULL);dCHK(err);
    stk->nominalRDeg = stk->constBDeg; /* The cheapest option, usually a good default */
    err = PetscOptionsInt("-nominal_rdeg","Nominal rule degree (will be larger if basis requires it)","",stk->nominalRDeg,&stk->nominalRDeg,NULL);dCHK(err);
    err = PetscOptionsTruth("-error_view","View errors","",stk->errorview,&stk->errorview,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_A","Rate factor (rheology)","",rheo->A,&rheo->A,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_eps","Regularization (rheology)","",rheo->eps,&rheo->eps,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_p","Power p=1+1/n where n is Glen exponent","",rheo->p,&rheo->p,NULL);dCHK(err);
    err = PetscOptionsInt("-exact","Exact solution choice","",exact,&exact,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_a","First scale parameter","",exc->a,&exc->a,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_b","Second scale parameter","",exc->b,&exc->b,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_c","Third scale parameter","",exc->c,&exc->c,NULL);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);

  switch (exact) {
    case 0:
      stk->exact.solution = StokesExact_0_Solution;
      stk->exact.forcing = StokesExact_0_Forcing;
      break;
    default: dERROR(1,"Exact solution %d not implemented");
  }

  err = dMeshCreate(stk->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);dCHK(err);
  stk->mesh = mesh;
  err = dMeshGetRoot(mesh,&domain);dCHK(err);

  err = dJacobiCreate(stk->comm,&jac);dCHK(err);
  err = dJacobiSetDegrees(jac,9,2);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);
  stk->jac = jac;

  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"stokes_rule_degree",stk->nominalRDeg,&rtag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"stokes_efs_velocity_degree",stk->constBDeg,&dtag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"stokes_efs_pressure_degree",stk->constBDeg-2,&dptag);dCHK(err);

  err = dFSCreate(stk->comm,&fsu);dCHK(err);
  err = dFSSetBlockSize(fsu,3);dCHK(err);
  err = dFSSetMesh(fsu,mesh,domain);dCHK(err);
  err = dFSSetRuleTag(fsu,jac,rtag);dCHK(err);
  err = dFSSetDegree(fsu,jac,dtag);dCHK(err);
  err = dFSRegisterBoundary(fsu,100,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSRegisterBoundary(fsu,200,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSRegisterBoundary(fsu,300,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = PetscObjectSetOptionsPrefix((dObject)fsu,"u");dCHK(err);
  err = dFSSetFromOptions(fsu);dCHK(err);
  stk->fsu = fsu;

  err = dFSCreate(stk->comm,&fsp);dCHK(err);
  err = dFSSetMesh(fsp,mesh,domain);dCHK(err);
  err = dFSSetRuleTag(fsp,jac,rtag);dCHK(err);
  err = dFSSetDegree(fsp,jac,dptag);dCHK(err);
  err = PetscObjectSetOptionsPrefix((dObject)fsp,"p");dCHK(err);
  /* No boundaries, the pressure space has Neumann conditions when Dirichlet velocity conditions are applied */
  err = dFSSetFromOptions(fsp);dCHK(err);
  stk->fsp = fsp;

  err = dFSCreateExpandedVector(fsu,&stk->xu);dCHK(err);
  err = VecDuplicate(stk->xu,&stk->yu);dCHK(err);

  err = dFSCreateExpandedVector(fsp,&stk->xp);dCHK(err);
  err = VecDuplicate(stk->xp,&stk->yp);dCHK(err);

  {                             /* Allocate space for stored values */
    dInt n,np;
    s_dRule *rule,*rulep;
    err = dFSGetElements(fsu,&n,NULL,&rule,NULL,NULL,NULL);dCHK(err);
    err = dFSGetElements(fsp,&np,NULL,&rulep,NULL,NULL,NULL);dCHK(err);
    if (n != np) dERROR(1,"pressure and velocity spaces have different number of elements");
    err = dMallocA(n+1,&stk->storeoff);dCHK(err);
    stk->storeoff[0] = 0;
    for (dInt i=0; i<n; i++) {
      dInt q,qp;
      err = dRuleGetSize(&rule[i],NULL,&q);dCHK(err);
      err = dRuleGetSize(&rulep[i],NULL,&qp);dCHK(err);
      if (q != qp) dERROR(1,"pressure and velocity spaces have different number of quadrature points on element %d",i);
      stk->storeoff[i+1] = stk->storeoff[i] + q;
    }
    err = dMallocA(stk->storeoff[n],&stk->store);dCHK(err);
    err = dMemzero(stk->store,stk->storeoff[n]*sizeof(stk->store[0]));dCHK(err);
    err = dFSRestoreElements(fsu,&n,NULL,&rule,NULL,NULL,NULL);dCHK(err);
    err = dFSRestoreElements(fsp,&np,NULL,&rulep,NULL,NULL,NULL);dCHK(err);
  }

  {
    dInt nu,np,rstart;
    IS   ublock,pblock;
    err = dFSCreateGlobalVector(stk->fsu,&stk->gvelocity);dCHK(err);
    err = VecDuplicate(stk->gvelocity,&stk->gvelocity_extra);dCHK(err);
    err = dFSCreateGlobalVector(stk->fsp,&stk->gpressure);dCHK(err);
    err = VecDuplicate(stk->gpressure,&stk->gpressure_extra);dCHK(err);
    err = VecGetLocalSize(stk->gvelocity,&nu);dCHK(err);
    err = VecGetLocalSize(stk->gpressure,&np);dCHK(err);
    err = VecCreateMPI(stk->comm,nu+np,PETSC_DETERMINE,&stk->gpacked);dCHK(err);
    err = VecGetOwnershipRange(stk->gpacked,&rstart,NULL);dCHK(err);
    err = ISCreateStride(stk->comm,nu,rstart,1,&ublock);dCHK(err);
    err = ISCreateStride(stk->comm,np,rstart+nu,1,&pblock);dCHK(err);
    err = VecScatterCreate(stk->gpacked,ublock,stk->gvelocity,NULL,&stk->extractVelocity);dCHK(err);
    err = VecScatterCreate(stk->gpacked,pblock,stk->gpressure,NULL,&stk->extractPressure);dCHK(err);
    err = ISDestroy(ublock);dCHK(err);
    err = ISDestroy(pblock);dCHK(err);
  }

  {                             /* Create shell matrices for use by the preconditioner */
    dInt nu,np;
    err = VecGetLocalSize(stk->gvelocity,&nu);dCHK(err);
    err = VecGetLocalSize(stk->gpressure,&np);dCHK(err);
    err = dFSGetMatrix(stk->fsu,MATSHELL,&stk->Au);dCHK(err);
    err = MatShellSetContext(stk->Au,stk);dCHK(err);
    err = MatShellSetOperation(stk->Au,MATOP_GET_VECS,(void(*)(void))MatGetVecs_Stokes);dCHK(err);
    err = MatShellSetOperation(stk->Au,MATOP_MULT,(void(*)(void))StokesShellMatMult_All);dCHK(err);
    err = MatShellSetOperation(stk->Au,MATOP_MULT_TRANSPOSE,(void(*)(void))StokesShellMatMult_All);dCHK(err); /* matrix is symmetric */
    err = MatShellSetOperation(stk->Au,MATOP_MULT_ADD,(void(*)(void))StokesShellMatMultAdd_All);dCHK(err);
    err = MatShellSetOperation(stk->Au,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))StokesShellMatMultAdd_All);dCHK(err);
    err = MatCreateShell(stk->comm,np,nu,PETSC_DETERMINE,PETSC_DETERMINE,stk,&stk->B);dCHK(err);
    err = MatShellSetOperation(stk->B,MATOP_GET_VECS,(void(*)(void))MatGetVecs_Stokes);dCHK(err);
    err = MatShellSetOperation(stk->B,MATOP_MULT,(void(*)(void))StokesShellMatMult_All);dCHK(err);
    err = MatShellSetOperation(stk->B,MATOP_MULT_TRANSPOSE,(void(*)(void))StokesShellMatMult_All);dCHK(err);
    err = MatShellSetOperation(stk->B,MATOP_MULT_ADD,(void(*)(void))StokesShellMatMultAdd_All);dCHK(err);
    err = MatShellSetOperation(stk->B,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))StokesShellMatMultAdd_All);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr StokesDestroy(Stokes stk)
{
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(stk->fsu);dCHK(err);
  err = dFSDestroy(stk->fsp);dCHK(err);
  err = dJacobiDestroy(stk->jac);dCHK(err);
  err = dMeshDestroy(stk->mesh);dCHK(err);
  err = dFree(stk->storeoff);dCHK(err);
  err = dFree(stk->store);dCHK(err);
#define _D(v)  do {if (v) {err = VecDestroy(v);dCHK(err);}} while (0)
  _D(stk->xu);
  _D(stk->yu);
  _D(stk->xp);
  _D(stk->yp);
  _D(stk->gvelocity);
  _D(stk->gpressure);
  _D(stk->gvelocity_extra);
  _D(stk->gpressure_extra);
  _D(stk->gpacked);
#undef _D
  err = VecScatterDestroy(stk->extractVelocity);dCHK(err);
  err = VecScatterDestroy(stk->extractPressure);dCHK(err);
  err = MatDestroy(stk->Au);dCHK(err);
  err = MatDestroy(stk->B);dCHK(err);
  if (stk->Pu) {err = MatDestroy(stk->Pu);dCHK(err);}
  if (stk->S)  {err = MatDestroy(stk->S);dCHK(err);}
  if (stk->kspA) {err = KSPDestroy(stk->kspA);dCHK(err);}
  if (stk->kspS) {err = KSPDestroy(stk->kspS);dCHK(err);}
  err = dFree(stk);dCHK(err);
  dFunctionReturn(0);
}

static inline void StokesPointwiseComputeStore(struct StokesRheology dUNUSED *rheo,const dReal dUNUSED x[3],const dScalar Du[],struct StokesStore *st)
{
  dScalar gamma_reg = dSqr(rheo->eps) + dColonSymScalar3(Du,Du);
  st->eta = rheo->A * pow(gamma_reg,rheo->p-2);
  st->deta = (rheo->p-2) * st->eta / gamma_reg;
  for (dInt i=0; i<6; i++) st->Du[i] = Du[i];
}

static inline void StokesPointwiseFunction(struct StokesRheology *rheo,struct StokesExact *exact,struct StokesExactCtx *exactctx,
                                           const dReal x[3],dReal weight,const dScalar Du[6],dScalar p,
                                           struct StokesStore *st,dScalar v[3],dScalar Dv[6],dScalar *q)
{
  dScalar fu[3],fp;
  StokesPointwiseComputeStore(rheo,x,Du,st);
  exact->forcing(exactctx,rheo,x,fu,&fp);
  for (dInt i=0; i<3; i++) v[i] = -weight * fu[i]; /* Coefficient of \a v in weak form, only appears in forcing term */
  *q   = -weight * (Du[0]+Du[1]+Du[2] + fp);       /* -q tr(Du) - forcing, note tr(Du) = div(u) */
  for (dInt i=0; i<3; i++) Dv[i] = weight * (st->eta * Du[i] - p); /* eta Dv:Du - p tr(Dv) */
  for (dInt i=3; i<6; i++) Dv[i] = weight * st->eta * Du[i];       /* eta Dv:Du */
}

static inline void StokesPointwiseJacobian(const struct StokesStore *restrict st,dReal weight,
                                           const dScalar Du[restrict static 6],dScalar p,
                                           dScalar Dv[restrict static 6],dScalar *restrict q)
{
                                /* Coefficients in weak form of Jacobian */
  const dScalar deta_colon = st->deta*dColonSymScalar3(st->Du,Du);                          /* eta' Dw:Du */
  for (dInt i=0; i<3; i++) Dv[i] = weight * (st->eta * Du[i] + deta_colon * st->Du[i] - p); /* eta Dv:Du + eta' (Dv:Dw)(Dw:Du) - p tr(Dv) */
  for (dInt i=3; i<6; i++) Dv[i] = weight * (st->eta * Du[i] + deta_colon * st->Du[i]);     /* eta Dv:Du + eta' (Dv:Dw)(Dw:Du) */
  *q = -weight*(Du[0]+Du[1]+Du[2]);                                                         /* -q tr(Du) */
}

static inline void StokesPointwiseJacobian_A(const struct StokesStore *restrict st,dReal weight,const dScalar Du[restrict static 6],dScalar Dv[restrict static 6])
{
  const dScalar deta_colon = st->deta*dColonSymScalar3(st->Du,Du);
  for (dInt i=0; i<6; i++) Dv[i] = weight * (st->eta*Du[i] + deta_colon*st->Du[i]);
}

static inline void StokesPointwiseJacobian_B(dReal weight,const dScalar Du[restrict static 6],dScalar *restrict q)
{
  *q = -weight*(Du[0]+Du[1]+Du[2]);
}

static inline void StokesPointwiseJacobian_Bt(dReal weight,dScalar p,dScalar Dv[restrict static 6])
{
  for (dInt i=0; i<3; i++) Dv[i] = -weight*p;
  for (dInt i=3; i<6; i++) Dv[i] = 0;
}

static dErr StokesFunction(SNES dUNUSED snes,Vec gx,Vec gy,void *ctx)
{
  dReal (*restrict geom)[3],(*restrict q)[3],(*restrict jinv)[3][3],*restrict jw;
  Stokes   stk = ctx;
  dFS      fsu = stk->fsu,fsp = stk->fsp;
  Vec      gxu,gxp;
  dInt     n,np,*off,*offp,*geomoff;
  s_dRule *rule,*rulep;
  s_dEFS  *efs,*efsp;
  dScalar *xu,*xp,*yu,*yp,*restrict vv,*restrict du,*restrict dv,*restrict pp,*restrict qq;
  dErr    err;

  dFunctionBegin;
  gxu = stk->gvelocity; gxp = stk->gpressure; /* Our work vectors */
  err = VecScatterBegin(stk->extractVelocity,gx,gxu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd  (stk->extractVelocity,gx,gxu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterBegin(stk->extractPressure,gx,gxp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd  (stk->extractPressure,gx,gxp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  /* dFS_INHOMOGENEOUS projects into inhomogeneous space (strongly enforcing boundary conditions) */
  err = dFSGlobalToExpanded(fsu,gxu,stk->xu,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err); /* velocity */
  err = dFSGlobalToExpanded(fsp,gxp,stk->xp,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err); /* pressure */
  err = VecGetArray(stk->xu,&xu);dCHK(err);
  err = VecGetArray(stk->xp,&xp);dCHK(err);
  err = VecGetArray(stk->yu,&yu);dCHK(err);
  err = VecGetArray(stk->yp,&yp);dCHK(err);
  err = dFSGetElements(fsu,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err); /* \note \a off is in terms of \e nodes, not \e dofs */
  err = dFSGetElements(fsp,&np,&offp,&rulep,&efsp,NULL,NULL);dCHK(err);
  if (n != np) dERROR(1,"number of elements in velocity and pressure spaces do not agree");
  err = dFSGetWorkspace(fsu,__func__,&q,&jinv,&jw,NULL,&vv,&du,&dv);dCHK(err);
  err = dFSGetWorkspace(fsp,__func__,NULL,NULL,NULL,&pp,&qq,NULL,NULL);dCHK(err); /* workspace for test and trial functions */
  for (dInt e=0; e<n; e++) {
    dInt Q;
    err = dRuleComputeGeometry(&rule[e],(const dReal(*)[3])(geom+geomoff[e]),q,jinv,jw);dCHK(err);
    err = dRuleGetSize(&rule[e],0,&Q);dCHK(err);
    {
      dInt Qp;
      err = dRuleGetSize(&rulep[e],0,&Qp);dCHK(err);
      if (Q != Qp) dERROR(1,"rule sizes on element %d do not agree",e);dCHK(err);
    }
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,xu+3*off[e],du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err); /* velocity gradients */
    err = dEFSApply(&efsp[e],(const dReal*)jinv,1,xp+offp[e],pp,dAPPLY_INTERP,INSERT_VALUES);dCHK(err); /* pressure values */
    for (dInt i=0; i<Q; i++) {
      struct StokesStore *restrict st = &stk->store[stk->storeoff[e]+i];
      dScalar Du[6],Dv[6];
      dTensorSymCompress3(&du[i*9],Du);
      StokesPointwiseFunction(&stk->rheo,&stk->exact,&stk->exactctx,q[i],jw[i],Du,pp[i],st,&vv[i*3],Dv,&qq[i]);
      dTensorSymUncompress3(Dv,&dv[i*9]);
    }
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,vv,yu+3*off[e],dAPPLY_INTERP_TRANSPOSE,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,dv,yu+3*off[e],dAPPLY_GRAD_TRANSPOSE,ADD_VALUES);dCHK(err);
    err = dEFSApply(&efsp[e],(const dReal*)jinv,1,qq,yp+offp[e],dAPPLY_INTERP_TRANSPOSE,INSERT_VALUES);dCHK(err);
  }
  err = dFSRestoreWorkspace(fsu,__func__,&q,&jinv,&jw,NULL,&vv,&du,&dv);dCHK(err);
  err = dFSRestoreWorkspace(fsp,__func__,NULL,NULL,NULL,&pp,&qq,NULL,NULL);dCHK(err);
  err = dFSRestoreElements(fsu,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSRestoreElements(fsp,&np,&offp,&rulep,&efsp,NULL,NULL);dCHK(err);
  err = VecRestoreArray(stk->xu,&xu);dCHK(err);
  err = VecRestoreArray(stk->xp,&xp);dCHK(err);
  err = VecRestoreArray(stk->yu,&yu);dCHK(err);
  err = VecRestoreArray(stk->yp,&yp);dCHK(err);
  {
    err = VecZeroEntries(gxu);dCHK(err);
    err = VecZeroEntries(gxp);dCHK(err);
    err = dFSExpandedToGlobal(fsu,stk->yu,gxu,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
    err = dFSExpandedToGlobal(fsp,stk->yp,gxp,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
    err = VecScatterBegin(stk->extractVelocity,gxu,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
    err = VecScatterEnd  (stk->extractVelocity,gxu,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
    err = VecScatterBegin(stk->extractPressure,gxp,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
    err = VecScatterEnd  (stk->extractPressure,gxp,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr StokesShellMatMult_J(Mat J,Vec gx,Vec gy)
{
  dReal (*restrict geom)[3],(*restrict q)[3],(*restrict jinv)[3][3],*restrict jw;
  Stokes   stk;
  dFS      fsu,fsp;
  Vec      gxu,gxp;
  dInt     n,np,*off,*offp,*geomoff;
  s_dRule *rule,*rulep;
  s_dEFS  *efs,*efsp;
  dScalar *xu,*xp,*yu,*yp,*restrict vv,*restrict du,*restrict dv,*restrict pp,*restrict qq;
  dErr    err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_StokesShellMult,J,gx,gy,0);dCHK(err);
  err = MatShellGetContext(J,(void**)&stk);dCHK(err);
  fsu = stk->fsu; fsp = stk->fsp;
  gxu = stk->gvelocity; gxp = stk->gpressure; /* Our work vectors */
  err = VecScatterBegin(stk->extractVelocity,gx,gxu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd  (stk->extractVelocity,gx,gxu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterBegin(stk->extractPressure,gx,gxp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd  (stk->extractPressure,gx,gxp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  /* dFS_HOMOGENEOUS projects into homogeneous space (because Dirichlet conditions are enforced strongly) */
  err = dFSGlobalToExpanded(fsu,gxu,stk->xu,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err); /* velocity */
  err = dFSGlobalToExpanded(fsp,gxp,stk->xp,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err); /* pressure */
  err = VecGetArray(stk->xu,&xu);dCHK(err);
  err = VecGetArray(stk->xp,&xp);dCHK(err);
  err = VecGetArray(stk->yu,&yu);dCHK(err);
  err = VecGetArray(stk->yp,&yp);dCHK(err);
  err = dFSGetElements(fsu,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err); /* \note \a off is in terms of \e nodes, not \e dofs */
  err = dFSGetElements(fsp,&np,&offp,&rulep,&efsp,NULL,NULL);dCHK(err);
  if (n != np) dERROR(1,"number of elements in velocity and pressure spaces do not agree");
  err = dFSGetWorkspace(fsu,__func__,&q,&jinv,&jw,NULL,&vv,&du,&dv);dCHK(err);
  err = dFSGetWorkspace(fsp,__func__,NULL,NULL,NULL,&pp,&qq,NULL,NULL);dCHK(err); /* workspace for test and trial functions */
  for (dInt e=0; e<n; e++) {
    dInt Q;
    err = dRuleComputeGeometry(&rule[e],(const dReal(*)[3])(geom+geomoff[e]),q,jinv,jw);dCHK(err);
    err = dRuleGetSize(&rule[e],0,&Q);dCHK(err);
    {
      dInt Qp;
      err = dRuleGetSize(&rulep[e],0,&Qp);dCHK(err);
      if (Q != Qp) dERROR(1,"rule sizes on element %d do not agree",e);dCHK(err);
    }
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,xu+3*off[e],du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err); /* velocity gradients */
    err = dEFSApply(&efsp[e],(const dReal*)jinv,1,xp+offp[e],pp,dAPPLY_INTERP,INSERT_VALUES);dCHK(err); /* pressure values */
    for (dInt i=0; i<Q; i++) {
      struct StokesStore *restrict st = &stk->store[stk->storeoff[e]+i];
      dScalar Du[6],Dv[6];
      dTensorSymCompress3(&du[i*9],Du);
      StokesPointwiseJacobian(st,jw[i],Du,pp[i],Dv,&qq[i]);
      dTensorSymUncompress3(Dv,&dv[i*9]);
    }
    err = dEFSApply(&efs[e],(const dReal*)jinv,3,dv,yu+3*off[e],dAPPLY_GRAD_TRANSPOSE,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efsp[e],(const dReal*)jinv,1,qq,yp+offp[e],dAPPLY_INTERP_TRANSPOSE,INSERT_VALUES);dCHK(err);
  }
  err = dFSRestoreWorkspace(fsu,__func__,&q,&jinv,&jw,NULL,&vv,&du,&dv);dCHK(err);
  err = dFSRestoreWorkspace(fsp,__func__,NULL,NULL,NULL,&pp,&qq,NULL,NULL);dCHK(err);
  err = dFSRestoreElements(fsu,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSRestoreElements(fsp,&np,&offp,&rulep,&efsp,NULL,NULL);dCHK(err);
  err = VecRestoreArray(stk->xu,&xu);dCHK(err);
  err = VecRestoreArray(stk->xp,&xp);dCHK(err);
  err = VecRestoreArray(stk->yu,&yu);dCHK(err);
  err = VecRestoreArray(stk->yp,&yp);dCHK(err);
  {
    err = VecZeroEntries(gxu);dCHK(err);
    err = VecZeroEntries(gxp);dCHK(err);
    err = dFSExpandedToGlobal(fsu,stk->yu,gxu,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err);
    err = dFSExpandedToGlobal(fsp,stk->yp,gxp,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err);
    err = VecScatterBegin(stk->extractVelocity,gxu,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
    err = VecScatterEnd  (stk->extractVelocity,gxu,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
    err = VecScatterBegin(stk->extractPressure,gxp,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
    err = VecScatterEnd  (stk->extractPressure,gxp,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  }
  err = PetscLogEventEnd(LOG_StokesShellMult,J,gx,gy,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesShellMatMult_J_block(Mat J,Vec gx,Vec gy)
{
  Stokes stk;
  Vec    gxu,gxp,gyu;
  dErr   err;

  dFunctionBegin;
  err = MatShellGetContext(J,(void**)&stk);dCHK(err);
  gxu = stk->gvelocity; gxp = stk->gpressure;
  gyu = stk->gvelocity_extra;
  err = VecScatterBegin(stk->extractVelocity,gx,gxu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd  (stk->extractVelocity,gx,gxu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = MatMult(stk->B,gxu,gxp);dCHK(err); /* Use pressure vector for output of p = B u */
  err = VecScatterBegin(stk->extractPressure,gxp,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (stk->extractPressure,gxp,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);

  err = MatMult(stk->Au,gxu,gyu);dCHK(err);
  err = VecScatterBegin(stk->extractPressure,gx,gxp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd  (stk->extractPressure,gx,gxp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = MatMultTransposeAdd(stk->B,gxp,gyu,gyu);dCHK(err);
  err = VecScatterBegin(stk->extractVelocity,gyu,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (stk->extractVelocity,gyu,gy,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  dFunctionReturn(0);
}

typedef enum {STOKES_MULT_A,STOKES_MULT_Bt,STOKES_MULT_B} StokesMultMode;

static dErr dUNUSED StokesShellMatMult_All_IorA(Mat A,Vec gx,Vec gy,Vec gz,InsertMode imode)
{
  Stokes          stk;
  dFS             fsx,fsy,fslarger;
  dInt            n,*off,*offy,*geomoff;
  s_dRule        *rule;
  s_dEFS         *efs,*efsy;
  Vec             X,Y;
  dScalar        *x,*y,*restrict uu,*restrict vv,*restrict du,*restrict dv;
  StokesMultMode  mmode;
  InsertMode      apply_imode;
  dReal (*restrict geom)[3],(*restrict q)[3],(*restrict jinv)[3][3],*restrict jw;
  dErr err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_StokesShellMult,A,gx,gy,gz);dCHK(err);
  err = MatShellGetContext(A,(void**)&stk);dCHK(err);
  {  /* Find out which block we have by comparing sizes */
    dInt nu,np,nx,ny;
    err = VecGetSize(stk->gvelocity,&nu);dCHK(err);
    err = VecGetSize(stk->gpressure,&np);dCHK(err);
    err = VecGetSize(gx,&nx);dCHK(err);
    err = VecGetSize(gy,&ny);dCHK(err);
    if (nx==nu && ny==nu) mmode = STOKES_MULT_A;
    else if (nx==np && ny==nu) mmode = STOKES_MULT_Bt;
    else if (nx==nu && ny==np) mmode = STOKES_MULT_B;
    else dERROR(1,"Sizes do not match, unknown mult operation");
  }
  switch (mmode) {
    case STOKES_MULT_A:  fsx = fsy = stk->fsu; X = stk->xu; Y = stk->yu; break;
    case STOKES_MULT_Bt: fsx = stk->fsp; fsy = stk->fsu; X = stk->xp; Y = stk->yu; break;
    case STOKES_MULT_B:  fsx = stk->fsu; fsy = stk->fsp; X = stk->xu; Y = stk->yp; break;
    default: dERROR(1,"should not happen");
  }
  apply_imode = (imode == ADD_VALUES && gy != gz) ? ADD_VALUES : INSERT_VALUES; /* \bug */
  err = dFSGlobalToExpanded(fsx,gx,X,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(X,&x);dCHK(err);
  err = VecGetArray(Y,&y);dCHK(err);
  err = dFSGetElements(fsx,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetElements(fsy,NULL,&offy,NULL,&efsy,NULL,NULL);dCHK(err);
  fslarger = off[n]>offy[n] ? fsx : fsy;
  err = dFSGetWorkspace(fslarger,__func__,&q,&jinv,&jw,&uu,&vv,&du,&dv);dCHK(err);
  for (dInt e=0; e<n; e++) {
    dInt Q;
    err = dRuleComputeGeometry(&rule[e],(const dReal(*)[3])(geom+geomoff[e]),q,jinv,jw);dCHK(err);
    err = dRuleGetSize(&rule[e],0,&Q);dCHK(err);
    switch (mmode) {
      case STOKES_MULT_A:
        err = dEFSApply(&efs[e],(const dReal*)jinv,3,x+3*off[e],du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
        for (dInt i=0; i<Q; i++) {
          struct StokesStore *restrict st = &stk->store[stk->storeoff[e]+i];
          dScalar Du[6],Dv[6],qq_unused[1];
          dTensorSymCompress3(&du[i*9],Du);
          StokesPointwiseJacobian(st,jw[i],Du,0,Dv,qq_unused);
          //StokesPointwiseJacobian_A(st,jw[i],Du,Dv);
          dTensorSymUncompress3(Dv,&dv[i*9]);
        }
        err = dEFSApply(&efsy[e],(const dReal*)jinv,3,dv,y+3*offy[e],dAPPLY_GRAD_TRANSPOSE,INSERT_VALUES);dCHK(err);
        break;
      case STOKES_MULT_Bt:
        err = dEFSApply(&efs[e],(const dReal*)jinv,1,x+off[e],uu,dAPPLY_INTERP,INSERT_VALUES);dCHK(err); /* pressure values */
        for (dInt i=0; i<Q; i++) {
          dScalar Dv[6];
          StokesPointwiseJacobian_Bt(jw[i],uu[i],Dv);
          dTensorSymUncompress3(Dv,&dv[i*9]);
        }
        err = dEFSApply(&efsy[e],(const dReal*)jinv,3,dv,y+3*offy[e],dAPPLY_GRAD_TRANSPOSE,INSERT_VALUES);dCHK(err);
        break;
      case STOKES_MULT_B:
        err = dEFSApply(&efs[e],(const dReal*)jinv,3,x+3*off[e],du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
        for (dInt i=0; i<Q; i++) {
          dScalar Du[6];
          dTensorSymCompress3(&du[i*9],Du);
          StokesPointwiseJacobian_B(jw[i],Du,&vv[i]); /* vv is pressure test function */
        }
        err = dEFSApply(&efsy[e],(const dReal*)jinv,1,vv,y+offy[e],dAPPLY_INTERP_TRANSPOSE,INSERT_VALUES);dCHK(err);
        break;
    }
  }
  err = dFSRestoreWorkspace(fslarger,__func__,&q,&jinv,&jw,&uu,&vv,&du,&dv);dCHK(err);
  err = dFSRestoreElements(fsx,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSRestoreElements(fsy,NULL,&offy,NULL,&efsy,NULL,NULL);dCHK(err);
  err = VecRestoreArray(X,&x);dCHK(err);
  err = VecRestoreArray(Y,&y);dCHK(err);
  switch (imode) {
    case INSERT_VALUES:
      if (gz) dERROR(1,"Cannot use INSERT_VALUES and set gz");
      gz = gy;
      err = VecZeroEntries(gz);dCHK(err);
      break;
    case ADD_VALUES:
      if (gz != gy) {
        err = VecCopy(gy,gz);dCHK(err);
      }
      break;
    default: dERROR(1,"unsupported imode");
  }
  err = dFSExpandedToGlobal(fsy,Y,gz,dFS_HOMOGENEOUS,imode);dCHK(err);
  err = PetscLogEventEnd(LOG_StokesShellMult,A,gx,gy,gz);dCHK(err);
  dFunctionReturn(0);
}

#if defined(ENABLE_PRECONDITIONING)
static dErr StokesJacobian(SNES dUNUSED snes,Vec gx,Mat *J,Mat *Jp,MatStructure *structure,void *ctx)
{
  Stokes stk = ctx;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*nx)[3];
  dScalar *x;
  dFS fs = stk->fs;
  dInt n,*off,*geomoff;
  dReal (*geom)[3];
  dErr err;

  dFunctionBegin;
  err = MatZeroEntries(*Jp);dCHK(err);
  err = dFSGlobalToExpanded(fs,gx,stk->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(stk->x,&x);dCHK(err);
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
            struct StokesStore st;
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
              StokesPointwiseComputeStore(&stk->rheo,qx[lq],st_u,&st_Du[0][0],&st);
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
                  StokesPointwiseJacobian(&stk->rheo,&st,jw[lq],u,&Du[0][0],v,&Dv[0][0]);
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
  err = VecRestoreArray(stk->x,&x);dCHK(err);

  /* These are both shell matrices, we call this so SNES knows the matrices have changed */
  err = MatAssemblyBegin(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  *structure = SAME_NONZERO_PATTERN;
  dFunctionReturn(0);
}

static dErr StokesErrorNorms(Stokes stk,Vec gx,dReal errorNorms[static 3],dReal gerrorNorms[static 3])
{
  dFS fs = stk->fs;
  dInt n,*off,*geomoff;
  s_dRule *rule;
  s_dEFS *efs;
  dReal (*geom)[3],(*q)[3],(*jinv)[3][3],*jw;
  dScalar *x,(*u)[3],(*du)[9];
  dErr err;

  dFunctionBegin;
  err = dMemzero(errorNorms,3*sizeof(errorNorms));dCHK(err);
  err = dMemzero(gerrorNorms,3*sizeof(gerrorNorms));dCHK(err);
  err = dFSGlobalToExpanded(fs,gx,stk->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(stk->x,&x);dCHK(err);
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
      stk->exact.solution(&stk->exactctx,&stk->rheo,q[i],uu,duu);
      for (dInt j=0; j<3; j++) {
        r[j] = u[i][j] - uu[j]; /* Function error at point */
        rsum += dSqr(r[j]);
        gr[j] = dSqrt(dSqr(du[i][j*3+0]-duu[j*3+0]) + dSqr(du[i][j*3+1]-duu[j*3+1]) + dSqr(du[i][j*3+2]-duu[j*3+2])); /* Gradient error at point */
        grsum += dSqr(gr[j]);
      }
      if (stk->errorview) {
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
  err = VecRestoreArray(stk->x,&x);dCHK(err);
  errorNorms[1] = dSqrt(errorNorms[1]);
  gerrorNorms[1] = dSqrt(gerrorNorms[1]);
  dFunctionReturn(0);
}
#endif /* defined(ENABLE_PRECONDITIONING) */

static dErr StokesPCSetUp(void *ctx)
{
  Stokes stk = ctx;
  PC pc;
  dErr err;

  dFunctionBegin;
  if (!stk->kspA) {
    err = KSPCreate(stk->comm,&stk->kspA);dCHK(err);
    err = KSPGetPC(stk->kspA,&pc);dCHK(err);
    err = PCSetType(pc,PCNONE);dCHK(err);
    err = KSPSetOptionsPrefix(stk->kspA,"saddle_A_");dCHK(err);
    err = KSPSetFromOptions(stk->kspA);dCHK(err);
  }
  if (1) {
    if (stk->Pu) {err = MatDestroy(stk->Pu);dCHK(err);}
    err = MatComputeExplicitOperator(stk->Au,&stk->Pu);dCHK(err);
  }
  err = KSPSetOperators(stk->kspA,stk->Au,stk->Pu?stk->Pu:stk->Au,SAME_NONZERO_PATTERN);dCHK(err);
  if (!stk->S) {
    Mat Bt;
    KSP ksp;
    MatNullSpace matnull;
    err = MatCreateTranspose(stk->B,&Bt);dCHK(err);
    err = MatCreateSchurComplement(stk->Au,Bt,stk->B,NULL,&stk->S);dCHK(err);
    err = MatSchurComplementGetKSP(stk->S,&ksp);dCHK(err);
    err = MatNullSpaceCreate(stk->comm,dTRUE,0,NULL,&matnull);dCHK(err);
    err = KSPSetNullSpace(ksp,matnull);dCHK(err);
    err = MatSetFromOptions(stk->S);dCHK(err);
    err = MatNullSpaceDestroy(matnull);dCHK(err);
    err = MatDestroy(Bt);dCHK(err);
  }
  if (!stk->kspS) {
    err = KSPCreate(stk->comm,&stk->kspS);dCHK(err);
    err = KSPGetPC(stk->kspS,&pc);dCHK(err);
    err = PCSetType(pc,PCNONE);dCHK(err);
    err = KSPSetOptionsPrefix(stk->kspS,"saddle_S_");dCHK(err);
    err = KSPSetFromOptions(stk->kspS);dCHK(err);
  }
  err = KSPSetOperators(stk->kspS,stk->S,stk->Pp?stk->Pp:stk->S,SAME_NONZERO_PATTERN);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesPCApply(void *ctx,Vec x,Vec y)
{
  Stokes stk = ctx;
  Vec xu,xp,yu,yp;
  dErr err;

  dFunctionBegin;
  xu = stk->gvelocity;
  xp = stk->gpressure;
  yu = stk->gvelocity_extra;
  yp = stk->gpressure_extra;
  err = VecScatterBegin(stk->extractPressure,x,xp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd  (stk->extractPressure,x,xp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = KSPSolve(stk->kspS,xp,yp);dCHK(err);
  err = VecScatterBegin(stk->extractPressure,yp,y,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (stk->extractPressure,yp,y,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = MatMultTranspose(stk->B,xp,xu);dCHK(err);
  err = VecScale(xu,-1);dCHK(err);
  err = VecScatterBegin(stk->extractVelocity,x,xu,ADD_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterEnd  (stk->extractVelocity,x,xu,ADD_VALUES,SCATTER_FORWARD);dCHK(err);
  err = KSPSolve(stk->kspA,xu,yu);dCHK(err);
  err = VecScatterBegin(stk->extractVelocity,yu,y,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (stk->extractVelocity,yu,y,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesGetSolutionField_All(Stokes stk,dFS fs,dTruth isvel,Vec *insoln)
{
  Vec      sol,xc,cvec;
  dScalar *x,*coords;
  dInt     n,bs;
  dErr     err;

  dFunctionBegin;
  *insoln = 0;
  err = dFSCreateGlobalVector(fs,&sol);dCHK(err);
  err = VecDohpGetClosure(sol,&xc);dCHK(err);
  err = dFSGetCoordinates(fs,&cvec);dCHK(err);
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
    dScalar u_unused[3],p_unused[1],du_unused[3*3],dp_unused[3];
    /* if \a isvel then \a x is the velocity field, otherwise it is the pressure field */
    stk->exact.solution(&stk->exactctx,&stk->rheo,&coords[3*i],isvel ? &x[i*bs] : u_unused,isvel ? p_unused : &x[i*bs],du_unused,dp_unused);
    /* printf("Node %3d: coords %+8f %+8f %+8f   exact %+8f %+8f %+8f\n",i,coords[3*i],coords[3*i+1],coords[3*i+2],x[3*i],x[3*i+1],x[3*i+2]); */
  }
  err = VecRestoreArray(xc,&x);dCHK(err);
  err = VecRestoreArray(cvec,&coords);dCHK(err);
  err = VecDestroy(cvec);dCHK(err);
  err = dFSInhomogeneousDirichletCommit(fs,xc);dCHK(err);
  err = VecDohpRestoreClosure(sol,&xc);dCHK(err);
  *insoln = sol;
  dFunctionReturn(0);
}

/** Creates a solution vector, commits the closure to each FS, returns packed solution vector */
static dErr StokesGetSolutionVector(Stokes stk,Vec *insoln)
{
  dErr err;
  Vec solu,solp,spacked;

  dFunctionBegin;
  *insoln = 0;
  err = StokesGetSolutionField_All(stk,stk->fsu,dTRUE,&solu);dCHK(err);
  err = StokesGetSolutionField_All(stk,stk->fsp,dFALSE,&solp);dCHK(err);
  err = VecDuplicate(stk->gpacked,&spacked);dCHK(err);
  err = VecScatterBegin(stk->extractVelocity,solu,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (stk->extractVelocity,solu,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(stk->extractPressure,solp,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (stk->extractPressure,solp,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecDestroy(solu);dCHK(err);
  err = VecDestroy(solp);dCHK(err);
  *insoln = spacked;
  dFunctionReturn(0);
}

static dErr StokesGetNullSpace(Stokes stk,MatNullSpace *matnull)
{
  dErr err;
  Vec r;

  dFunctionBegin;
  err = VecDuplicate(stk->gpacked,&r);dCHK(err);
  err = VecZeroEntries(r);dCHK(err);
  err = VecSet(stk->gpressure,1);dCHK(err);
  err = VecScatterBegin(stk->extractPressure,stk->gpressure,r,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (stk->extractPressure,stk->gpressure,r,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecNormalize(r,PETSC_NULL);dCHK(err);
  err = MatNullSpaceCreate(stk->comm,dFALSE,1,&r,matnull);dCHK(err);
  err = VecDestroy(r);dCHK(err);
  dFunctionReturn(0);
}


int main(int argc,char *argv[])
{
  char mtype[256] = MATBAIJ;
  Stokes stk;
  MPI_Comm comm;
  PetscViewer viewer;
  Mat J,Jp;
  MatFDColoring fdcolor = 0;
  Vec r,x,soln;
  SNES snes;
  dTruth nocheck,compute_explicit,use_jblock;
  dErr err;

  err = PetscInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = PetscLogEventRegister("StokesShellMult",MAT_COOKIE,&LOG_StokesShellMult);dCHK(err);

  err = StokesCreate(comm,&stk);dCHK(err);
  err = StokesSetFromOptions(stk);dCHK(err);

  err = VecDuplicate(stk->gpacked,&r);dCHK(err);

  err = PetscOptionsGetString(NULL,"-q1mat_type",mtype,sizeof(mtype),NULL);dCHK(err);
  if (0) {
    err = dFSGetMatrix(stk->fsu,mtype,&Jp);dCHK(err);
    err = MatSetOptionsPrefix(Jp,"q1");dCHK(err);
    err = MatSeqAIJSetPreallocation(Jp,27,NULL);dCHK(err);
  }

  err = PetscOptionsBegin(stk->comm,NULL,"Stokes solver options",__FILE__);dCHK(err); {
    err = PetscOptionsName("-nocheck_error","Do not compute errors","",&nocheck);dCHK(err);
    err = PetscOptionsName("-compute_explicit","Compute explicit Jacobian (only very small sizes)","",&compute_explicit);dCHK(err);
    err = PetscOptionsName("-use_jblock","Use blocks to apply Jacobian instead of unified (more efficient) version","",&use_jblock);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);
  {
    dInt m;
    err = VecGetLocalSize(r,&m);dCHK(err);
    err = MatCreateShell(comm,m,m,PETSC_DETERMINE,PETSC_DETERMINE,stk,&J);dCHK(err);
    err = MatSetOptionsPrefix(J,"j");dCHK(err);
    if (use_jblock) {
      err = MatShellSetOperation(J,MATOP_MULT,(void(*)(void))StokesShellMatMult_J_block);dCHK(err);
    } else {
      err = MatShellSetOperation(J,MATOP_MULT,(void(*)(void))StokesShellMatMult_J);dCHK(err);
    }
    err = MatCreate(comm,&Jp);dCHK(err);
    err = MatSetSizes(Jp,m,m,PETSC_DETERMINE,PETSC_DETERMINE);dCHK(err);
    err = MatSetOptionsPrefix(Jp,"jp");dCHK(err);
    err = MatSetFromOptions(Jp);dCHK(err);
    for (dInt i=0; i<m; i++) {
      err = MatSetValue(Jp,i,i,0,INSERT_VALUES);dCHK(err);
    }
    err = MatAssemblyBegin(Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd(Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  }
  err = SNESCreate(comm,&snes);dCHK(err);
  err = SNESSetFunction(snes,r,StokesFunction,stk);dCHK(err);
  switch (2) {
    case 1:
      err = SNESSetJacobian(snes,J,Jp,SNESDefaultComputeJacobian,stk);dCHK(err); break;
    case 2: {
      ISColoring iscolor;
      err = MatGetColoring(Jp,MATCOLORING_ID,&iscolor);dCHK(err);
      err = MatFDColoringCreate(Jp,iscolor,&fdcolor);dCHK(err);
      err = ISColoringDestroy(iscolor);dCHK(err);
      err = MatFDColoringSetFunction(fdcolor,(PetscErrorCode(*)(void))StokesFunction,stk);dCHK(err);
      err = MatFDColoringSetFromOptions(fdcolor);dCHK(err);
      err = SNESSetJacobian(snes,J,Jp,SNESDefaultComputeJacobianColor,fdcolor);dCHK(err);
    } break;
    case 3:
      //err = SNESSetJacobian(snes,J,Jp,StokesJacobian,stk);dCHK(err);
    default: dERROR(1,"Not supported");
  }
  {
    KSP ksp;
    PC pc;
    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPGetPC(ksp,&pc);dCHK(err);
    err = PCSetType(pc,PCSHELL);dCHK(err);
    err = PCShellSetContext(pc,stk);dCHK(err);
    err = PCShellSetApply(pc,StokesPCApply);dCHK(err);
    err = PCShellSetSetUp(pc,StokesPCSetUp);dCHK(err);
  }
  err = SNESSetFromOptions(snes);dCHK(err);
  err = StokesGetSolutionVector(stk,&soln);dCHK(err);
  {
    dReal nrm;
    err = SNESComputeFunction(snes,soln,r);dCHK(err);
    err = VecNorm(r,NORM_2,&nrm);dCHK(err);
    err = PetscPrintf(comm,"Norm of discrete residual for exact solution %g\n",nrm);dCHK(err);
  }

  if (compute_explicit) {
    Mat expmat,expmat_full,Jfull;
    dInt m,n;
    dTruth flg;
    err = MatGetLocalSize(J,&m,&n);dCHK(err);
    err = MatCreateShell(comm,m,n,PETSC_DETERMINE,PETSC_DETERMINE,stk,&Jfull);dCHK(err);
    err = MatShellSetOperation(Jfull,MATOP_MULT,(void(*)(void))StokesShellMatMult_J);dCHK(err);

    err = MatComputeExplicitOperator(J,&expmat);dCHK(err);
    err = MatComputeExplicitOperator(Jfull,&expmat_full);dCHK(err);
    err = MatDestroy(Jfull);dCHK(err);
    //err = MatAXPY(expmat,-1,expmat_full,SAME_NONZERO_PATTERN);dCHK(err);
    err = MatSetOptionsPrefix(expmat,"explicit_");dCHK(err);
    err = MatSetFromOptions(expmat);dCHK(err);
    flg = dFALSE;
    err = PetscOptionsGetTruth(NULL,"-explicit_mat_view",&flg,NULL);dCHK(err);
    if (flg) {err = MatView(expmat,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);}
    flg = dFALSE;
    err = PetscOptionsGetTruth(NULL,"-explicit_mat_view_draw",&flg,NULL);dCHK(err);
    if (flg) {err = MatView(expmat,PETSC_VIEWER_DRAW_WORLD);dCHK(err);}
    err = MatDestroy(expmat);dCHK(err);
    err = MatDestroy(expmat_full);dCHK(err);
  }

  {                             /* Set null space */
    KSP ksp;
    MatNullSpace matnull;
    Mat mffd;
    dTruth isnull;
    Vec U,F;
    err = StokesGetNullSpace(stk,&matnull);dCHK(err);
    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPSetNullSpace(ksp,matnull);dCHK(err);
    err = MatNullSpaceRemove(matnull,soln,NULL);dCHK(err);
    /* The following is a real test of whether the null space is correct */
    err = MatCreateSNESMF(snes,&mffd);dCHK(err);
    err = MatSetFromOptions(mffd);dCHK(err);
    err = VecDuplicate(r,&U);dCHK(err);
    err = VecDuplicate(r,&F);dCHK(err);
    err = VecSet(U,0);dCHK(err);
    err = SNESComputeFunction(snes,U,F);dCHK(err);
    err = MatMFFDSetBase(mffd,U,F);dCHK(err);
    err = MatNullSpaceTest(matnull,mffd,&isnull);dCHK(err);
    err = MatDestroy(mffd);dCHK(err);
    if (!isnull) dERROR(1,"Vector is not in the null space of the MFFD operator");dCHK(err);
    err = MatNullSpaceTest(matnull,J,&isnull);dCHK(err);
    err = MatMult(J,r,F);dCHK(err);
    if (!isnull) dERROR(1,"Vector is not in the null space of J");dCHK(err);
    err = MatNullSpaceTest(matnull,Jp,&isnull);dCHK(err);
    if (!isnull) dERROR(1,"Vector is not in the null space of Jp");dCHK(err);
    err = VecDestroy(U);dCHK(err);
    err = VecDestroy(F);dCHK(err);
    err = MatNullSpaceDestroy(matnull);dCHK(err);
  }
  err = VecDuplicate(r,&x);dCHK(err);
  err = VecZeroEntries(r);dCHK(err);
  err = VecZeroEntries(x);dCHK(err);
  err = SNESSolve(snes,NULL,x);dCHK(err); /* ###  SOLVE  ### */
  if (0) {
    MatNullSpace matnull;
    KSP ksp;
    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPGetNullSpace(ksp,&matnull);dCHK(err); /* does not reference */
    err = MatNullSpaceRemove(matnull,x,NULL);dCHK(err);
  }
  if (!nocheck) {
    dReal anorm[2],anorminf,inorm[3];//,enorm[3],gnorm[3];
    //err = StokesErrorNorms(stk,x,enorm,gnorm);dCHK(err);
    err = VecNorm(r,NORM_1_AND_2,anorm);dCHK(err);
    err = VecNorm(r,NORM_INFINITY,&anorminf);dCHK(err);
    err = VecWAXPY(r,-1,soln,x);dCHK(err);
    err = VecNorm(r,NORM_1_AND_2,inorm);dCHK(err);
    err = VecNorm(r,NORM_INFINITY,&inorm[2]);dCHK(err);
    err = dPrintf(comm,"Algebraic residual        |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",anorm[0],anorm[1],anorminf);dCHK(err);
    err = dPrintf(comm,"Interpolation residual    |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",inorm[0],inorm[1],inorm[2]);dCHK(err);
    //err = dPrintf(comm,"Pointwise solution error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",enorm[0],enorm[1],enorm[2]);dCHK(err);
    //err = dPrintf(comm,"Pointwise gradient error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",gnorm[0],gnorm[1],gnorm[2]);dCHK(err);
  }

  err = VecDestroy(r);dCHK(err);
  err = VecDestroy(x);dCHK(err);
  err = VecDestroy(soln);dCHK(err);
  err = SNESDestroy(snes);dCHK(err);
  if (fdcolor) {err = MatFDColoringDestroy(fdcolor);dCHK(err);}
  if (J != Jp) {err = MatDestroy(J);dCHK(err);}
  err = MatDestroy(Jp);dCHK(err);
  err = StokesDestroy(stk);dCHK(err);
  err = PetscFinalize();dCHK(err);
  return 0;
}
