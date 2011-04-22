static const char help[] = "Solve nonlinear elasticity using dual order hp elements.\n"
  "The model problem is\n"
  "  -div(F S) = f\n"
  "where\n"
  "  F = 1 + Du                  (deformation gradient)\n"
  "  S = lambda tr(E) I + 2 mu E (Second Piola-Kirchoff tensor, Saint-Venant Kirchoff constitutive model)\n"
  "  E = (Du + Du^T + Du^T Du)/2 (Green-Lagrange tensor)\n"
  "  D is the symmetric gradient operator\n"
  "  mu and lambda are the Lame parameters\n\n";

#include <petscsnes.h>
#include <dohpfs.h>
#include <dohpviewer.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>

#include <elastexact.h>         /* Generated by elastexact.py */

static PetscLogEvent LOG_ElastShellMult;

struct ElastExact {
  void (*solution)(const struct ElastExactCtx*,const dReal x[3],dScalar u[3],dScalar du[9]);
  void (*forcing)(const struct ElastExactCtx*,const struct ElastParam*,const dReal x[3],dScalar f[3]);
};

struct ElastStore {
  dReal F[9];                /* Deformation gradient */
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
  dBool                 errorview,bdy100;
  dReal                 E0,nu;  /* Material parameters, here only for diagnostics */
  dQuadratureMethod     function_qmethod,jacobian_qmethod;
  dRulesetIterator      regioniter[EVAL_UB];
};

static dErr ElastSetMaterial(Elast elt,dReal E0,dReal nu)
{
  dFunctionBegin;
  elt->E0 = E0;
  elt->nu = nu;
  elt->param.lambda = E0*nu / ((1 + nu)*(1-2*nu));
  elt->param.mu     = E0 / (2*(1+nu));
  dFunctionReturn(0);
}

static dErr ElastCreate(MPI_Comm comm,Elast *elast)
{
  Elast elt;
  dErr err;

  dFunctionBegin;
  *elast = 0;
  err = dNew(struct ElastCtx,&elt);dCHK(err);
  elt->comm = comm;

  elt->constBDeg = 3;
  err = ElastSetMaterial(elt,1,0.3);dCHK(err);
  elt->bdy100 = dFALSE;
  elt->function_qmethod = dQUADRATURE_METHOD_FAST;
  elt->jacobian_qmethod = dQUADRATURE_METHOD_SPARSE;

  *elast = elt;
  dFunctionReturn(0);
}

static dErr ElastSetFromOptions(Elast elt)
{
  struct ElastExactCtx *exc = &elt->exactctx;
  dMesh mesh;
  dFS fs;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag dtag;
  dInt exact;
  dErr err;

  dFunctionBegin;
  exact = 0; exc->a = exc->b = exc->c = 1, exc->scale = 1;
  err = PetscOptionsBegin(elt->comm,NULL,"Elasticity options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-const_bdeg","Use constant isotropic degree on all elements","",elt->constBDeg,&elt->constBDeg,NULL);dCHK(err);
    err = PetscOptionsBool("-error_view","View errors","",elt->errorview,&elt->errorview,NULL);dCHK(err);
    err = PetscOptionsReal("-elast_E0","Young's modulus","",elt->E0,&elt->E0,NULL);dCHK(err);
    err = PetscOptionsReal("-elast_nu","Poisson ratio","",elt->nu,&elt->nu,NULL);dCHK(err);
    err = ElastSetMaterial(elt,elt->E0,elt->nu);dCHK(err);
    err = PetscOptionsBool("-bdy100","Only use boundary 100","",elt->bdy100,&elt->bdy100,NULL);dCHK(err);
    err = PetscOptionsEnum("-elast_f_qmethod","Quadrature method for residual evaluation/matrix-free","",dQuadratureMethods,(PetscEnum)elt->function_qmethod,(PetscEnum*)&elt->function_qmethod,NULL);dCHK(err);
    err = PetscOptionsEnum("-elast_jac_qmethod","Quadrature to use for Jacobian assembly","",dQuadratureMethods,(PetscEnum)elt->jacobian_qmethod,(PetscEnum*)&elt->jacobian_qmethod,NULL);dCHK(err);
    err = PetscOptionsInt("-exact","Exact solution choice","",exact,&exact,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_a","First scale parameter","",exc->a,&exc->a,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_b","Second scale parameter","",exc->b,&exc->b,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_c","Third scale parameter","",exc->c,&exc->c,NULL);dCHK(err);
    err = PetscOptionsReal("-exact_scale","Scale the amount of deformation [0,1]","",exc->scale,&exc->scale,NULL);dCHK(err);
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
  if (!elt->bdy100) {
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
    err = dRulesetDestroy(&ruleset);dCHK(err); /* Give ownership to iterator */
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
  err = dFSDestroy(&elt->fs);dCHK(err);
  err = dJacobiDestroy(&elt->jac);dCHK(err);
  err = dMeshDestroy(&elt->mesh);dCHK(err);
  err = VecDestroy(&elt->x);dCHK(err);
  err = VecDestroy(&elt->y);dCHK(err);
  for (dInt i=0; i<EVAL_UB; i++) {
    dRulesetIteratorDestroy(&elt->regioniter[i]);dCHK(err);
  }
  err = dFree(elt);dCHK(err);
  dFunctionReturn(0);
}

static inline dScalar DotColumn(const dScalar a[9],const dScalar b[9],dInt i,dInt j)
{return a[0*3+i]*b[0*3+j] + a[1*3+i]*b[1*3+j] + a[2*3+i]*b[2*3+j];}
static inline void DeformationGradient(dScalar F[9],const dScalar H[9])
{for (dInt i=0; i<3; i++) for (dInt j=0; j<3; j++) F[i*3+j] = (i==j) + H[i*3+j];}
static inline void ElastPointwiseComputeStore(struct ElastParam dUNUSED *prm,const dReal dUNUSED x[3],const dScalar dUNUSED u[3],const dScalar Du[9],struct ElastStore *st)
{DeformationGradient(st->F,Du);}

static inline void ElastPiolaKirchoff2(struct ElastParam *prm,const dScalar F[9],dScalar S[6])
{
  dScalar E[6],trace;
  // Green-Lagrangian strain tensor: E = (H + H' + H'*H)/2 = (F'*F - 1)/2
  E[0] = 0.5*(DotColumn(F,F,0,0) - 1);
  E[1] = 0.5*(DotColumn(F,F,1,1) - 1);
  E[2] = 0.5*(DotColumn(F,F,2,2) - 1);
  E[3] = 0.5*(DotColumn(F,F,0,1));
  E[4] = 0.5*(DotColumn(F,F,0,2));
  E[5] = 0.5*(DotColumn(F,F,1,2));
  trace = E[0] + E[1] + E[2];
  for (dInt i=0; i<6; i++)      /* Second Piola-Kirchoff stress tensor: S = lambda*tr(E)*1 + 2*mu*E */
    S[i] = prm->lambda*trace*(i<3) + 2*prm->mu*E[i];
}
static inline void ElastPiolaKirchoff2Jacobian(struct ElastParam *prm,const dScalar F[9],const dScalar dF[9],dScalar dS[6])
{
  dScalar dE[6],dtrace;
  // Green-Lagrangian strain tensor: E = (H + H' + H'*H)/2 = (F'*F - 1)/2
  // The perturbation of E in direction dF is: dE = (dH + dH' + dH'*H + H'*dH)/2
  dE[0] = 0.5*(DotColumn(F,dF,0,0) + DotColumn(dF,F,0,0));
  dE[1] = 0.5*(DotColumn(F,dF,1,1) + DotColumn(dF,F,1,1));
  dE[2] = 0.5*(DotColumn(F,dF,2,2) + DotColumn(dF,F,2,2));
  dE[3] = 0.5*(DotColumn(F,dF,0,1) + DotColumn(dF,F,0,1));
  dE[4] = 0.5*(DotColumn(F,dF,0,2) + DotColumn(dF,F,0,2));
  dE[5] = 0.5*(DotColumn(F,dF,1,2) + DotColumn(dF,F,1,2));
  dtrace = dE[0] + dE[1] + dE[2];
  for (dInt i=0; i<6; i++)      // Saint-Venant Kirchoff model: S = lambda*tr(E)*1 + 2*mu*E
    dS[i] = prm->lambda*dtrace*(i<3) + 2*prm->mu*dE[i];
}
static inline void ElastPointwiseFunction(struct ElastParam *prm,struct ElastExact *exact,struct ElastExactCtx *exactctx,
                                          const dReal x[3],dReal weight,const dScalar u[3],const dScalar Du[9],
                                          struct ElastStore *st,dScalar v[3],dScalar Dv[9])
{
  dScalar f[3],S[6];
  exact->forcing(exactctx,prm,x,f);
  for (dInt i=0; i<3; i++) v[i] = -weight * f[i]; // Body force
  ElastPointwiseComputeStore(prm,x,u,Du,st);
  ElastPiolaKirchoff2(prm,st->F,S);
  dTensorMultGESY3(Dv,st->F,S); // Weak form integrand is:   Dv:Pi, Pi = F*S, F = 1+H, H = Dv
  for (dInt i=0; i<9; i++) Dv[i] *= weight;
}

static inline void ElastPointwiseJacobian(struct ElastParam *prm,const struct ElastStore *restrict st,dReal weight,
                                          const dScalar dUNUSED u[restrict static 3],const dScalar Du[restrict static 9],
                                          dScalar v[restrict static 3],dScalar Dv[restrict static 9])
{
  dScalar S[6],dS[6];
  v[0] = v[1] = v[2] = 0;
  ElastPiolaKirchoff2(prm,st->F,S);             // Recompute second Piola-Kirchoff tensor
  ElastPiolaKirchoff2Jacobian(prm,st->F,Du,dS); // Compute derivative of S in direction Du
  dTensorMultGESY3(Dv,st->F,dS);    // Dv : F dS
  dTensorMultAddGESY3(Dv,Du,S);     // Dv : dF S
  for (dInt i=0; i<9; i++) Dv[i] *= weight;
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
      //err = dRealTableView(P*3,P*3,&K[0][0][0][0],PETSC_VIEWER_STDOUT_WORLD,"K");dCHK(err);
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
      elt->exact.solution(&elt->exactctx,x[i],uu,&duu[0][0]);
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
    elt->exact.solution(&elt->exactctx,&coords[3*i],&x[i*bs],du_unused);
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
  Mat J,Jp;
  Vec r,x,soln;
  SNES snes;
  dBool  nojshell = dFALSE,nocheck = dFALSE,viewdhm = dFALSE;
  dErr err;

  err = dInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
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
    err = PetscOptionsBool("-viewdhm","View the solution","",viewdhm,&viewdhm,NULL);dCHK(err);
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
  if (viewdhm) {
    dViewer view;
    err = PetscViewerCreate(comm,&view);dCHK(err);
    err = PetscViewerSetType(view,PETSCVIEWERDHM);dCHK(err);
    err = PetscViewerFileSetName(view,"elast.dhm");dCHK(err);
    err = PetscViewerFileSetMode(view,FILE_MODE_WRITE);dCHK(err);
    err = dViewerDHMSetTimeUnits(view,"hour",PETSC_PI*1e7/3600);dCHK(err);
    err = dViewerDHMSetTime(view,0.1);dCHK(err);
    err = VecView(x,view);dCHK(err);
    err = PetscViewerDestroy(&view);dCHK(err);
  }

  err = VecDestroy(&r);dCHK(err);
  err = VecDestroy(&x);dCHK(err);
  err = VecDestroy(&soln);dCHK(err);
  err = SNESDestroy(&snes);dCHK(err);
  if (J != Jp) {err = MatDestroy(&J);dCHK(err);}
  err = MatDestroy(&Jp);dCHK(err);
  err = ElastDestroy(elt);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
