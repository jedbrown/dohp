static const char help[] = "Solve non-Newtonian Stokes problem using dual order hp elements.\n"
  "The model problem is\n"
  "  -div(eta Du) + grad(p) = f\n"
  "                  div(u) = g\n"
  "where\n"
  "  D is the symmetric gradient operator\n"
  "  eta(gamma) = B (0.5*eps^2 + gamma)^{(p-2)/2}\n"
  "  gamma = Du : Du/2\n"
  "The weak form is\n"
  "  int_Omega eta Dv:Du - p div(v) - q div(u) - f_u.v - f_p.q = 0\n"
  "with Jacobian\n"
  "  int_Omega eta Dv:Du + eta' (Dv:Dw)(Dw:Du) - p div(v) - q div(u) = 0\n"
  "The problem is linear for p=2, an incompressible for g=0\n\n";

#include <petscsnes.h>
#include <dohpstring.h>
#include <dohpviewer.h>

#include "stokesimpl.h"

PetscFunctionList StokesCaseList = NULL;

#define StokesCaseType char*

dErr StokesCaseRegister(const char *name,StokesCaseCreateFunction screate)
{
  dErr err;
  dFunctionBegin;
  err = PetscFunctionListAdd(PETSC_COMM_WORLD,&StokesCaseList,name,"",(void(*)(void))screate);dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseFind(const char *name,StokesCaseCreateFunction *screate)
{
  dErr err;

  dFunctionBegin;
  err = PetscFunctionListFind(PETSC_COMM_WORLD,StokesCaseList,name,PETSC_FALSE,(void(**)(void))screate);dCHK(err);
  if (!*screate) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Stokes Case \"%s\" could not be found",name);
  dFunctionReturn(0);
}
static dErr StokesCaseSetType(StokesCase scase,const StokesCaseType type)
{
  dErr err;
  StokesCaseCreateFunction f;

  dFunctionBegin;
  err = StokesCaseFind(type,&f);dCHK(err);
  err = (*f)(scase);dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseSetFromOptions(StokesCase scase)
{
  struct StokesRheology *rheo = &scase->rheo;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsBegin(scase->comm,NULL,"StokesCase_Exact options",__FILE__);dCHK(err); {
    err = PetscOptionsReal("-rheo_B","Rate factor (rheology)","",rheo->B,&rheo->B,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_eps","Regularization (rheology)","",rheo->eps,&rheo->eps,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_p","Power p=1+1/n where n is Glen exponent","",rheo->p,&rheo->p,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_nu","Poisson ratio to use in Jacobian","",rheo->nu,&rheo->nu,NULL);dCHK(err);
    err = PetscOptionsReal("-gravity","Nondimensional gravitational force","",scase->gravity,&scase->gravity,NULL);dCHK(err);
    if (scase->setfromoptions) {err = (*scase->setfromoptions)(scase);dCHK(err);}
  } err = PetscOptionsEnd();dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseDestroy(StokesCase *scase)
{
  dErr err;

  dFunctionBegin;
  if ((*scase)->destroy) {err = ((*scase)->destroy)(*scase);dCHK(err);}
  err = dFree(*scase);dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseRegisterAll(void)
{
  dErr err;

  dFunctionBegin;
  err = StokesCaseRegisterAll_Exact();dCHK(err);
#if defined dHAVE_GDAL
  err = StokesCaseRegisterAll_Jako();dCHK(err);
#endif
  dFunctionReturn(0);
}

#define ALEN(a) ((dInt)(sizeof(a)/sizeof(a)[0]))
static PetscLogEvent LOG_StokesShellMult;

static dErr StokesGetNullSpace(Stokes stk,MatNullSpace *matnull);
static dErr StokesShellMatMult_All_IorA(Mat A,Vec gx,Vec gy,Vec gz,InsertMode,StokesMultMode);
static dErr MatMult_Nest_StokesCoupled(Mat J,Vec gx,Vec gy);
static dErr StokesShellMatMult_A(Mat A,Vec gx,Vec gy) {return StokesShellMatMult_All_IorA(A,gx,gy,NULL,INSERT_VALUES,STOKES_MULT_A);}
static dErr StokesShellMatMult_Bt(Mat A,Vec gx,Vec gy) {return StokesShellMatMult_All_IorA(A,gx,gy,NULL,INSERT_VALUES,STOKES_MULT_Bt);}
static dErr StokesShellMatMult_B(Mat A,Vec gx,Vec gy) {return StokesShellMatMult_All_IorA(A,gx,gy,NULL,INSERT_VALUES,STOKES_MULT_B);}
static dErr StokesShellMatMultAdd_A(Mat A,Vec gx,Vec gy,Vec gz) {return StokesShellMatMult_All_IorA(A,gx,gy,gz,ADD_VALUES,STOKES_MULT_A);}
static dErr StokesShellMatMultAdd_Bt(Mat A,Vec gx,Vec gy,Vec gz) {return StokesShellMatMult_All_IorA(A,gx,gy,gz,ADD_VALUES,STOKES_MULT_Bt);}
static dErr StokesShellMatMultAdd_B(Mat A,Vec gx,Vec gy,Vec gz) {return StokesShellMatMult_All_IorA(A,gx,gy,gz,ADD_VALUES,STOKES_MULT_B);}
static dErr MatGetVecs_Stokes(Mat,Vec*,Vec*);

static dErr StokesCreate(MPI_Comm comm,Stokes *stokes)
{
  Stokes stk;
  dErr err;

  dFunctionBegin;
  *stokes = 0;
  err = dNew(struct _n_Stokes,&stk);dCHK(err);
  stk->comm = comm;

  stk->constBDeg     = 3;
  stk->pressureCodim = 2;
  stk->dirichlet[0]  = 100;
  stk->dirichlet[1]  = 200;
  stk->dirichlet[2]  = 300;
  stk->alldirichlet  = dTRUE;
  stk->function_qmethod = dQUADRATURE_METHOD_FAST;
  stk->jacobian_qmethod = dQUADRATURE_METHOD_SPARSE;

  err = dCalloc(sizeof(*stk->scase),&stk->scase);dCHK(err);
  stk->scase->rheo.B   = 1;
  stk->scase->rheo.eps = 1;
  stk->scase->rheo.p   = 2;
  stk->scase->rheo.nu  = 0;

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
  if (nu==np) dERROR(PETSC_COMM_SELF,1,"Degenerate case, don't know which space to copy");
  if (x) {
    if (n == nu) {
      err = VecDuplicate(stk->gvelocity,x);dCHK(err);
    } else if (n == np) {
      err = VecDuplicate(stk->gpressure,x);dCHK(err);
    } else dERROR(PETSC_COMM_SELF,1,"sizes do not agree with either space");
  }
  if (y) {
    if (n == nu) {
      err = VecDuplicate(stk->gvelocity,y);dCHK(err);
    } else if (n == np) {
      err = VecDuplicate(stk->gpressure,y);dCHK(err);
    } else dERROR(PETSC_COMM_SELF,1,"sizes do not agree with either space");
  }
  dFunctionReturn(0);
}

static dErr StokesSetFromOptions(Stokes stk)
{
  char scasename[256] = "Exact1";
  dMesh mesh;
  dFS fsu,fsp;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag dtag,dptag;
  dErr err;

  dFunctionBegin;
  err = dStrcpyS(stk->mattype_Ap,sizeof(stk->mattype_Ap),MATBAIJ);dCHK(err);
  err = dStrcpyS(stk->mattype_Dp,sizeof(stk->mattype_Dp),MATAIJ);dCHK(err);
  err = PetscOptionsBegin(stk->comm,NULL,"Stokesicity options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-const_bdeg","Use constant isotropic degree on all elements","",stk->constBDeg,&stk->constBDeg,NULL);dCHK(err);
    err = PetscOptionsInt("-pressure_codim","Reduce pressure space by this factor","",stk->pressureCodim,&stk->pressureCodim,NULL);dCHK(err);
    err = PetscOptionsBool("-cardinal_mass","Assemble diagonal mass matrix","",stk->cardinalMass,&stk->cardinalMass,NULL);dCHK(err);
    err = PetscOptionsList("-stokes_Ap_mat_type","Matrix type for velocity operator","",MatList,stk->mattype_Ap,stk->mattype_Ap,sizeof(stk->mattype_Ap),NULL);dCHK(err);
    err = PetscOptionsList("-stokes_Dp_mat_type","Matrix type for pressure operator","",MatList,stk->mattype_Dp,stk->mattype_Dp,sizeof(stk->mattype_Dp),NULL);dCHK(err);
    err = PetscOptionsEnum("-stokes_f_qmethod","Quadrature method for residual evaluation/matrix-free","",dQuadratureMethods,(PetscEnum)stk->function_qmethod,(PetscEnum*)&stk->function_qmethod,NULL);dCHK(err);
    err = PetscOptionsEnum("-stokes_jac_qmethod","Quadrature to use for Jacobian assembly","",dQuadratureMethods,(PetscEnum)stk->jacobian_qmethod,(PetscEnum*)&stk->jacobian_qmethod,NULL);dCHK(err);
    {
      dBool flg; dInt n = ALEN(stk->dirichlet);
      err = PetscOptionsIntArray("-dirichlet","List of boundary sets on which to impose Dirichlet conditions","",stk->dirichlet,&n,&flg);dCHK(err);
      if (flg) {
        for (dInt i=n; i<ALEN(stk->dirichlet); i++) stk->dirichlet[i] = 0; /* Clear out any leftover values */
        if (n < 3) stk->alldirichlet = dFALSE;                             /* @bug More work to determine independent of the mesh whether all the boundaries are Dirichlet */
      }
    }
    err = PetscOptionsList("-stokes_case","Which sort of case to run","",StokesCaseList,scasename,scasename,sizeof(scasename),NULL);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);

  err = dMeshCreate(stk->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);
  err = dMeshGetRoot(mesh,&domain);dCHK(err); /* Need a taggable set */
  err = dMeshSetDuplicateEntsOnly(mesh,domain,&domain);dCHK(err);
  err = PetscObjectSetName((PetscObject)mesh,"dMesh_0");dCHK(err);

  err = dJacobiCreate(stk->comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);

  err = dMeshCreateRuleTagIsotropic(mesh,domain,"stokes_efs_velocity_degree",stk->constBDeg,&dtag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,"stokes_efs_pressure_degree",stk->constBDeg-stk->pressureCodim,&dptag);dCHK(err);

  err = dFSCreate(stk->comm,&fsu);dCHK(err);
  err = dFSSetBlockSize(fsu,3);dCHK(err);
  err = dFSSetMesh(fsu,mesh,domain);dCHK(err);
  err = dFSSetDegree(fsu,jac,dtag);dCHK(err);
  for (dInt i=0; i<ALEN(stk->dirichlet) && stk->dirichlet[i]>0; i++) {
    err = dFSRegisterBoundary(fsu,stk->dirichlet[i],dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  }
  err = PetscObjectSetOptionsPrefix((dObject)fsu,"u");dCHK(err);
  err = dFSSetFromOptions(fsu);dCHK(err);
  err = PetscObjectSetName((PetscObject)fsu,"dFS_U_0");dCHK(err);
  stk->fsu = fsu;

  err = dFSCreate(stk->comm,&fsp);dCHK(err);
  err = dFSSetMesh(fsp,mesh,domain);dCHK(err);
  err = dFSSetDegree(fsp,jac,dptag);dCHK(err);
  err = PetscObjectSetOptionsPrefix((dObject)fsp,"p");dCHK(err);
  /* No boundaries, the pressure space has Neumann conditions when Dirichlet velocity conditions are applied */
  err = dFSSetFromOptions(fsp);dCHK(err);
  err = PetscObjectSetName((PetscObject)fsp,"dFS_P_0");dCHK(err);
  stk->fsp = fsp;

  err = dFSCreateExpandedVector(fsu,&stk->xu);dCHK(err);
  err = VecDuplicate(stk->xu,&stk->yu);dCHK(err);

  err = dFSCreateExpandedVector(fsp,&stk->xp);dCHK(err);
  err = VecDuplicate(stk->xp,&stk->yp);dCHK(err);

  {
    dInt nu,np,rstart,nul,npl;
    IS   ublock,pblock;
    Vec  Vc,Vgh,Pc,Pgh;
    err = dFSCreateGlobalVector(stk->fsu,&stk->gvelocity);dCHK(err);
    err = dFSCreateGlobalVector(stk->fsp,&stk->gpressure);dCHK(err);
    err = PetscObjectSetName((PetscObject)stk->gvelocity,"Velocity");dCHK(err);
    err = PetscObjectSetName((PetscObject)stk->gpressure,"Pressure");dCHK(err);
    err = VecGetLocalSize(stk->gvelocity,&nu);dCHK(err);
    err = VecGetLocalSize(stk->gpressure,&np);dCHK(err);
    err = VecCreateMPI(stk->comm,nu+np,PETSC_DETERMINE,&stk->gpacked);dCHK(err);
    err = VecGetOwnershipRange(stk->gpacked,&rstart,NULL);dCHK(err);
    err = ISCreateStride(stk->comm,nu,rstart,1,&ublock);dCHK(err);
    err = ISCreateStride(stk->comm,np,rstart+nu,1,&pblock);dCHK(err);
    err = ISSetBlockSize(ublock,3);dCHK(err);
    err = VecScatterCreate(stk->gpacked,ublock,stk->gvelocity,NULL,&stk->extractVelocity);dCHK(err);
    err = VecScatterCreate(stk->gpacked,pblock,stk->gpressure,NULL,&stk->extractPressure);dCHK(err);
    stk->ublock = ublock;
    stk->pblock = pblock;
    /* Create local index sets */
    err = VecDohpGetClosure(stk->gvelocity,&Vc);dCHK(err);
    err = VecDohpGetClosure(stk->gpressure,&Pc);dCHK(err);
    err = VecGhostGetLocalForm(Vc,&Vgh);dCHK(err);
    err = VecGhostGetLocalForm(Pc,&Pgh);dCHK(err);
    err = VecGetLocalSize(Vgh,&nul);dCHK(err);
    err = VecGetLocalSize(Pgh,&npl);dCHK(err);
    err = VecGhostRestoreLocalForm(Vc,&Vgh);dCHK(err);
    err = VecGhostRestoreLocalForm(Pc,&Pgh);dCHK(err);
    err = VecDohpRestoreClosure(stk->gvelocity,&Vc);dCHK(err);
    err = VecDohpRestoreClosure(stk->gpressure,&Pc);dCHK(err);
    err = ISCreateStride(PETSC_COMM_SELF,nul,0,1,&stk->lublock);dCHK(err);
    err = ISCreateStride(PETSC_COMM_SELF,npl,nul,1,&stk->lpblock);dCHK(err);
    err = ISSetBlockSize(stk->lublock,3);dCHK(err);
  }
  err = dJacobiDestroy(&jac);dCHK(err);
  err = dMeshDestroy(&mesh);dCHK(err);

  err = StokesCaseSetType(stk->scase,scasename);dCHK(err);
  err = dFSGetBoundingBox(stk->fsu,stk->scase->bbox);dCHK(err);
  err = StokesCaseSetFromOptions(stk->scase);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesGetRegionIterator(Stokes stk,StokesEvaluation eval,dRulesetIterator *riter)
{
  dErr err;

  dFunctionBegin;
  if (!stk->regioniter[eval]) {
    dRulesetIterator iter;
    dRuleset ruleset;
    dFS cfs;
    dMeshESH domain;
    dQuadratureMethod qmethod;
    switch (eval) {
    case EVAL_FUNCTION: qmethod = stk->function_qmethod; break;
    case EVAL_JACOBIAN: qmethod = stk->jacobian_qmethod; break;
    default: dERROR(stk->comm,PETSC_ERR_ARG_OUTOFRANGE,"Unknown evaluation context");
    }
    err = dFSGetDomain(stk->fsu,&domain);dCHK(err);
    err = dFSGetPreferredQuadratureRuleSet(stk->fsu,domain,dTYPE_REGION,dTOPO_ALL,qmethod,&ruleset);dCHK(err);
    err = dFSGetCoordinateFS(stk->fsu,&cfs);dCHK(err);
    err = dRulesetCreateIterator(ruleset,cfs,&iter);dCHK(err);
    err = dRulesetDestroy(&ruleset);dCHK(err); /* Give ownership to iterator */
    err = dRulesetIteratorAddFS(iter,stk->fsu);dCHK(err);
    err = dRulesetIteratorAddFS(iter,stk->fsp);dCHK(err);
    if (eval == EVAL_FUNCTION) {err = dRulesetIteratorAddStash(iter,0,sizeof(struct StokesStore));dCHK(err);}
    stk->regioniter[eval] = iter;
  }
  *riter = stk->regioniter[eval];
  dFunctionReturn(0);
}

static dErr StokesExtractGlobalSplit(Stokes stk,Vec gx,Vec *gxu,Vec *gxp)
{
  dErr err;

  dFunctionBegin;
  if (gxu) {
    *gxu = stk->gvelocity;
    err = VecScatterBegin(stk->extractVelocity,gx,*gxu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (stk->extractVelocity,gx,*gxu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  if (gxp) {
    *gxp = stk->gpressure;
    err = VecScatterBegin(stk->extractPressure,gx,*gxp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (stk->extractPressure,gx,*gxp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr StokesCommitGlobalSplit(Stokes stk,Vec *gxu,Vec *gxp,Vec gy,InsertMode imode)
{
  dErr err;

  dFunctionBegin;
  dASSERT(*gxu == stk->gvelocity);
  dASSERT(*gxp == stk->gpressure);
  err = VecScatterBegin(stk->extractVelocity,*gxu,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (stk->extractVelocity,*gxu,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(stk->extractPressure,*gxp,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (stk->extractPressure,*gxp,gy,imode,SCATTER_REVERSE);dCHK(err);
  *gxu = NULL;
  *gxp = NULL;
  dFunctionReturn(0);
}

static dErr StokesDestroy(Stokes *instk)
{
  Stokes stk = *instk;
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(&stk->fsu);dCHK(err);
  err = dFSDestroy(&stk->fsp);dCHK(err);
  err = VecDestroy(&stk->xu);dCHK(err);
  err = VecDestroy(&stk->yu);dCHK(err);
  err = VecDestroy(&stk->xp);dCHK(err);
  err = VecDestroy(&stk->yp);dCHK(err);
  err = VecDestroy(&stk->gvelocity);dCHK(err);
  err = VecDestroy(&stk->gpressure);dCHK(err);
  err = VecDestroy(&stk->gpacked);dCHK(err);
  err = VecScatterDestroy(&stk->extractVelocity);dCHK(err);
  err = VecScatterDestroy(&stk->extractPressure);dCHK(err);
  err = ISDestroy(&stk->ublock);dCHK(err);
  err = ISDestroy(&stk->pblock);dCHK(err);
  err = ISDestroy(&stk->lublock);dCHK(err);
  err = ISDestroy(&stk->lpblock);dCHK(err);
  for (dInt i=0; i<EVAL_UB; i++) {err = dRulesetIteratorDestroy(&stk->regioniter[i]);dCHK(err);}
  err = StokesCaseDestroy(&stk->scase);dCHK(err);
  err = dFree(*instk);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesGetMatrices(Stokes stk,dBool use_jblock,Mat *J,Mat *Jp)
{
  dErr err;
  dInt m,nu,np;
  Mat A,B,Bt,Ap,Dp;
  IS splitis[2];

  dFunctionBegin;
  err = VecGetLocalSize(stk->gpacked,&m);dCHK(err);
  err = VecGetLocalSize(stk->gvelocity,&nu);dCHK(err);
  err = VecGetLocalSize(stk->gpressure,&np);dCHK(err);

  /* Create high-order matrix for diagonal velocity block, with context \a stk */
  err = MatCreateShell(stk->comm,nu,nu,PETSC_DETERMINE,PETSC_DETERMINE,stk,&A);dCHK(err);
  err = MatShellSetOperation(A,MATOP_GET_VECS,(void(*)(void))MatGetVecs_Stokes);dCHK(err);
  err = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))StokesShellMatMult_A);dCHK(err);
  err = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void))StokesShellMatMult_A);dCHK(err);
  err = MatShellSetOperation(A,MATOP_MULT_ADD,(void(*)(void))StokesShellMatMultAdd_A);dCHK(err);
  err = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))StokesShellMatMultAdd_A);dCHK(err);
  err = MatSetOptionsPrefix(A,"A_");dCHK(err);

  /* Create off-diagonal high-order matrix, with context \a stk */
  err = MatCreateShell(stk->comm,np,nu,PETSC_DETERMINE,PETSC_DETERMINE,stk,&B);dCHK(err);
  err = MatShellSetOperation(B,MATOP_GET_VECS,(void(*)(void))MatGetVecs_Stokes);dCHK(err);
  err = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))StokesShellMatMult_B);dCHK(err);
  err = MatShellSetOperation(B,MATOP_MULT_TRANSPOSE,(void(*)(void))StokesShellMatMult_Bt);dCHK(err);
  err = MatShellSetOperation(B,MATOP_MULT_ADD,(void(*)(void))StokesShellMatMultAdd_B);dCHK(err);
  err = MatShellSetOperation(B,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))StokesShellMatMultAdd_Bt);dCHK(err);
  err = MatCreateTranspose(B,&Bt);dCHK(err);
  err = MatSetOptionsPrefix(B,"B_");dCHK(err);
  err = MatSetOptionsPrefix(Bt,"Bt_");dCHK(err);

  splitis[0] = stk->ublock;
  splitis[1] = stk->pblock;
  /* Create the matrix-free operator */
  err = MatCreateNest(stk->comm,2,splitis,2,splitis,((Mat[]){A,Bt,B,NULL}),J);dCHK(err);
  err = MatSetOptionsPrefix(*J,"J_");dCHK(err);
  err = MatSetFromOptions(*J);dCHK(err);
  if (!use_jblock) {
    err = MatShellSetOperation(*J,MATOP_MULT,(void(*)(void))MatMult_Nest_StokesCoupled);dCHK(err);
    err = MatShellSetOperation(*J,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_Nest_StokesCoupled);dCHK(err);
  }

  err = MatDestroy(&A);dCHK(err);
  err = MatDestroy(&Bt);dCHK(err);
  err = MatDestroy(&B);dCHK(err);

  /* Create real matrix to be used for preconditioning */
  err = dFSCreateMatrix(stk->fsu,stk->mattype_Ap,&Ap);dCHK(err);
  err = dFSCreateMatrix(stk->fsp,stk->mattype_Dp,&Dp);dCHK(err);
  err = MatSetOptionsPrefix(Ap,"Ap_");dCHK(err);
  err = MatSetOptionsPrefix(Dp,"Dp_");dCHK(err);
  err = MatSetOption(Ap,MAT_SYMMETRIC,PETSC_TRUE);dCHK(err);
  err = MatSetOption(Dp,MAT_SYMMETRIC,PETSC_TRUE);dCHK(err);
  err = MatSetFromOptions(Ap);dCHK(err);
  err = MatSetFromOptions(Dp);dCHK(err);

  err = MatCreateNest(stk->comm,2,splitis,2,splitis,((Mat[]){Ap,NULL,NULL,Dp}),Jp);dCHK(err);
  err = MatSetOptionsPrefix(*Jp,"Jp_");dCHK(err);
  err = MatSetFromOptions(*Jp);dCHK(err);

  {                             /* Allocate for the pressure Poisson, used by PCLSC */
    Mat L;
    Vec Mdiag;
    err = dFSCreateMatrix(stk->fsp,stk->mattype_Dp,&L);dCHK(err);
    err = MatSetOptionsPrefix(L,"stokes_L_");dCHK(err);
    err = MatSetFromOptions(L);dCHK(err);
    err = PetscObjectCompose((dObject)Dp,"LSC_L",(dObject)L);dCHK(err);
    err = PetscObjectCompose((dObject)Dp,"LSC_Lp",(dObject)L);dCHK(err);
    err = MatDestroy(&L);dCHK(err); /* don't keep a reference */
    err = VecDuplicate(stk->gvelocity,&Mdiag);dCHK(err);
    err = PetscObjectCompose((dObject)Dp,"LSC_M_diag",(dObject)Mdiag);
    err = VecDestroy(&Mdiag);dCHK(err); /* don't keep a reference */
  }

  err = MatDestroy(&Ap);dCHK(err); /* release reference to Jp */
  err = MatDestroy(&Dp);dCHK(err); /* release reference to Jp */
  dFunctionReturn(0);
}

static inline void StokesPointwiseComputeStore(const struct StokesRheology *rheo,const dReal dUNUSED x[3],const dScalar Du[],struct StokesStore *st)
{
  dScalar gamma_reg = 0.5*dSqr(rheo->eps) + 0.5*dColonSymScalar3(Du,Du);
  st->eta = rheo->B * pow(gamma_reg,0.5*(rheo->p-2));
  st->deta = 0.5*(rheo->p-2) * st->eta / gamma_reg;
  for (dInt i=0; i<6; i++) st->Du[i] = Du[i];
}

static inline dReal StokesLambda(const struct StokesRheology *rheo,const struct StokesStore *st)
{
  return ((1 + 0*st->eta) * rheo->nu) / (1 - 2*rheo->nu);
}

static inline void StokesPointwiseFunction(StokesCase scase,
                                           const dReal x[3],dReal weight,const dScalar Du[6],dScalar p,
                                           struct StokesStore *st,dScalar v[3],dScalar Dv[6],dScalar *q)
{
  const struct StokesRheology *rheo = &scase->rheo;
  dScalar fu[3],fp,trDu = Du[0]+Du[1]+Du[2];
  dReal lambda;
  StokesPointwiseComputeStore(rheo,x,Du,st);
  lambda = 0*StokesLambda(rheo,st);
  scase->forcing(scase,x,fu,&fp);
  for (dInt i=0; i<3; i++) v[i] = -weight * fu[i]; /* Coefficient of \a v in weak form, only appears in forcing term */
  *q   = -weight * (Du[0]+Du[1]+Du[2] + fp);       /* -q tr(Du) - forcing, note tr(Du) = div(u) */
  for (dInt i=0; i<3; i++) Dv[i] = weight * (st->eta * Du[i] + lambda*trDu - p); /* eta Dv:Du - p tr(Dv) */
  for (dInt i=3; i<6; i++) Dv[i] = weight * st->eta * Du[i];       /* eta Dv:Du */
}

static inline void StokesPointwiseJacobian(const struct StokesRheology *rheo,const struct StokesStore *restrict st,dReal weight,
                                           const dScalar Du[restrict static 6],dScalar p,
                                           dScalar Dv[restrict static 6],dScalar *restrict q)
{
  const dReal lambda = 0*StokesLambda(rheo,st);
                                /* Coefficients in weak form of Jacobian */
  const dScalar deta_colon = st->deta*dColonSymScalar3(st->Du,Du);                      // eta' Dw:Du
  for (dInt i=0; i<6; i++) Dv[i] = weight * (st->eta * Du[i] + deta_colon * st->Du[i]   // eta Dv:Du + eta' (Dv:Dw)(Dw:Du)
                                             + lambda * (i<3) * (Du[0] + Du[1] + Du[2]) // Penalty: lambda div(v) div(u)
                                             - p*(i<3));                                // - p tr(Dv)
  *q = -weight*(Du[0]+Du[1]+Du[2]);                                                     // -q tr(Du)
}

static inline void StokesPointwiseJacobian_A(const struct StokesRheology *rheo,const struct StokesStore *restrict st,dReal weight,const dScalar Du[restrict static 6],dScalar Dv[restrict static 6])
{
  const dScalar deta_colon = st->deta*dColonSymScalar3(st->Du,Du);
  const dReal lambda = StokesLambda(rheo,st);
  for (dInt i=0; i<6; i++) Dv[i] = weight * (st->eta*Du[i] + deta_colon*st->Du[i] + lambda*(i<3)*(Du[0]+Du[1]+Du[2]));
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
  Stokes           stk = ctx;
  dErr             err;
  Vec              Coords,gxu,gxp;
  dRulesetIterator iter;

  dFunctionBegin;
  err = StokesExtractGlobalSplit(stk,gx,&gxu,&gxp);dCHK(err);
  err = StokesGetRegionIterator(stk,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(stk->fsu,&Coords);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, gxu,dFS_INHOMOGENEOUS,gxu,dFS_INHOMOGENEOUS, gxp,dFS_INHOMOGENEOUS,gxp,dFS_INHOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*v)[3],(*dv)[9],*p,*q;
    dInt Q;
    struct StokesStore *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,&v,&dv, &p,NULL,&q,NULL);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar Du[6],Dv[6];
      dTensorSymCompress3(du[i],Du);
      StokesPointwiseFunction(stk->scase,x[i],jw[i],Du,p[i],&stash[i],v[i],Dv,&q[i]);
      dTensorSymUncompress3(Dv,dv[i]);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, v,dv, q,NULL);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  //err = VecView(gxp,0);dCHK(err);
  err = StokesCommitGlobalSplit(stk,&gxu,&gxp,gy,INSERT_VALUES);dCHK(err);
  dFunctionReturn(0);
}

static dErr MatMult_Nest_StokesCoupled(Mat J,Vec gx,Vec gy)
{
  Stokes           stk;
  Vec              Coords,gxu,gxp;
  dRulesetIterator iter;
  dErr             err;
  Mat              A;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_StokesShellMult,J,gx,gy,0);dCHK(err);
  err = MatNestGetSubMat(J,0,0,&A);dCHK(err);
  err = MatShellGetContext(A,(void**)&stk);dCHK(err);
  err = StokesExtractGlobalSplit(stk,gx,&gxu,&gxp);dCHK(err);
  err = StokesGetRegionIterator(stk,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(stk->fsu,&Coords);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, gxu,dFS_HOMOGENEOUS,gxu,dFS_HOMOGENEOUS, gxp,dFS_HOMOGENEOUS,gxp,dFS_HOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*dv)[9],*p,*q;
    dInt Q;
    struct StokesStore *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,&dv, &p,NULL,&q,NULL);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar Du[6],Dv[6];
      dTensorSymCompress3(du[i],Du);
      StokesPointwiseJacobian(&stk->scase->rheo,&stash[i],jw[i],Du,p[i],Dv,&q[i]);
      dTensorSymUncompress3(Dv,dv[i]);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,dv, q,NULL);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = StokesCommitGlobalSplit(stk,&gxu,&gxp,gy,INSERT_VALUES);dCHK(err);
  err = PetscLogEventEnd(LOG_StokesShellMult,J,gx,gy,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesShellMatMult_All_IorA(Mat A,Vec gx,Vec gy,Vec gz,InsertMode imode,StokesMultMode mmode)
{
  Stokes           stk;
  dRulesetIterator iter;
  Vec              Coords;
  dErr             err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_StokesShellMult,A,gx,gy,gz);dCHK(err);
  err = MatShellGetContext(A,(void**)&stk);dCHK(err);
  {  /* Check that we have correct sizes */
    dInt nu,np,nx,ny;
    err = VecGetSize(stk->gvelocity,&nu);dCHK(err);
    err = VecGetSize(stk->gpressure,&np);dCHK(err);
    err = VecGetSize(gx,&nx);dCHK(err);
    err = VecGetSize(gy,&ny);dCHK(err);
    switch (mmode) {
    case STOKES_MULT_A: dASSERT(nx==nu && ny==nu); break;
    case STOKES_MULT_Bt: dASSERT(nx==np && ny==nu); break;
    case STOKES_MULT_B: dASSERT(nx==nu && ny==np); break;
    default: dERROR(PETSC_COMM_SELF,1,"Sizes do not match, unknown mult operation");
    }
  }

  switch (imode) {
  case INSERT_VALUES:
    if (gz) dERROR(stk->comm,PETSC_ERR_ARG_INCOMP,"Cannot use INSERT_VALUES and set gz");
    gz = gy;
    err = VecZeroEntries(gz);dCHK(err);
    break;
  case ADD_VALUES:
    if (gz != gy) {
      err = VecCopy(gy,gz);dCHK(err);
    }
    break;
  default: dERROR(stk->comm,PETSC_ERR_ARG_OUTOFRANGE,"unsupported imode");
  }

  err = StokesGetRegionIterator(stk,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(stk->fsu,&Coords);dCHK(err);
  switch (mmode) {
  case STOKES_MULT_A:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, gx,dFS_HOMOGENEOUS,gz,dFS_HOMOGENEOUS, NULL,              NULL);dCHK(err);
    break;
  case STOKES_MULT_Bt:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,              gz,dFS_HOMOGENEOUS, gx,dFS_HOMOGENEOUS,NULL);dCHK(err);
    break;
  case STOKES_MULT_B:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, gx,dFS_HOMOGENEOUS,NULL,               NULL,              gz,dFS_HOMOGENEOUS);dCHK(err);
    break;
  default: dERROR(stk->comm,PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmode");
  }
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*du)[9],(*dv)[9],*p,*q;
    dInt Q;
    struct StokesStore *stash;
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    switch (mmode) {
    case STOKES_MULT_A:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,&du,NULL,&dv, NULL,NULL,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {
        dScalar Du[6],Dv[6];
        dTensorSymCompress3(du[i],Du);
        StokesPointwiseJacobian_A(&stk->scase->rheo,&stash[i],jw[i],Du,Dv);
        dTensorSymUncompress3(Dv,dv[i]);
      }
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,dv, NULL,NULL);dCHK(err);
      break;
    case STOKES_MULT_Bt:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,NULL,NULL,&dv, &p,NULL,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {
        dScalar Dv[6];
        StokesPointwiseJacobian_Bt(jw[i],p[i],Dv);
        dTensorSymUncompress3(Dv,dv[i]);
      }
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,dv, NULL,NULL);dCHK(err);
      break;
    case STOKES_MULT_B:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,&du,NULL,NULL, NULL,NULL,&q,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {
        dScalar Du[6];
        dTensorSymCompress3(du[i],Du);
        StokesPointwiseJacobian_B(jw[i],Du,&q[i]);
      }
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,NULL, q,NULL);dCHK(err);
      break;
    default: dERROR(stk->comm,PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmode");
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = PetscLogEventEnd(LOG_StokesShellMult,A,gx,gy,gz);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesJacobianAssemble_Velocity(Stokes stk,Mat Ap,Vec Mdiag,Vec gx)
{
  dRulesetIterator iter;
  Vec Coords,gxu;
  dScalar *Kflat;
  dErr err;

  dFunctionBegin;
  err = VecZeroEntries(Mdiag);dCHK(err);
  err = StokesExtractGlobalSplit(stk,gx,&gxu,NULL);dCHK(err);
  err = dFSGetGeometryVectorExpanded(stk->fsu,&Coords);dCHK(err);
  err = StokesGetRegionIterator(stk,EVAL_JACOBIAN,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, gxu,dFS_INHOMOGENEOUS,Mdiag,dFS_HOMOGENEOUS, NULL,NULL);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSplit(iter, NULL,NULL,NULL, NULL,&Kflat,NULL, NULL,NULL,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw,*interp_flat,*deriv_flat;
    const dInt *rowcol;
    dScalar (*x)[3],(*dx)[3][3],(*du)[9],(*v)[3];
    dInt Q,P;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,&du,&v,NULL, NULL,NULL,NULL,NULL);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat, NULL,NULL,NULL,NULL);dCHK(err);
    {                           /* Scope so that we can declare new VLA pointers for convenient assembly */
      const dReal (*interp)[P] = (const dReal(*)[P])interp_flat;
      const dReal (*deriv)[P][3] = (const dReal(*)[P][3])deriv_flat;
      dScalar (*K)[3][P][3] = (dScalar(*)[3][P][3])Kflat;
      err = PetscMemzero(K,P*3*P*3*sizeof(K[0][0][0][0]));dCHK(err);
      for (dInt q=0; q<Q; q++) {
        struct StokesStore store;
        dScalar Dusym[6];
        dTensorSymCompress3(du[q],Dusym);
        StokesPointwiseComputeStore(&stk->scase->rheo,x[q],Dusym,&store);
        for (dInt j=0; j<P; j++) { /* trial functions */
          for (dInt fj=0; fj<3; fj++) {
            dScalar duu[3][3] = {{0},{0},{0}},dv[3][3],Duusym[6],Dvsym[6];
            duu[fj][0] = deriv[q][j][0];
            duu[fj][1] = deriv[q][j][1];
            duu[fj][2] = deriv[q][j][2];
            dTensorSymCompress3(&duu[0][0],Duusym);
            StokesPointwiseJacobian_A(&stk->scase->rheo,&store,jw[q],Duusym,Dvsym);
            dTensorSymUncompress3(Dvsym,&dv[0][0]);
            for (dInt i=0; i<P; i++) {
              for (dInt fi=0; fi<3; fi++) {
                K[i][fi][j][fj] += (+ deriv[q][i][0] * dv[fi][0]
                                    + deriv[q][i][1] * dv[fi][1]
                                    + deriv[q][i][2] * dv[fi][2]);
              }
            }
          }
        }
      }
      err = dFSMatSetValuesBlockedExpanded(stk->fsu,Ap,P,rowcol,P,rowcol,&K[0][0][0][0],ADD_VALUES);dCHK(err);
      for (dInt i=0; i<P; i++) {
        dScalar Mentry = 0;
        for (dInt q=0; q<Q; q++) Mentry += interp[q][i] * jw[q] * interp[q][i]; /* Integrate the diagonal entry over this element */
        v[i][0] += Mentry;
        v[i][1] += Mentry;
        v[i][2] += Mentry;
      }
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, v,NULL, NULL,NULL);dCHK(err);
    err = dRulesetIteratorRestorePatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat, NULL,NULL,NULL,NULL);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesJacobianAssemble_Pressure(Stokes stk,Mat D,Mat Daux,Vec gx)
{
  dRulesetIterator iter;
  Vec              Coords,gxu;
  dScalar          *Kflat,*Kflat_aux;
  const dInt       *Ksizes;
  dErr             err;

  dFunctionBegin;
  /* It might seem weird to be getting velocity in the pressure assembly.  The reason is that this preconditioner
  * (indeed the entire problem) is always linear in pressure.  It \e might be nonlinear in velocity. */
  err = StokesExtractGlobalSplit(stk,gx,&gxu,NULL);dCHK(err);
  err = dFSGetGeometryVectorExpanded(stk->fsu,&Coords);dCHK(err);
  err = StokesGetRegionIterator(stk,EVAL_JACOBIAN,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, gxu,dFS_INHOMOGENEOUS,NULL, NULL,NULL);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSplit(iter, NULL,NULL,NULL, NULL,NULL,NULL, NULL,NULL,&Kflat);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSizes(iter,NULL,NULL,&Ksizes);dCHK(err);
  err = dMallocA(Ksizes[8],&Kflat_aux);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw,*interp_flat,*deriv_flat;
    const dInt *rowcol;
    dScalar (*x)[3],(*dx)[3][3],(*du)[9];
    dInt Q,P;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,&du,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat);dCHK(err);
    {
      const dReal (*interp)[P] = (const dReal(*)[P])interp_flat;
      const dReal (*deriv)[P][3] = (const dReal(*)[P][3])deriv_flat;
      dScalar (*K)[P] = (dScalar(*)[P])Kflat,(*Ka)[P] = (dScalar(*)[P])Kflat_aux;
      err = PetscMemzero(K,P*P*sizeof(K[0][0]));dCHK(err);
      err = PetscMemzero(Ka,P*P*sizeof(K[0][0]));dCHK(err);
      for (dInt q=0; q<Q; q++) {
        struct StokesStore store;
        dScalar Dusym[6];
        dTensorSymCompress3(du[q],Dusym);
        StokesPointwiseComputeStore(&stk->scase->rheo,x[q],Dusym,&store);dCHK(err);
        for (dInt j=0; j<P; j++) { /* trial functions */
          for (dInt i=0; i<P; i++) {
            /* Scaled mass matrx */
            K[i][j] += interp[q][i] * jw[q] * (1./store.eta) * interp[q][j];
            /* Neumann Laplacian */
            Ka[i][j] += (+ deriv[q][i][0] * jw[q] * deriv[q][j][0]
                         + deriv[q][i][1] * jw[q] * deriv[q][j][1]
                         + deriv[q][i][2] * jw[q] * deriv[q][j][2]);
          }
        }
      }
      err = dFSMatSetValuesBlockedExpanded(stk->fsp,D,P,rowcol,P,rowcol,&K[0][0],ADD_VALUES);dCHK(err);
      if (Daux) {err = dFSMatSetValuesBlockedExpanded(stk->fsp,Daux,P,rowcol,P,rowcol,&Ka[0][0],ADD_VALUES);dCHK(err);}
    }
    err = dRulesetIteratorRestorePatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat, NULL,NULL,NULL,NULL);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = dFree(Kflat_aux);dCHK(err);
  dFunctionReturn(0);
}


static dErr StokesJacobian(SNES dUNUSED snes,Vec gx,Mat *J,Mat *Jp,MatStructure *structure,void *ctx)
{
  Stokes stk = ctx;
  dErr err;
  Mat A,D,Daux;
  Vec Mdiag;

  dFunctionBegin;
  err = MatGetLocalSubMatrix(*Jp,stk->lublock,stk->lublock,&A);dCHK(err);
  err = MatGetLocalSubMatrix(*Jp,stk->lpblock,stk->lpblock,&D);dCHK(err);
  err = PetscObjectQuery((dObject)D,"LSC_M_diag",(dObject*)&Mdiag);dCHK(err);
  err = PetscObjectQuery((dObject)D,"LSC_L",(dObject*)&Daux);dCHK(err);
  err = MatZeroEntries(*Jp);dCHK(err);
  if (Daux) {err = MatZeroEntries(Daux);dCHK(err);}
  err = StokesJacobianAssemble_Velocity(stk,A,Mdiag,gx);dCHK(err);
  err = StokesJacobianAssemble_Pressure(stk,D,Daux,gx);dCHK(err);
  if (Daux) {
    err = MatAssemblyBegin(Daux,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd  (Daux,MAT_FINAL_ASSEMBLY);dCHK(err);
  }
  err = MatRestoreLocalSubMatrix(*Jp,stk->lublock,stk->lublock,&A);dCHK(err);
  err = MatRestoreLocalSubMatrix(*Jp,stk->lpblock,stk->lpblock,&D);dCHK(err);

  /* MatNest calls assembly on the constituent pieces */
  err = MatAssemblyBegin(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  if (*J != *Jp) {
    err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  }
  *structure = SAME_NONZERO_PATTERN;
  dFunctionReturn(0);
}

static dErr StokesErrorNorms(Stokes stk,Vec gx,dReal errorNorms[3],dReal gerrorNorms[3],dReal perrorNorms[3])
{
  dErr             err;
  Vec              Coords,gxu,gxp;
  PetscScalar      volume = 0,pressureshift = 0;
  dRulesetIterator iter;

  dFunctionBegin;
  err = dNormsStart(errorNorms,gerrorNorms);dCHK(err);
  err = dNormsStart(perrorNorms,NULL);dCHK(err);
  err = StokesExtractGlobalSplit(stk,gx,&gxu,&gxp);dCHK(err);
  err = StokesGetRegionIterator(stk,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(stk->fsu,&Coords);dCHK(err);
  if (stk->alldirichlet) { // Do a volume integral of the exact solution to that we can remove the constant pressure mode
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,NULL, gxp,dFS_INHOMOGENEOUS,NULL);dCHK(err);
    while (dRulesetIteratorHasPatch(iter)) {
      const dReal *jw;
      const dScalar (*x)[3],*p;
      dInt Q;
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,NULL,NULL,NULL, NULL,NULL,NULL,NULL, &p,NULL,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {
        dScalar uu[3],duu[9],pp[1],dpp[3];
        err = stk->scase->solution(stk->scase,x[i],uu,duu,pp,dpp);dCHK(err);
        volume += jw[i];
        pressureshift += (pp[0] - p[i]) * jw[i]; // The computed pressure sum is zero, but the continuous integral may not be
      }
      err = dRulesetIteratorNextPatch(iter);dCHK(err);
    }
    err = dRulesetIteratorFinish(iter);dCHK(err);
    pressureshift /= volume;
  }
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, gxu,dFS_INHOMOGENEOUS,NULL, gxp,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw;
    const dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*p)[1];
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, (const dScalar**)&u,(const dScalar**)&du,NULL,NULL, (const dScalar**)&p,NULL,NULL,NULL);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[3],duu[9],pp[1],dpp[3];
      err = stk->scase->solution(stk->scase,x[i],uu,duu,pp,dpp);dCHK(err);
      pp[0] -= pressureshift;
      err = dNormsUpdate(errorNorms,gerrorNorms,jw[i],3,uu,u[i],duu,du[i]);dCHK(err);
      err = dNormsUpdate(perrorNorms,NULL,jw[i],1,pp,p[i],NULL,NULL);dCHK(err);
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = dNormsFinish(errorNorms,gerrorNorms);dCHK(err);
  err = dNormsFinish(perrorNorms,NULL);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesGetSolutionField_All(Stokes stk,dFS fs,dBool isvel,Vec *insoln)
{
  Vec      sol,xc,cvecg,cvec;
  dScalar *x;
  const dScalar *coords;
  dInt     n,bs;
  dErr     err;

  dFunctionBegin;
  *insoln = 0;
  err = dFSCreateGlobalVector(fs,&sol);dCHK(err);
  err = VecDohpGetClosure(sol,&xc);dCHK(err);
  err = dFSGetNodalCoordinatesGlobal(fs,&cvecg);dCHK(err);
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
    dScalar u_unused[3],p_unused[1],du_unused[3*3],dp_unused[3];
    /* if \a isvel then \a x is the velocity field, otherwise it is the pressure field */
    err = stk->scase->solution(stk->scase,&coords[3*i],isvel ? &x[i*bs] : u_unused,du_unused,isvel ? p_unused : &x[i*bs],dp_unused);dCHK(err);
    /* printf("Node %3d: coords %+8f %+8f %+8f   exact %+8f %+8f %+8f\n",i,coords[3*i],coords[3*i+1],coords[3*i+2],x[3*i],x[3*i+1],x[3*i+2]); */
  }
  err = VecRestoreArray(xc,&x);dCHK(err);
  err = VecRestoreArrayRead(cvec,&coords);dCHK(err);
  err = VecDohpRestoreClosure(cvecg,&cvec);dCHK(err);
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
  err = VecDestroy(&solu);dCHK(err);
  err = VecDestroy(&solp);dCHK(err);
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
  err = VecDestroy(&r);dCHK(err);
  dFunctionReturn(0);
}

static dErr CheckNullSpace(SNES snes,Vec residual,dBool compute_explicit)
{
  Mat          mffd,J,Jp;
  dBool        isnull;
  Vec          U,F;
  MatStructure mstruct;
  MatNullSpace matnull;
  KSP          ksp;
  dErr         err;

  dFunctionBegin;
  err = SNESGetKSP(snes,&ksp);dCHK(err);
  err = KSPGetNullSpace(ksp,&matnull);dCHK(err);
  err = MatCreateSNESMF(snes,&mffd);dCHK(err);
  err = MatSetFromOptions(mffd);dCHK(err);
  {
    err = VecDuplicate(residual,&U);dCHK(err);
    err = VecDuplicate(residual,&F);dCHK(err);
  }
  err = SNESGetJacobian(snes,&J,&Jp,NULL,NULL);dCHK(err);
  err = VecSet(U,0);dCHK(err);
  err = SNESComputeFunction(snes,U,F);dCHK(err); /* Need base for MFFD */
  err = MatMFFDSetBase(mffd,U,F);dCHK(err);
  err = MatNullSpaceTest(matnull,mffd,&isnull);dCHK(err);
  if (!isnull) dERROR(PETSC_COMM_SELF,1,"Vector is not in the null space of the MFFD operator");dCHK(err);
  err = MatNullSpaceTest(matnull,J,&isnull);dCHK(err);
  if (!isnull) dERROR(PETSC_COMM_SELF,1,"Vector is not in the null space of J");dCHK(err);
  err = SNESComputeJacobian(snes,U,&J,&Jp,&mstruct);dCHK(err); /* To assemble blocks of Jp */
  err = MatNullSpaceTest(matnull,Jp,&isnull);dCHK(err);
  // At present, Jp intentionally contains an auxilliary matrix in the (p,p) block. It does not have the same null space as the Jacobian so we disable the error below.
  if (false && !isnull) dERROR(PETSC_COMM_SELF,1,"Vector is not in the null space of Jp");dCHK(err);
  err = MatDestroy(&mffd);dCHK(err);
  if (compute_explicit) {
    Mat expmat,expmat_fd;
    dInt m,n;
    dBool contour = dFALSE;
    err = MatGetLocalSize(J,&m,&n);dCHK(err);
    err = MatComputeExplicitOperator(J,&expmat);dCHK(err);
    err = MatDuplicate(expmat,MAT_DO_NOT_COPY_VALUES,&expmat_fd);dCHK(err);
    err = SNESComputeJacobianDefault(snes,U,&expmat_fd,&expmat_fd,&mstruct,NULL);dCHK(err);
    err = MatSetOptionsPrefix(expmat,"explicit_");dCHK(err);
    err = MatSetOptionsPrefix(expmat_fd,"explicit_fd_");dCHK(err);
    err = MatSetFromOptions(expmat);dCHK(err);
    err = MatSetFromOptions(expmat_fd);dCHK(err);

    err = PetscOptionsGetBool(NULL,"-mat_view_contour",&contour,NULL);dCHK(err);
    if (contour) {err = PetscViewerPushFormat(PETSC_VIEWER_DRAW_WORLD,PETSC_VIEWER_DRAW_CONTOUR);dCHK(err);}
    {
      dBool flg = dFALSE;
      err = PetscOptionsGetBool(NULL,"-explicit_mat_view",&flg,NULL);dCHK(err);
      if (flg) {
        err = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD,"###  Explicit matrix using mat-free implementation of J\n");dCHK(err);
        err = MatView(expmat,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
      }
      flg = dFALSE;
      err = PetscOptionsGetBool(NULL,"-explicit_mat_view_draw",&flg,NULL);dCHK(err);
      if (flg) {err = MatView(expmat,PETSC_VIEWER_DRAW_WORLD);dCHK(err);}
    }

    {
      dBool flg = dFALSE;
      err = PetscOptionsGetBool(NULL,"-explicit_fd_mat_view",&flg,NULL);dCHK(err);
      if (flg) {
        err = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD,"###  Explicit matrix using FD\n");dCHK(err);
        err = MatView(expmat_fd,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
      }
      flg = dFALSE;
      err = PetscOptionsGetBool(NULL,"-explicit_fd_mat_view_draw",&flg,NULL);dCHK(err);
      if (flg) {err = MatView(expmat_fd,PETSC_VIEWER_DRAW_WORLD);dCHK(err);}
    }

    err = MatAXPY(expmat,-1,expmat_fd,SAME_NONZERO_PATTERN);dCHK(err);
    {
      dBool flg = dFALSE;
      err = PetscOptionsGetBool(NULL,"-explicit_diff_mat_view",&flg,NULL);dCHK(err);
      if (flg) {
        err = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD,"###  Difference between mat-free implementation of J and FD\n");dCHK(err);
        err = MatView(expmat,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
      }
      flg = dFALSE;
      err = PetscOptionsGetBool(NULL,"-explicit_diff_mat_view_draw",&flg,NULL);dCHK(err);
      if (flg) {err = MatView(expmat,PETSC_VIEWER_DRAW_WORLD);dCHK(err);}
    }
    if (contour) {err = PetscViewerPopFormat(PETSC_VIEWER_DRAW_WORLD);dCHK(err);}
    err = MatDestroy(&expmat);dCHK(err);
    err = MatDestroy(&expmat_fd);dCHK(err);
  }
  err = VecDestroy(&U);dCHK(err);
  err = VecDestroy(&F);dCHK(err);
  dFunctionReturn(0);
}


int main(int argc,char *argv[])
{
  Stokes stk;
  MPI_Comm comm;
  Mat J,Jp;
  MatFDColoring fdcolor = NULL;
  Vec r,x,soln = NULL;
  SNES snes;
  dBool check_error,check_null,compute_explicit,use_jblock,viewdhm;
  dErr err;

  err = dInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  err = PetscLogEventRegister("StokesShellMult",MAT_CLASSID,&LOG_StokesShellMult);dCHK(err);

  err = StokesCaseRegisterAll();dCHK(err);
  err = StokesCreate(comm,&stk);dCHK(err);
  err = StokesSetFromOptions(stk);dCHK(err);

  err = VecDuplicate(stk->gpacked,&r);dCHK(err);
  err = VecDuplicate(r,&x);dCHK(err);

  err = PetscOptionsBegin(stk->comm,NULL,"Stokes solver options",__FILE__);dCHK(err); {
    check_error = stk->scase->reality ? dFALSE : dTRUE;
    err = PetscOptionsBool("-check_error","Compute errors","",check_error,&check_error,NULL);dCHK(err);
    err = PetscOptionsBool("-use_jblock","Use blocks to apply Jacobian instead of unified (more efficient) version","",use_jblock=dFALSE,&use_jblock,NULL);dCHK(err);
    err = PetscOptionsBool("-viewdhm","View the solution","",viewdhm=dFALSE,&viewdhm,NULL);dCHK(err);
    err = PetscOptionsBool("-check_null","Check that constant pressure really is in the null space","",check_null=dFALSE,&check_null,NULL);dCHK(err);
    if (check_null) {
      err = PetscOptionsBool("-compute_explicit","Compute explicit Jacobian (only very small sizes)","",compute_explicit=dFALSE,&compute_explicit,NULL);dCHK(err);
    }
  } err = PetscOptionsEnd();dCHK(err);
  err = StokesGetMatrices(stk,use_jblock,&J,&Jp);dCHK(err);
  err = SNESCreate(comm,&snes);dCHK(err);
  err = SNESSetFunction(snes,r,StokesFunction,stk);dCHK(err);
  switch (3) {
    case 1:
      err = SNESSetJacobian(snes,J,Jp,SNESComputeJacobianDefault,stk);dCHK(err); break;
    case 2: {
      ISColoring iscolor;
      err = MatGetColoring(Jp,MATCOLORINGID,&iscolor);dCHK(err);
      err = MatFDColoringCreate(Jp,iscolor,&fdcolor);dCHK(err);
      err = ISColoringDestroy(&iscolor);dCHK(err);
      err = MatFDColoringSetFunction(fdcolor,(PetscErrorCode(*)(void))StokesFunction,stk);dCHK(err);
      err = MatFDColoringSetFromOptions(fdcolor);dCHK(err);
      err = SNESSetJacobian(snes,J,Jp,SNESComputeJacobianDefaultColor,fdcolor);dCHK(err);
    } break;
    case 3:
      err = SNESSetJacobian(snes,J,Jp,StokesJacobian,stk);dCHK(err);
      break;
    default: dERROR(PETSC_COMM_SELF,1,"Not supported");
  }
  err = SNESSetFromOptions(snes);dCHK(err);
  {
    KSP    ksp;
    PC     pc;

    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPGetPC(ksp,&pc);dCHK(err);
    err = PCFieldSplitSetIS(pc,"u",stk->ublock);dCHK(err);
    err = PCFieldSplitSetIS(pc,"p",stk->pblock);dCHK(err);
  }
  err = StokesGetSolutionVector(stk,&soln);dCHK(err);
  if (!stk->scase->reality) {
    dReal nrm;
    MatStructure mstruct;
    Vec b;
    err = VecDuplicate(x,&b);dCHK(err);
    err = VecZeroEntries(x);dCHK(err);
    err = SNESComputeFunction(snes,x,b);dCHK(err); /* -f */
    err = SNESComputeFunction(snes,soln,r);dCHK(err);
    err = VecNorm(r,NORM_2,&nrm);dCHK(err);
    err = dPrintf(comm,"Norm of discrete residual for exact solution %g\n",nrm);dCHK(err);
    err = SNESComputeJacobian(snes,soln,&J,&Jp,&mstruct);dCHK(err);
    err = MatMult(J,soln,r);dCHK(err);
    err = VecAXPY(r,1,b);dCHK(err); /* Jx - f */
    err = VecNorm(r,NORM_2,&nrm);dCHK(err);
    err = dPrintf(comm,"Norm of discrete linear residual at exact solution %g\n",nrm);dCHK(err);
    err = VecDestroy(&b);dCHK(err);
  }

  if (stk->alldirichlet) {                             /* Set null space */
    KSP ksp;
    MatNullSpace matnull;
    err = StokesGetNullSpace(stk,&matnull);dCHK(err);
    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPSetNullSpace(ksp,matnull);dCHK(err);
    if (soln) {err = MatNullSpaceRemove(matnull,soln,NULL);dCHK(err);}
    err = MatNullSpaceDestroy(&matnull);dCHK(err);
  }
  if (check_null) {
    err = CheckNullSpace(snes,r,compute_explicit);dCHK(err);
  }
  err = VecZeroEntries(r);dCHK(err);
  err = VecZeroEntries(x);dCHK(err);
  err = SNESSolve(snes,NULL,x);dCHK(err); /* ###  SOLVE  ### */
  if (stk->alldirichlet) {
    MatNullSpace matnull;
    KSP ksp;
    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPGetNullSpace(ksp,&matnull);dCHK(err); /* does not reference */
    err = MatNullSpaceRemove(matnull,x,NULL);dCHK(err);
  }
  if (check_error) {
    dReal anorm[3],inorm[3],enorm[3],gnorm[3],epnorm[3];
    err = StokesErrorNorms(stk,x,enorm,gnorm,epnorm);dCHK(err);
    err = dNormsAlgebraicScaled(anorm,r);dCHK(err);
    err = VecWAXPY(r,-1,soln,x);dCHK(err);
    err = dNormsAlgebraicScaled(inorm,r);dCHK(err);
    err = dPrintf(comm,"Algebraic residual        |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",anorm[0],anorm[1],anorm[2]);dCHK(err);
    err = dPrintf(comm,"Interpolation residual    |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",inorm[0],inorm[1],inorm[2]);dCHK(err);
    err = dPrintf(comm,"Pointwise velocity error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",enorm[0],enorm[1],enorm[2]);dCHK(err);
    err = dPrintf(comm,"Pointwise gradient error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",gnorm[0],gnorm[1],gnorm[2]);dCHK(err);
    err = dPrintf(comm,"Pointwise pressure error  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",epnorm[0],epnorm[1],epnorm[2]);dCHK(err);
  }
  if (viewdhm) {
    Vec Xu,Xp;
    dViewer view;
    err = PetscViewerCreate(comm,&view);dCHK(err);
    err = PetscViewerSetType(view,PETSCVIEWERDHM);dCHK(err);
    err = PetscViewerFileSetName(view,"stokes.dhm");dCHK(err);
    err = PetscViewerFileSetMode(view,FILE_MODE_WRITE);dCHK(err);
    err = StokesExtractGlobalSplit(stk,x,&Xu,&Xp);dCHK(err);
    err = dFSDirichletProject(stk->fsu,Xu,dFS_INHOMOGENEOUS);dCHK(err);
    err = dFSDirichletProject(stk->fsp,Xp,dFS_INHOMOGENEOUS);dCHK(err);
    err = VecView(Xu,view);dCHK(err);
    err = VecView(Xp,view);dCHK(err);
    err = PetscViewerDestroy(&view);dCHK(err);
  }

  err = VecDestroy(&r);dCHK(err);
  err = VecDestroy(&x);dCHK(err);
  err = VecDestroy(&soln);dCHK(err);
  err = SNESDestroy(&snes);dCHK(err);
  err = MatFDColoringDestroy(&fdcolor);dCHK(err);
  if (J != Jp) {err = MatDestroy(&J);dCHK(err);}
  err = MatDestroy(&Jp);dCHK(err);
  err = StokesDestroy(&stk);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
