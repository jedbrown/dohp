static const char help[] = "Solve viscous flow coupled to a heat transport problem using dual order elements.\n"
  "The model problem is\n"
  "  -div(eta Du) + grad(p) = f\n"
  "                  div(u) = g\n"
  "    div(u T) - eta Du:Du = h\n"
  "where\n"
  "  D is the symmetric gradient operator\n"
  "  eta(gamma,T) = B(T) (0.5*eps^2 + gamma)^{(p-2)/2}\n"
  "  gamma = Du : Du/2\n"
  "  B(T) = B_0 exp(Q/(n R T))\n"
  "The weak form is\n"
  "  int_Omega eta Dv:Du - p div(v) - q div(u) - v.f - q g -  = 0\n"
  "with Jacobian\n"
  "  int_Omega eta Dv:Du + eta' (Dv:Dw)(Dw:Du) - p div(v) - q div(u) = 0\n"
  "The problem is linear for p=2, an incompressible for g=0\n\n";

#include <petscts.h>
#include <dohpstring.h>
#include <dohpviewer.h>

#include "vhtimpl.h"

PetscFList VHTCaseList = NULL;

#define VHTCaseType char*

dErr VHTCaseRegister(const char *name,VHTCaseCreateFunction screate)
{
  dErr err;
  dFunctionBegin;
  err = PetscFListAdd(&VHTCaseList,name,"",(void(*)(void))screate);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseFind(const char *name,VHTCaseCreateFunction *screate)
{
  dErr err;

  dFunctionBegin;
  err = PetscFListFind(VHTCaseList,PETSC_COMM_WORLD,name,PETSC_FALSE,(void(**)(void))screate);dCHK(err);
  if (!*screate) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"VHT Case \"%s\" could not be found",name);
  dFunctionReturn(0);
}
static dErr VHTCaseSetType(VHTCase scase,const VHTCaseType type)
{
  dErr err;
  VHTCaseCreateFunction f;

  dFunctionBegin;
  err = VHTCaseFind(type,&f);dCHK(err);
  err = (*f)(scase);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseSetFromOptions(VHTCase scase)
{
  struct VHTRheology *rheo = &scase->rheo;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsBegin(scase->comm,NULL,"VHTCase_Exact options",__FILE__);dCHK(err); {
    err = PetscOptionsReal("-rheo_B0","Rate factor (rheology)","",rheo->B0,&rheo->B0,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_R","Ideal gas constant","",rheo->R,&rheo->R,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_Q","Activation Energy","",rheo->Q,&rheo->Q,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_eps","Regularization (rheology)","",rheo->eps,&rheo->eps,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_p","Power p=1+1/n where n is Glen exponent","",rheo->pe,&rheo->pe,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_kappa0","Thermal conductivity","",rheo->kappa0,&rheo->kappa0,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_kappa1","Thermal conductivity","",rheo->kappa1,&rheo->kappa1,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_T0","Reference temperature (corresponds to enthalpy=0)","",rheo->T0,&rheo->T0,NULL);dCHK(err);
    err = PetscOptionsReal("-gravity","Nondimensional gravitational force","",scase->gravity,&scase->gravity,NULL);dCHK(err);
    if (scase->setfromoptions) {err = (*scase->setfromoptions)(scase);dCHK(err);}
  } err = PetscOptionsEnd();dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseDestroy(VHTCase *scase)
{
  dErr err;

  dFunctionBegin;
  if ((*scase)->destroy) {err = ((*scase)->destroy)(*scase);dCHK(err);}
  err = dFree(*scase);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseRegisterAll(void)
{
  dErr err;

  dFunctionBegin;
  err = VHTCaseRegisterAll_Exact();dCHK(err);
  dFunctionReturn(0);
}

#define ALEN(a) ((dInt)(sizeof(a)/sizeof(a)[0]))
static PetscLogEvent LOG_VHTShellMult;

static dErr VHTGetNullSpace(VHT vht,MatNullSpace *matnull);
static dErr MatMultXIorA_VHT_stokes(Mat A,Vec gx,Vec gy,Vec gz,InsertMode,VHTMultMode);
static dErr MatMult_Nest_VHT_all(Mat J,Vec gx,Vec gy);
static dErr MatMult_VHT_uu(Mat A,Vec gx,Vec gy) {return MatMultXIorA_VHT_stokes(A,gx,gy,NULL,INSERT_VALUES,VHT_MULT_UU);}
static dErr MatMult_VHT_up(Mat A,Vec gx,Vec gy) {return MatMultXIorA_VHT_stokes(A,gx,gy,NULL,INSERT_VALUES,VHT_MULT_UP);}
static dErr MatMult_VHT_pu(Mat A,Vec gx,Vec gy) {return MatMultXIorA_VHT_stokes(A,gx,gy,NULL,INSERT_VALUES,VHT_MULT_PU);}
static dErr MatMultAdd_VHT_uu(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultXIorA_VHT_stokes(A,gx,gy,gz,ADD_VALUES,VHT_MULT_UU);}
static dErr MatMultAdd_VHT_up(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultXIorA_VHT_stokes(A,gx,gy,gz,ADD_VALUES,VHT_MULT_UP);}
static dErr MatMultAdd_VHT_pu(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultXIorA_VHT_stokes(A,gx,gy,gz,ADD_VALUES,VHT_MULT_PU);}
static dErr MatMult_VHT_ee(Mat A,Vec gx,Vec gy);
static dErr MatGetVecs_VHT_stokes(Mat,Vec*,Vec*);
static dErr MatGetVecs_VHT_ee(Mat,Vec*,Vec*);

static dErr VHTCreate(MPI_Comm comm,VHT *invht)
{
  VHT vht;
  dErr err;

  dFunctionBegin;
  *invht = 0;
  err = dNew(struct _n_VHT,&vht);dCHK(err);
  vht->comm = comm;

  vht->velocityBDeg  = 3;
  vht->pressureCodim = 1;
  vht->enthalpyBDeg  = 3;
  vht->dirichlet[0]  = 100;
  vht->dirichlet[1]  = 200;
  vht->dirichlet[2]  = 300;
  vht->alldirichlet  = dTRUE;
  vht->function_qmethod = dQUADRATURE_METHOD_FAST;
  vht->jacobian_qmethod = dQUADRATURE_METHOD_SPARSE;

  err = dCalloc(sizeof(*vht->scase),&vht->scase);dCHK(err);

  vht->scase->rheo.B0     = 1;
  vht->scase->rheo.R      = 1;
  vht->scase->rheo.Q      = 1;
  vht->scase->rheo.eps    = 1;
  vht->scase->rheo.pe     = 2;
  vht->scase->rheo.kappa0 = 10;
  vht->scase->rheo.kappa1 = 5;
  vht->scase->rheo.T0     = 10;

  *invht = vht;
  dFunctionReturn(0);
}

static dErr MatGetVecs_VHT_stokes(Mat A,Vec *x,Vec *y)
{
  VHT vht;
  dInt m,n,nu,np;
  dErr err;

  dFunctionBegin;
  err = MatShellGetContext(A,(void**)&vht);dCHK(err);
  err = MatGetLocalSize(A,&m,&n);dCHK(err);
  err = VecGetLocalSize(vht->gvelocity,&nu);dCHK(err);
  err = VecGetLocalSize(vht->gpressure,&np);dCHK(err);
  if (nu==np) dERROR(PETSC_COMM_SELF,1,"Degenerate case, don't know which space to copy");
  if (x) {
    if (n == nu) {
      err = VecDuplicate(vht->gvelocity,x);dCHK(err);
    } else if (n == np) {
      err = VecDuplicate(vht->gpressure,x);dCHK(err);
    } else dERROR(PETSC_COMM_SELF,1,"sizes do not agree with either space");
  }
  if (y) {
    if (n == nu) {
      err = VecDuplicate(vht->gvelocity,y);dCHK(err);
    } else if (n == np) {
      err = VecDuplicate(vht->gpressure,y);dCHK(err);
    } else dERROR(PETSC_COMM_SELF,1,"sizes do not agree with either space");
  }
  dFunctionReturn(0);
}

static dErr MatGetVecs_VHT_ee(Mat A,Vec *x,Vec *y)
{
  VHT vht;
  dErr err;

  dFunctionBegin;
  err = MatShellGetContext(A,(void**)&vht);dCHK(err);
  if (x) {err = VecDuplicate(vht->genthalpy,x);dCHK(err);}
  if (y) {err = VecDuplicate(vht->genthalpy,y);dCHK(err);}
  dFunctionReturn(0);
}

static dErr VHTSetFromOptions(VHT vht)
{
  char scasename[256] = "Exact0";
  dMesh mesh;
  dFS fsu,fsp,fse;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag dutag,dptag,detag;
  dErr err;

  dFunctionBegin;
  err = dStrcpyS(vht->mattype_Buu,sizeof(vht->mattype_Buu),MATBAIJ);dCHK(err);
  err = dStrcpyS(vht->mattype_Bpp,sizeof(vht->mattype_Bpp),MATAIJ);dCHK(err);
  err = dStrcpyS(vht->mattype_Bee,sizeof(vht->mattype_Bee),MATAIJ);dCHK(err);
  err = PetscOptionsBegin(vht->comm,NULL,"Viscous Heat Transport options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-vht_u_bdeg","Constant isotropic degree to use for velocity","",vht->velocityBDeg,&vht->velocityBDeg,NULL);dCHK(err);
    err = PetscOptionsInt("-vht_p_codim","Reduce pressure space by this factor","",vht->pressureCodim,&vht->pressureCodim,NULL);dCHK(err);
    err = PetscOptionsInt("-vht_e_bdeg","Constant isotropic degree to use for enthalpy","",vht->enthalpyBDeg,&vht->enthalpyBDeg,NULL);dCHK(err);
    err = PetscOptionsBool("-vht_cardinal_mass","Assemble diagonal mass matrix","",vht->cardinalMass,&vht->cardinalMass,NULL);dCHK(err);
    err = PetscOptionsList("-vht_Buu_mat_type","Matrix type for velocity-velocity operator","",MatList,vht->mattype_Buu,vht->mattype_Buu,sizeof(vht->mattype_Buu),NULL);dCHK(err);
    err = PetscOptionsList("-vht_Bpp_mat_type","Matrix type for pressure-pressure operator","",MatList,vht->mattype_Bpp,vht->mattype_Bpp,sizeof(vht->mattype_Bpp),NULL);dCHK(err);
    err = PetscOptionsList("-vht_Bee_mat_type","Matrix type for enthalpy-enthalpy operator","",MatList,vht->mattype_Bee,vht->mattype_Bee,sizeof(vht->mattype_Bee),NULL);dCHK(err);
    err = PetscOptionsEnum("-vht_f_qmethod","Quadrature method for residual evaluation/matrix-free","",dQuadratureMethods,(PetscEnum)vht->function_qmethod,(PetscEnum*)&vht->function_qmethod,NULL);dCHK(err);
    err = PetscOptionsEnum("-vht_jac_qmethod","Quadrature to use for Jacobian assembly","",dQuadratureMethods,(PetscEnum)vht->jacobian_qmethod,(PetscEnum*)&vht->jacobian_qmethod,NULL);dCHK(err);
    {
      dBool flg; dInt n = ALEN(vht->dirichlet);
      err = PetscOptionsIntArray("-dirichlet","List of boundary sets on which to impose Dirichlet conditions","",vht->dirichlet,&n,&flg);dCHK(err);
      if (flg) {
        for (dInt i=n; i<ALEN(vht->dirichlet); i++) vht->dirichlet[i] = 0; /* Clear out any leftover values */
        if (n < 3) vht->alldirichlet = dFALSE;                             /* @bug More work to determine independent of the mesh whether all the boundaries are Dirichlet */
      }
    }
    err = PetscOptionsList("-vht_case","Which sort of case to run","",VHTCaseList,scasename,scasename,sizeof(scasename),NULL);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);

  err = dMeshCreate(vht->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);
  err = dMeshGetRoot(mesh,&domain);dCHK(err); /* Need a taggable set */
  err = dMeshSetDuplicateEntsOnly(mesh,domain,&domain);dCHK(err);
  err = PetscObjectSetName((PetscObject)mesh,"dMesh_0");dCHK(err);

  err = dJacobiCreate(vht->comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);

  err = dMeshCreateRuleTagIsotropic(mesh,domain,"vht_efs_velocity_degree",vht->velocityBDeg,&dutag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,"vht_efs_pressure_degree",vht->velocityBDeg-vht->pressureCodim,&dptag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,"vht_efs_enthalpy_degree",vht->enthalpyBDeg,&detag);dCHK(err);

  err = dFSCreate(vht->comm,&fsu);dCHK(err);
  err = dFSSetBlockSize(fsu,3);dCHK(err);
  err = dFSSetMesh(fsu,mesh,domain);dCHK(err);
  err = dFSSetDegree(fsu,jac,dutag);dCHK(err);
  for (dInt i=0; i<ALEN(vht->dirichlet) && vht->dirichlet[i]>0; i++) {
    err = dFSRegisterBoundary(fsu,vht->dirichlet[i],dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  }
  err = PetscObjectSetOptionsPrefix((dObject)fsu,"u");dCHK(err);
  err = dFSSetFromOptions(fsu);dCHK(err);
  err = PetscObjectSetName((PetscObject)fsu,"dFS_U_0");dCHK(err);
  vht->fsu = fsu;

  err = dFSCreate(vht->comm,&fsp);dCHK(err);
  err = dFSSetMesh(fsp,mesh,domain);dCHK(err);
  err = dFSSetDegree(fsp,jac,dptag);dCHK(err);
  err = PetscObjectSetOptionsPrefix((dObject)fsp,"p");dCHK(err);
  /* No boundaries, the pressure space has Neumann conditions when Dirichlet velocity conditions are applied */
  err = dFSSetFromOptions(fsp);dCHK(err);
  err = PetscObjectSetName((PetscObject)fsp,"dFS_P_0");dCHK(err);
  vht->fsp = fsp;

  err = dFSCreate(vht->comm,&fse);dCHK(err);
  err = dFSSetMesh(fse,mesh,domain);dCHK(err);
  err = dFSSetDegree(fse,jac,detag);dCHK(err);
  err = dFSRegisterBoundary(fse,100,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSRegisterBoundary(fse,200,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = dFSRegisterBoundary(fse,300,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  err = PetscObjectSetOptionsPrefix((dObject)fse,"e");dCHK(err);
  err = dFSSetFromOptions(fse);dCHK(err);
  err = PetscObjectSetName((PetscObject)fse,"dFS_E_0");dCHK(err);
  vht->fse = fse;

  err = dFSCreateExpandedVector(fsu,&vht->xu);dCHK(err);
  err = VecDuplicate(vht->xu,&vht->yu);dCHK(err);

  err = dFSCreateExpandedVector(fsp,&vht->xp);dCHK(err);
  err = VecDuplicate(vht->xp,&vht->yp);dCHK(err);

  err = dFSCreateExpandedVector(fsu,&vht->xe);dCHK(err);
  err = VecDuplicate(vht->xe,&vht->ye);dCHK(err);

  {
    dInt nu,np,ne,nul,npl,nel;
    err = dFSCreateGlobalVector(vht->fsu,&vht->gvelocity);dCHK(err);
    err = dFSCreateGlobalVector(vht->fsp,&vht->gpressure);dCHK(err);
    err = dFSCreateGlobalVector(vht->fse,&vht->genthalpy);dCHK(err);
    err = PetscObjectSetName((PetscObject)vht->gvelocity,"Velocity");dCHK(err);
    err = PetscObjectSetName((PetscObject)vht->gpressure,"Pressure");dCHK(err);
    err = PetscObjectSetName((PetscObject)vht->genthalpy,"Enthalpy");dCHK(err);
    err = VecGetLocalSize(vht->gvelocity,&nu);dCHK(err);
    err = VecGetLocalSize(vht->gpressure,&np);dCHK(err);
    err = VecGetLocalSize(vht->genthalpy,&ne);dCHK(err);

    {                           /* Get local sizes of the closure */
      Vec  Vc,Vgh,Pc,Pgh,Ec,Egh;
      err = VecDohpGetClosure(vht->gvelocity,&Vc);dCHK(err);
      err = VecDohpGetClosure(vht->gpressure,&Pc);dCHK(err);
      err = VecDohpGetClosure(vht->genthalpy,&Ec);dCHK(err);
      err = VecGhostGetLocalForm(Vc,&Vgh);dCHK(err);
      err = VecGhostGetLocalForm(Pc,&Pgh);dCHK(err);
      err = VecGhostGetLocalForm(Ec,&Egh);dCHK(err);
      err = VecGetLocalSize(Vgh,&nul);dCHK(err);
      err = VecGetLocalSize(Pgh,&npl);dCHK(err);
      err = VecGetLocalSize(Egh,&nel);dCHK(err);
      err = VecGhostRestoreLocalForm(Vc,&Vgh);dCHK(err);
      err = VecGhostRestoreLocalForm(Pc,&Pgh);dCHK(err);
      err = VecGhostRestoreLocalForm(Ec,&Egh);dCHK(err);
      err = VecDohpRestoreClosure(vht->gvelocity,&Vc);dCHK(err);
      err = VecDohpRestoreClosure(vht->gpressure,&Pc);dCHK(err);
      err = VecDohpRestoreClosure(vht->genthalpy,&Ec);dCHK(err);
    }

    {                           /* Set up the Stokes sub-problem */
      IS   ublock,pblock;
      dInt rstart;
      err = VecCreateMPI(vht->comm,nu+np,PETSC_DETERMINE,&vht->stokes.x);dCHK(err);
      err = VecDuplicate(vht->stokes.x,&vht->stokes.y);dCHK(err);
      err = VecGetOwnershipRange(vht->stokes.x,&rstart,NULL);dCHK(err);
      err = ISCreateStride(vht->comm,nu,rstart,1,&ublock);dCHK(err);
      err = ISCreateStride(vht->comm,np,rstart+nu,1,&pblock);dCHK(err);
      err = ISSetBlockSize(ublock,3);dCHK(err);
      err = VecScatterCreate(vht->stokes.x,ublock,vht->gvelocity,NULL,&vht->stokes.extractVelocity);dCHK(err);
      err = VecScatterCreate(vht->stokes.x,pblock,vht->gpressure,NULL,&vht->stokes.extractPressure);dCHK(err);
      vht->stokes.ublock = ublock;
      vht->stokes.pblock = pblock;
      /* Create local index sets */
      err = ISCreateStride(PETSC_COMM_SELF,nul,0,1,&vht->stokes.lublock);dCHK(err);
      err = ISCreateStride(PETSC_COMM_SELF,npl,nul,1,&vht->stokes.lpblock);dCHK(err);
      err = ISSetBlockSize(vht->stokes.lublock,3);dCHK(err);
    }
    {                           /* Set up the Stokes sub-problem */
      IS   ublock,pblock,eblock;
      dInt rstart;
      err = VecCreateMPI(vht->comm,nu+np+ne,PETSC_DETERMINE,&vht->gpacked);dCHK(err);
      err = VecGetOwnershipRange(vht->gpacked,&rstart,NULL);dCHK(err);
      err = ISCreateStride(vht->comm,nu,rstart,1,&ublock);dCHK(err);
      err = ISCreateStride(vht->comm,np,rstart+nu,1,&pblock);dCHK(err);
      err = ISCreateStride(vht->comm,ne,rstart+nu+np,1,&eblock);dCHK(err);
      err = ISSetBlockSize(ublock,3);dCHK(err);
      err = VecScatterCreate(vht->gpacked,ublock,vht->gvelocity,NULL,&vht->all.extractVelocity);dCHK(err);
      err = VecScatterCreate(vht->gpacked,pblock,vht->gpressure,NULL,&vht->all.extractPressure);dCHK(err);
      err = VecScatterCreate(vht->gpacked,eblock,vht->genthalpy,NULL,&vht->all.extractEnthalpy);dCHK(err);
      vht->all.ublock = ublock;
      vht->all.pblock = pblock;
      vht->all.eblock = eblock;
      /* Create local index sets */
      err = ISCreateStride(PETSC_COMM_SELF,nul,0,1,&vht->all.lublock);dCHK(err);
      err = ISCreateStride(PETSC_COMM_SELF,npl,nul,1,&vht->all.lpblock);dCHK(err);
      err = ISCreateStride(PETSC_COMM_SELF,nel,nul+npl,1,&vht->all.leblock);dCHK(err);
      err = ISSetBlockSize(vht->all.lublock,3);dCHK(err);
    }
  }
  err = dJacobiDestroy(&jac);dCHK(err);
  err = dMeshDestroy(&mesh);dCHK(err);

  err = VHTCaseSetType(vht->scase,scasename);dCHK(err);
  err = dFSGetBoundingBox(vht->fsu,vht->scase->bbox);dCHK(err);
  err = VHTCaseSetFromOptions(vht->scase);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTGetRegionIterator(VHT vht,VHTEvaluation eval,dRulesetIterator *riter)
{
  dErr err;

  dFunctionBegin;
  if (!vht->regioniter[eval]) {
    dRulesetIterator iter;
    dRuleset ruleset;
    dFS cfs;
    dMeshESH domain;
    dQuadratureMethod qmethod;
    switch (eval) {
    case EVAL_FUNCTION: qmethod = vht->function_qmethod; break;
    case EVAL_JACOBIAN: qmethod = vht->jacobian_qmethod; break;
    default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Unknown evaluation context");
    }
    err = dFSGetDomain(vht->fsu,&domain);dCHK(err);
    err = dFSGetPreferredQuadratureRuleSet(vht->fsu,domain,dTYPE_REGION,dTOPO_ALL,qmethod,&ruleset);dCHK(err);
    err = dFSGetCoordinateFS(vht->fsu,&cfs);dCHK(err);
    err = dRulesetCreateIterator(ruleset,cfs,&iter);dCHK(err);
    err = dRulesetDestroy(&ruleset);dCHK(err); /* Give ownership to iterator */
    err = dRulesetIteratorAddFS(iter,vht->fsu);dCHK(err);
    err = dRulesetIteratorAddFS(iter,vht->fsp);dCHK(err);
    err = dRulesetIteratorAddFS(iter,vht->fse);dCHK(err);
    if (eval == EVAL_FUNCTION) {err = dRulesetIteratorAddStash(iter,0,sizeof(struct VHTStash));dCHK(err);}
    vht->regioniter[eval] = iter;
  }
  *riter = vht->regioniter[eval];
  dFunctionReturn(0);
}

static dErr VHTExtractGlobalSplit(VHT vht,Vec X,Vec *Xu,Vec *Xp,Vec *Xe)
{
  dErr err;

  dFunctionBegin;
  if (Xu) {
    *Xu = vht->gvelocity;
    err = VecScatterBegin(vht->all.extractVelocity,X,*Xu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (vht->all.extractVelocity,X,*Xu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  if (Xp) {
    *Xp = vht->gpressure;
    err = VecScatterBegin(vht->all.extractPressure,X,*Xp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (vht->all.extractPressure,X,*Xp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  if (Xe) {
    *Xe = vht->genthalpy;
    err = VecScatterBegin(vht->all.extractEnthalpy,X,*Xe,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (vht->all.extractEnthalpy,X,*Xe,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr VHTCommitGlobalSplit(VHT vht,Vec *gxu,Vec *gxp,Vec *gxe,Vec gy,InsertMode imode)
{
  dErr err;

  dFunctionBegin;
  dASSERT(*gxu == vht->gvelocity);
  dASSERT(*gxp == vht->gpressure);
  dASSERT(*gxe == vht->genthalpy);
  err = VecScatterBegin(vht->all.extractVelocity,*gxu,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractVelocity,*gxu,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(vht->all.extractPressure,*gxp,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractPressure,*gxp,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(vht->all.extractEnthalpy,*gxe,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractEnthalpy,*gxe,gy,imode,SCATTER_REVERSE);dCHK(err);
  *gxu = NULL;
  *gxp = NULL;
  dFunctionReturn(0);
}

static dErr VHTDestroy(VHT *invht)
{
  VHT vht = *invht;
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(&vht->fsu);dCHK(err);
  err = dFSDestroy(&vht->fsp);dCHK(err);
  err = dFSDestroy(&vht->fse);dCHK(err);
  err = VecDestroy(&vht->xu);dCHK(err);
  err = VecDestroy(&vht->yu);dCHK(err);
  err = VecDestroy(&vht->xp);dCHK(err);
  err = VecDestroy(&vht->yp);dCHK(err);
  err = VecDestroy(&vht->xe);dCHK(err);
  err = VecDestroy(&vht->ye);dCHK(err);
  err = VecDestroy(&vht->gvelocity);dCHK(err);
  err = VecDestroy(&vht->gpressure);dCHK(err);
  err = VecDestroy(&vht->genthalpy);dCHK(err);
  err = VecDestroy(&vht->gpacked);dCHK(err);
  {
    err = ISDestroy(&vht->stokes.ublock);dCHK(err);
    err = ISDestroy(&vht->stokes.pblock);dCHK(err);
    err = ISDestroy(&vht->stokes.lublock);dCHK(err);
    err = ISDestroy(&vht->stokes.lpblock);dCHK(err);
    err = VecScatterDestroy(&vht->stokes.extractVelocity);dCHK(err);
    err = VecScatterDestroy(&vht->stokes.extractPressure);dCHK(err);
    err = VecDestroy(&vht->stokes.x);dCHK(err);
    err = VecDestroy(&vht->stokes.y);dCHK(err);
  }
  {
    err = ISDestroy(&vht->all.ublock);dCHK(err);
    err = ISDestroy(&vht->all.pblock);dCHK(err);
    err = ISDestroy(&vht->all.eblock);dCHK(err);
    err = ISDestroy(&vht->all.lublock);dCHK(err);
    err = ISDestroy(&vht->all.lpblock);dCHK(err);
    err = ISDestroy(&vht->all.leblock);dCHK(err);
    err = VecScatterDestroy(&vht->all.extractVelocity);dCHK(err);
    err = VecScatterDestroy(&vht->all.extractPressure);dCHK(err);
    err = VecScatterDestroy(&vht->all.extractEnthalpy);dCHK(err);
    err = VecScatterDestroy(&vht->all.extractStokes);dCHK(err);
  }
  for (dInt i=0; i<EVAL_UB; i++) {err = dRulesetIteratorDestroy(&vht->regioniter[i]);dCHK(err);}
  err = VHTCaseDestroy(&vht->scase);dCHK(err);
  err = dFree(*invht);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTGetMatrices(VHT vht,dBool use_jblock,Mat *J,Mat *P)
{
  dErr err;
  dInt m,nu,np,ne;
  Mat Juu,Jup,Jue,Jpu,Jpp,Jpe,Jeu,Jep,Jee,Buu,Bpp,Bee;
  IS splitis[3];

  dFunctionBegin;
  err = VecGetLocalSize(vht->gpacked,&m);dCHK(err);
  err = VecGetLocalSize(vht->gvelocity,&nu);dCHK(err);
  err = VecGetLocalSize(vht->gpressure,&np);dCHK(err);
  err = VecGetLocalSize(vht->genthalpy,&ne);dCHK(err);

  /* Create high-order matrix for diagonal velocity block, with context \a vht */
  err = MatCreateShell(vht->comm,nu,nu,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Juu);dCHK(err);
  err = MatShellSetOperation(Juu,MATOP_GET_VECS,(void(*)(void))MatGetVecs_VHT_stokes);dCHK(err);
  err = MatShellSetOperation(Juu,MATOP_MULT,(void(*)(void))MatMult_VHT_uu);dCHK(err);
  err = MatShellSetOperation(Juu,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_VHT_uu);dCHK(err);
  err = MatShellSetOperation(Juu,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_VHT_uu);dCHK(err);
  err = MatShellSetOperation(Juu,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))MatMultAdd_VHT_uu);dCHK(err);
  err = MatSetOptionsPrefix(Juu,"Juu_");dCHK(err);

  /* Create off-diagonal high-order matrix, with context \a vht */
  err = MatCreateShell(vht->comm,np,nu,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Jpu);dCHK(err);
  err = MatShellSetOperation(Jpu,MATOP_GET_VECS,(void(*)(void))MatGetVecs_VHT_stokes);dCHK(err);
  err = MatShellSetOperation(Jpu,MATOP_MULT,(void(*)(void))MatMult_VHT_pu);dCHK(err);
  err = MatShellSetOperation(Jpu,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_VHT_up);dCHK(err);
  err = MatShellSetOperation(Jpu,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_VHT_pu);dCHK(err);
  err = MatShellSetOperation(Jpu,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))MatMultAdd_VHT_up);dCHK(err);
  err = MatCreateTranspose(Jpu,&Jup);dCHK(err);
  err = MatSetOptionsPrefix(Jpu,"Jpu_");dCHK(err);
  err = MatSetOptionsPrefix(Jup,"Jup_");dCHK(err);

  /* These entries are really zero */
  Jpp = NULL;
  Jpe = NULL;
  Jep = NULL;

  /* @todo These off-diagonal blocks are not actually zero. Assume coupled application of the Jacobian and additive fieldsplit at this point */
  Jue = NULL;
  Jeu = NULL;

  /* Enthalpy-enthalpy coupling */
  err = MatCreateShell(vht->comm,ne,ne,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Jee);dCHK(err);
  err = MatShellSetOperation(Jee,MATOP_GET_VECS,(void(*)(void))MatGetVecs_VHT_ee);dCHK(err);
  err = MatShellSetOperation(Jee,MATOP_MULT,(void(*)(void))MatMult_VHT_ee);dCHK(err);
  err = MatSetOptionsPrefix(Jee,"Jee_");dCHK(err);

  splitis[0] = vht->all.ublock;
  splitis[1] = vht->all.pblock;
  splitis[2] = vht->all.eblock;
  /* Create the matrix-free operator */
  err = MatCreateNest(vht->comm,3,splitis,3,splitis,((Mat[]){Juu,Jup,Jue, Jpu,Jpp,Jpe, Jeu,Jep,Jee}),J);dCHK(err);
  err = MatSetOptionsPrefix(*J,"J_");dCHK(err);
  err = MatSetFromOptions(*J);dCHK(err);
  if (!use_jblock) {
    err = MatShellSetOperation(*J,MATOP_MULT,(void(*)(void))MatMult_Nest_VHT_all);dCHK(err);
  }

  err = MatDestroy(&Juu);dCHK(err);
  err = MatDestroy(&Jup);dCHK(err);
  err = MatDestroy(&Jue);dCHK(err);
  err = MatDestroy(&Jpu);dCHK(err);
  err = MatDestroy(&Jpp);dCHK(err);
  err = MatDestroy(&Jpe);dCHK(err);
  err = MatDestroy(&Jeu);dCHK(err);
  err = MatDestroy(&Jep);dCHK(err);
  err = MatDestroy(&Jee);dCHK(err);

  /* Create real matrix to be used for preconditioning */
  err = dFSGetMatrix(vht->fsu,vht->mattype_Buu,&Buu);dCHK(err);
  err = dFSGetMatrix(vht->fsp,vht->mattype_Bpp,&Bpp);dCHK(err);
  err = dFSGetMatrix(vht->fse,vht->mattype_Bee,&Bee);dCHK(err);
  err = MatSetOptionsPrefix(Buu,"Buu_");dCHK(err);
  err = MatSetOptionsPrefix(Bpp,"Bpp_");dCHK(err);
  err = MatSetOptionsPrefix(Bee,"Bee_");dCHK(err);
  err = MatSetOption(Buu,MAT_SYMMETRIC,PETSC_TRUE);dCHK(err);
  err = MatSetOption(Bpp,MAT_SYMMETRIC,PETSC_TRUE);dCHK(err);
  err = MatSetFromOptions(Buu);dCHK(err);
  err = MatSetFromOptions(Bpp);dCHK(err);
  err = MatSetFromOptions(Bee);dCHK(err);
  err = MatCreateNest(vht->comm,3,splitis,3,splitis,((Mat[]){Buu,NULL,NULL, NULL,Bpp,NULL, NULL,NULL,Bee}),P);dCHK(err);
  err = MatSetOptionsPrefix(*P,"B_");dCHK(err);
  err = MatSetFromOptions(*P);dCHK(err);

  err = MatDestroy(&Buu);dCHK(err);
  err = MatDestroy(&Bpp);dCHK(err);
  err = MatDestroy(&Bee);dCHK(err);
  dFunctionReturn(0);
}

// The "physics" functions below propagate derivatives
static inline void VHTRheoTransition(dScalar a,dScalar b,dReal width,dScalar x,dScalar *y,dScalar *dy)
{ // Smooth transition from state a to state b over width
  *y = a + (b-a)*0.5*(1 + tanh(x/width));
  *dy = (b-a)/width * 0.5 * (1 - dSqr(tanh(x/width)));
}
static inline void VHTRheoKappa(struct VHTRheology *rheo,dScalar e,dScalar *kappa,dScalar *kappa_e)
{
  VHTRheoTransition(rheo->kappa0,rheo->kappa1,1.0,e,kappa,kappa_e);
}
static inline void VHTTemperature(struct VHTRheology *rheo,dScalar e,dScalar *T,dScalar *T_e)
{
  *T = rheo->T0 + e;
  *T_e = 1;
}
static inline void VHTArrhenius(struct VHTRheology *rheo,dScalar e,dScalar *B,dScalar *B_e)
{
  dScalar n = 1./(rheo->pe-1),T,T_e,exparg,exparg_e;
  VHTTemperature(rheo,e,&T,&T_e);
  exparg   = rheo->Q / (rheo->R * T * n);
  exparg_e = -exparg / T * T_e;
  *B       = rheo->B0 * exp(exparg);
  *B_e     = rheo->B0 * exp(exparg) * exparg_e;
}
static inline void VHTRheoViscosity(struct VHTRheology *rheo,const dScalar Du[6],dScalar e,dScalar *eta,dScalar *eta_gamma,dScalar *eta_e)
{
  const dScalar
    p = rheo->pe,
    gamma_reg = 0.5*dSqr(rheo->eps) + 0.5*dColonSymScalar3(Du,Du),
    power = pow(gamma_reg,0.5*(p-2)),
    power_gamma = 0.5*(p-2) * power / gamma_reg;
  dScalar B,B_e;
  VHTArrhenius(rheo,e,&B,&B_e);
  *eta = B * power;
  *eta_gamma = B * power_gamma;
  *eta_e = B_e * power;
}

static inline void VHTPointwiseComputeStash(struct VHTRheology *rheo,const dReal dUNUSED x[3],const dScalar u[3],const dScalar Du[6],const dScalar e[1],const dScalar de[3],struct VHTStash *st)
{
  VHTRheoViscosity(rheo,Du,e[0],&st->eta,&st->eta_gamma,&st->eta_e);
  for (dInt i=0; i<3; i++) st->u[i] = u[i];
  for (dInt i=0; i<6; i++) st->Du[i] = Du[i];
  VHTRheoKappa(rheo,e[0],&st->kappa,&st->kappa_e);
  st->e = e[0];
  for (dInt i=0; i<3; i++) st->de[i] = de[i];
}

static inline void VHTPointwiseFunction(VHTCase scase,const dReal x[3],dReal weight,
                                        const dScalar u[3],const dScalar Du[6],const dScalar p[1],const dScalar e[1],const dScalar de[3],
                                        struct VHTStash *st,
                                        dScalar u_[3],dScalar Du_[6],dScalar p_[1],dScalar e_[1],dScalar de_[3])
{
  dScalar fu[3],fp[1],fe[1],Sigma;
  VHTPointwiseComputeStash(&scase->rheo,x,u,Du,e,de,st);
  scase->forcing(scase,x,fu,fp,fe);
  for (dInt i=0; i<3; i++) u_[i] = -weight * fu[i];                            // Momentum forcing term
  for (dInt i=0; i<6; i++) Du_[i] = weight * (st->eta * Du[i] - (i<3)*p[0]);   // eta Dv:Du - p tr(Dv)
  p_[0] = -weight * (Du[0]+Du[1]+Du[2] + fp[0]);                               // -q tr(Du) - forcing, note tr(Du) = div(u)
  Sigma = st->eta * dColonSymScalar3(Du,Du);                                   // Strain heating (Sigma)
  e_[0] = -weight * (Sigma + fe[0]);                                           // Strain heating and thermal forcing
  for (dInt i=0; i<3; i++) de_[i] = weight * (-u[i]*e[0] + st->kappa * de[i]); // Transport and diffusion
}

static inline void VHTPointwiseJacobian(const struct VHTStash *restrict st,dReal weight,
                                        const dScalar u[3],const dScalar Du[6],const dScalar p[1],const dScalar e[1],const dScalar de[3],
                                        dScalar u_[3],dScalar Du_[6],dScalar p_[1],dScalar e_[1],dScalar de_[3])
{
  const dScalar deta_colon = st->eta_gamma*dColonSymScalar3(st->Du,Du);                     // eta' Du:Dw
  for (dInt i=0; i<3; i++) u_[i] = 0;
  for (dInt i=0; i<6; i++) Du_[i] = weight * (+ st->eta * Du[i]              // eta Dv:Dw
                                              + deta_colon * st->Du[i]       // Dv : [eta' Du (x) Du] : Dw
                                              + st->eta_e * e[1] * st->Du[i] // eta_e e Du
                                              - p[0]*(i<3));                 // tr(Dv) p
  p_[0] = -weight*(Du[0]+Du[1]+Du[2]);                                       // -q tr(Du)
  e_[0] = weight * (- 2 * st->eta * dColonSymScalar3(st->Du,Du)              // Differentiate Sigma = eta(Du,e) Du : Du
                    - deta_colon * dColonSymScalar3(st->Du,st->Du)
                    - st->eta_e * e[0] * dColonSymScalar3(st->Du,st->Du));
  for (dInt i=0; i<3; i++) de_[i] = weight * (- st->u[i] *     e[0]
                                              -     u[i] * st->e
                                              + st->kappa * de[i]
                                              + st->kappa_e * e[0] * st->de[i]);
}

static inline void VHTPointwiseJacobian_uu(const struct VHTStash *st,dReal weight,const dScalar Du[6],dScalar Dv[6])
{
  const dScalar deta_colon = st->eta_gamma*dColonSymScalar3(st->Du,Du);
  for (dInt i=0; i<6; i++) Dv[i] = weight * (st->eta*Du[i] + deta_colon*st->Du[i]);
}

static inline void VHTPointwiseJacobian_pu(dReal weight,const dScalar Du[6],dScalar q[1])
{
  q[0] = -weight*(Du[0]+Du[1]+Du[2]);
}

static inline void VHTPointwiseJacobian_up(dReal weight,dScalar p,dScalar Dv[6])
{
  for (dInt i=0; i<6; i++) Dv[i] = -weight*p*(i<3);
}

static inline void VHTPointwiseJacobian_ee(const struct VHTStash *st,dReal weight,const dScalar e[1],const dScalar de[3],dScalar e_[1],dScalar de_[3])
{
  e_[0] = weight * st->eta_e * e[0] * dColonSymScalar3(st->Du,st->Du);
  for (dInt i=0; i<3; i++) de_[i] = weight * (- st->u[i] * e[0]
                                              + st->kappa * de[i]
                                              + st->kappa_e * e[0] * st->de[i]);
}

static dErr VHTFunction(SNES dUNUSED snes,Vec X,Vec Y,void *ctx)
{
  VHT              vht = ctx;
  dErr             err;
  Vec              Coords,Xu,Xp,Xe;
  dRulesetIterator iter;

  dFunctionBegin;
  err = VHTExtractGlobalSplit(vht,X,&Xu,&Xp,&Xe);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,Xu,dFS_INHOMOGENEOUS, Xp,dFS_INHOMOGENEOUS,Xp,dFS_INHOMOGENEOUS, Xe,dFS_INHOMOGENEOUS,Xe,dFS_INHOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*p)[1],(*e)[1],(*de)[3];
    dScalar (*u_)[3],(*du_)[9],(*p_)[1],(*e_)[1],(*de_)[3];
    dInt Q;
    struct VHTStash *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,&u_,&du_, &p,NULL,&p_,NULL, &e,&de,&e_,&de_);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar Du[6],Dv[6];
      dTensorSymCompress3(du[i],Du);
      VHTPointwiseFunction(vht->scase,x[i],jw[i], u[i],Du,p[i],e[i],de[i], &stash[i], u_[i],Dv,p_[i],e_[i],de_[i]);
      dTensorSymUncompress3(Dv,du_[i]);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, u_,du_, p_,NULL, e_,de_);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = VHTCommitGlobalSplit(vht,&Xu,&Xp,&Xe,Y,INSERT_VALUES);dCHK(err);
  dFunctionReturn(0);
}

static dErr MatMult_Nest_VHT_all(Mat J,Vec X,Vec Y)
{
  VHT              vht;
  Vec              Coords,Xu,Xp,Xe;
  dRulesetIterator iter;
  dErr             err;
  Mat              A;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_VHTShellMult,J,X,Y,0);dCHK(err);
  err = MatNestGetSubMat(J,0,0,&A);dCHK(err);
  err = MatShellGetContext(A,(void**)&vht);dCHK(err);
  err = VHTExtractGlobalSplit(vht,X,&Xu,&Xp,&Xe);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_HOMOGENEOUS,Xu,dFS_HOMOGENEOUS, Xp,dFS_HOMOGENEOUS,Xp,dFS_HOMOGENEOUS, Xe,dFS_HOMOGENEOUS,Xe,dFS_HOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*p)[1],(*e)[1],(*de)[3];
    dScalar (*u_)[3],(*du_)[9],(*p_)[1],(*e_)[1],(*de_)[3];
    dInt Q;
    struct VHTStash *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,&u_,&du_, &p,NULL,&p_,NULL, &e,&de,&e_,&de_);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar Du[6],Dv[6];
      dTensorSymCompress3(du[i],Du);
      VHTPointwiseJacobian(&stash[i],jw[i],u[i],Du,p[i],e[i],de[i],u_[i],Dv,p_[i],e_[i],de_[i]);
      dTensorSymUncompress3(Dv,du_[i]);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, u_,du_, p_,NULL, e_,de_);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = VHTCommitGlobalSplit(vht,&Xu,&Xp,&Xe,Y,INSERT_VALUES);dCHK(err);
  err = PetscLogEventEnd(LOG_VHTShellMult,J,X,Y,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr MatMultXIorA_VHT_stokes(Mat A,Vec X,Vec Y,Vec Z,InsertMode imode,VHTMultMode mmode)
{
  VHT           vht;
  dRulesetIterator iter;
  Vec              Coords;
  dErr             err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_VHTShellMult,A,X,Y,Z);dCHK(err);
  err = MatShellGetContext(A,(void**)&vht);dCHK(err);
  {  /* Check that we have correct sizes */
    dInt nu,np,nx,ny;
    err = VecGetSize(vht->gvelocity,&nu);dCHK(err);
    err = VecGetSize(vht->gpressure,&np);dCHK(err);
    err = VecGetSize(X,&nx);dCHK(err);
    err = VecGetSize(Y,&ny);dCHK(err);
    switch (mmode) {
    case VHT_MULT_UU: dASSERT(nx==nu && ny==nu); break;
    case VHT_MULT_UP: dASSERT(nx==np && ny==nu); break;
    case VHT_MULT_PU: dASSERT(nx==nu && ny==np); break;
    default: dERROR(PETSC_COMM_SELF,1,"Sizes do not match, unknown mult operation");
    }
  }

  switch (imode) {
  case INSERT_VALUES:
    if (Z) dERROR(vht->comm,PETSC_ERR_ARG_INCOMP,"Cannot use INSERT_VALUES and set gz");
    Z = Y;
    err = VecZeroEntries(Z);dCHK(err);
    break;
  case ADD_VALUES:
    if (Z != Y) {
      err = VecCopy(Y,Z);dCHK(err);
    }
    break;
  default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"unsupported imode");
  }

  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  switch (mmode) {
  case VHT_MULT_UU:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, X,dFS_HOMOGENEOUS,Z,dFS_HOMOGENEOUS, NULL,             NULL,              NULL,NULL);dCHK(err);
    break;
  case VHT_MULT_UP:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,             Z,dFS_HOMOGENEOUS, X,dFS_HOMOGENEOUS,NULL,              NULL,NULL);dCHK(err);
    break;
  case VHT_MULT_PU:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, X,dFS_HOMOGENEOUS,NULL,              NULL,             Z,dFS_HOMOGENEOUS, NULL,NULL);dCHK(err);
    break;
  default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmode");
  }
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*du)[9],(*dv)[9],*p,*q;
    dInt Q;
    struct VHTStash *stash;
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    switch (mmode) {
    case VHT_MULT_UU:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,&du,NULL,&dv, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {
        dScalar Du[6],Dv[6];
        dTensorSymCompress3(du[i],Du);
        VHTPointwiseJacobian_uu(&stash[i],jw[i],Du,Dv);
        dTensorSymUncompress3(Dv,dv[i]);
      }
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,dv, NULL,NULL, NULL,NULL);dCHK(err);
      break;
    case VHT_MULT_UP:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,NULL,NULL,&dv, &p,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {
        dScalar Dv[6];
        VHTPointwiseJacobian_up(jw[i],p[i],Dv);
        dTensorSymUncompress3(Dv,dv[i]);
      }
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,dv, NULL,NULL, NULL,NULL);dCHK(err);
      break;
    case VHT_MULT_PU:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,&du,NULL,NULL, NULL,NULL,&q,NULL, NULL,NULL,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {
        dScalar Du[6];
        dTensorSymCompress3(du[i],Du);
        VHTPointwiseJacobian_pu(jw[i],Du,&q[i]);
      }
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,NULL, q,NULL, NULL,NULL);dCHK(err);
      break;
    default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmode");
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = PetscLogEventEnd(LOG_VHTShellMult,A,X,Y,Z);dCHK(err);
  dFunctionReturn(0);
}

static dErr MatMult_VHT_ee(Mat A,Vec X,Vec Y)
{
  VHT              vht;
  dRulesetIterator iter;
  Vec              Coords;
  dErr             err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_VHTShellMult,A,X,Y,0);dCHK(err);
  err = MatShellGetContext(A,(void**)&vht);dCHK(err);
  err = VecZeroEntries(Y);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,NULL, NULL,NULL, X,Y);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*e)[1],(*de)[3],(*e_)[1],(*de_)[3];
    dInt Q;
    struct VHTStash *stash;
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL, &e,&de,&e_,&de_);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      VHTPointwiseJacobian_ee(&stash[i],jw[i],e[i],de[i],e_[i],de_[i]);dCHK(err);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,NULL, NULL,NULL, e_,de_);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = PetscLogEventEnd(LOG_VHTShellMult,A,X,Y,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTJacobianAssemble_Velocity(VHT vht,Mat Buu,Vec Mdiag,Vec X)
{
  dRulesetIterator iter;
  Vec Coords,Xu;
  dScalar *Kflat;
  dErr err;

  dFunctionBegin;
  err = VHTExtractGlobalSplit(vht,X,&Xu,NULL,NULL);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_JACOBIAN,&iter);dCHK(err);
  if (Mdiag) {
    err = VecZeroEntries(Mdiag);dCHK(err);
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,Mdiag,dFS_HOMOGENEOUS, NULL,NULL, NULL,NULL);dCHK(err);
  } else {
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,NULL, NULL,NULL, NULL,NULL);dCHK(err);
  }
  err = dRulesetIteratorGetMatrixSpaceSplit(iter, NULL,NULL,NULL,NULL, NULL,&Kflat,NULL,NULL, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw,*interp_flat,*deriv_flat;
    const dInt *rowcol;
    dScalar (*x)[3],(*dx)[3][3],(*u)[3],(*du)[9],(*v)[3],(*e)[1],(*de)[3];
    dInt Q,P;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,&v,NULL, NULL,NULL,NULL,NULL, &e,&de,NULL,NULL);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
    {                           /* Scope so that we can declare new VLA pointers for convenient assembly */
      const dReal (*interp)[P] = (const dReal(*)[P])interp_flat;
      const dReal (*deriv)[P][3] = (const dReal(*)[P][3])deriv_flat;
      dScalar (*K)[3][P][3] = (dScalar(*)[3][P][3])Kflat;
      err = PetscMemzero(K,P*3*P*3*sizeof(K[0][0][0][0]));dCHK(err);
      for (dInt q=0; q<Q; q++) {
        struct VHTStash stash;
        dScalar Dusym[6];
        dTensorSymCompress3(du[q],Dusym);
        VHTPointwiseComputeStash(&vht->scase->rheo,x[q],u[q],Dusym,e[q],de[q],&stash);
        for (dInt j=0; j<P; j++) { /* trial functions */
          for (dInt fj=0; fj<3; fj++) {
            dScalar duu[3][3] = {{0},{0},{0}},dv[3][3],Duusym[6],Dvsym[6];
            duu[fj][0] = deriv[q][j][0];
            duu[fj][1] = deriv[q][j][1];
            duu[fj][2] = deriv[q][j][2];
            dTensorSymCompress3(&duu[0][0],Duusym);
            VHTPointwiseJacobian_uu(&stash,jw[q],Duusym,Dvsym);
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
      err = dFSMatSetValuesBlockedExpanded(vht->fsu,Buu,8,rowcol,8,rowcol,&K[0][0][0][0],ADD_VALUES);dCHK(err);
      for (dInt i=0; i<P; i++) {
        dScalar Mentry = 0;
        for (dInt q=0; q<Q; q++) Mentry += interp[q][i] * jw[q] * interp[q][i]; /* Integrate the diagonal entry over this element */
        v[i][0] += Mentry;
        v[i][1] += Mentry;
        v[i][2] += Mentry;
      }
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, (dScalar**)&v,NULL, NULL,NULL, NULL,NULL);dCHK(err);
    err = dRulesetIteratorRestorePatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTJacobianAssemble_PressureEnthalpy(VHT vht,Mat Bpp,Mat Daux,Mat Bee,Vec X)
{
  dRulesetIterator iter;
  Vec              Coords,Xu,Xe;
  dScalar          *Kpp_flat,*Kppaux_flat,*Kee_flat;
  const dInt       *Ksizes;
  dErr             err;

  dFunctionBegin;
  /* It might seem weird to be getting velocity and enthalpy in the pressure assembly.  The reason is that this preconditioner
  * (indeed the entire problem) is always linear in pressure.  It \e might be nonlinear in velocity and enthalpy. */
  err = VHTExtractGlobalSplit(vht,X,&Xu,NULL,&Xe);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_JACOBIAN,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,NULL, NULL,NULL, Xe,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSplit(iter, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL, NULL,NULL,&Kpp_flat,NULL, NULL,NULL,NULL,&Kee_flat);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSizes(iter,NULL,NULL,&Ksizes);dCHK(err);
  err = dMallocA(Ksizes[2*4+2],&Kppaux_flat);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw,*interpp_flat,*derivp_flat,*interpe_flat,*derive_flat;
    const dInt *rowcolp,*rowcole;
    dScalar (*x)[3],(*dx)[3][3],(*u)[3],(*du)[9],(*e)[1],(*de)[3];
    dInt Q,Pp,Pe;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,NULL, NULL,NULL,NULL,NULL, &e,&de,NULL,NULL);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL, &Pp,&rowcolp,&interpp_flat,&derivp_flat, &Pe,&rowcole,&interpe_flat,&derive_flat);dCHK(err);
    {
      const dReal (*interpp)[Pp] = (const dReal(*)[Pp])interpp_flat;
      const dReal (*derivp)[Pp][3] = (const dReal(*)[Pp][3])derivp_flat;
      const dReal (*interpe)[Pe] = (const dReal(*)[Pe])interpe_flat;
      const dReal (*derive)[Pe][3] = (const dReal(*)[Pe][3])derive_flat;
      dScalar (*Kpp)[Pp] = (dScalar(*)[Pp])Kpp_flat,(*Kppaux)[Pp] = (dScalar(*)[Pp])Kppaux_flat;
      dScalar (*Kee)[Pe] = (dScalar(*)[Pe])Kee_flat;
      err = PetscMemzero(Kpp,Pp*Pp*sizeof(Kpp[0][0]));dCHK(err);
      err = PetscMemzero(Kppaux,Pp*Pp*sizeof(Kppaux[0][0]));dCHK(err);
      err = PetscMemzero(Kee,Pe*Pe*sizeof(Kee[0][0]));dCHK(err);
      for (dInt q=0; q<Q; q++) {
        struct VHTStash stash;
        dScalar Dusym[6];
        dTensorSymCompress3(du[q],Dusym);
        VHTPointwiseComputeStash(&vht->scase->rheo,x[q],u[q],Dusym,e[q],de[q],&stash);dCHK(err);
        /* Pressure-pressure Jacobians  */
        for (dInt j=0; j<Pp; j++) { /* trial functions */
          for (dInt i=0; i<Pp; i++) {
            /* Scaled mass matrx */
            Kpp[i][j] += interpp[q][i] * jw[q] * (1./stash.eta) * interpp[q][j];
            /* Neumann Laplacian */
            Kppaux[i][j] += (+ derivp[q][i][0] * jw[q] * derivp[q][j][0]
                             + derivp[q][i][1] * jw[q] * derivp[q][j][1]
                             + derivp[q][i][2] * jw[q] * derivp[q][j][2]);
          }
        }
        /* Enthalpy-enthalpy Jacobian */
        for (dInt j=0; j<Pe; j++) {
          const dScalar ez[1] = {interpe[q][j]},dez[3] = {derive[q][j][0],derive[q][j][1],derive[q][j][2]};
          dScalar e_[1],de_[3];
          VHTPointwiseJacobian_ee(&stash,jw[q],ez,dez,e_,de_);
          for (dInt i=0; i<Pe; i++) {
            Kee[i][j] += (interpe[q][i] * e_[0]
                          + derive[q][i][0] * de_[0]
                          + derive[q][i][1] * de_[1]
                          + derive[q][i][2] * de_[2]);
          }
        }
      }
      err = dFSMatSetValuesBlockedExpanded(vht->fsp,Bpp,Pp,rowcolp,Pp,rowcolp,&Kpp[0][0],ADD_VALUES);dCHK(err);
      if (Daux) {err = dFSMatSetValuesBlockedExpanded(vht->fsp,Daux,Pp,rowcolp,Pp,rowcolp,&Kppaux[0][0],ADD_VALUES);dCHK(err);}
      err = dFSMatSetValuesBlockedExpanded(vht->fse,Bee,Pe,rowcole,Pe,rowcole,&Kee[0][0],ADD_VALUES);dCHK(err);
    }
    err = dRulesetIteratorRestorePatchAssembly(iter, NULL,NULL,NULL,NULL, &Pp,&rowcolp,&interpp_flat,&derivp_flat, &Pe,&rowcole,&interpe_flat,&derive_flat);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = dFree(Kppaux_flat);dCHK(err);
  dFunctionReturn(0);
}


static dErr VHTJacobian(SNES dUNUSED snes,Vec X,Mat *J,Mat *B,MatStructure *structure,void *ctx)
{
  VHT vht = ctx;
  dErr err;
  Mat Buu,Bpp,Daux,Bee;
  Vec Mdiag;

  dFunctionBegin;
  err = MatGetLocalSubMatrix(*B,vht->all.lublock,vht->all.lublock,&Buu);dCHK(err);
  err = MatGetLocalSubMatrix(*B,vht->all.lpblock,vht->all.lpblock,&Bpp);dCHK(err);
  err = MatGetLocalSubMatrix(*B,vht->all.leblock,vht->all.leblock,&Bee);dCHK(err);
  err = PetscObjectQuery((dObject)Bpp,"LSC_M_diag",(dObject*)&Mdiag);dCHK(err);
  err = PetscObjectQuery((dObject)Bpp,"LSC_L",(dObject*)&Daux);dCHK(err);
  err = MatZeroEntries(*B);dCHK(err);
  if (Daux) {err = MatZeroEntries(Daux);dCHK(err);}
  err = VHTJacobianAssemble_Velocity(vht,Buu,Mdiag,X);dCHK(err);
  err = VHTJacobianAssemble_PressureEnthalpy(vht,Bpp,Daux,Bee,X);dCHK(err);
  if (Daux) {
    err = MatAssemblyBegin(Daux,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd  (Daux,MAT_FINAL_ASSEMBLY);dCHK(err);
  }
  err = MatRestoreLocalSubMatrix(*B,vht->all.lublock,vht->all.lublock,&Buu);dCHK(err);
  err = MatRestoreLocalSubMatrix(*B,vht->all.lpblock,vht->all.lpblock,&Bpp);dCHK(err);
  err = MatRestoreLocalSubMatrix(*B,vht->all.leblock,vht->all.leblock,&Bee);dCHK(err);

  /* MatNest calls assembly on the constituent pieces */
  err = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);dCHK(err);
  if (*J != *B) {
    err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  }
  *structure = SAME_NONZERO_PATTERN;
  dFunctionReturn(0);
}

static dErr VHTGetPressureShift(VHT vht,Vec Xp,dScalar *pressureshift)
{
  dErr             err;
  dScalar          volume = 0,shift = 0;
  dRulesetIterator iter;
  Vec              Coords;

  dFunctionBegin;
  *pressureshift = 0;
  if (!vht->alldirichlet) dFunctionReturn(0);
   // Do a volume integral of the exact solution to that we can remove the constant pressure mode
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,NULL, Xp,dFS_INHOMOGENEOUS,NULL, NULL,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw;
    const dScalar (*x)[3],*p;
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,NULL,NULL,NULL, NULL,NULL,NULL,NULL, &p,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[3],duu[9],pp[1],dpp[3],ee[1],dee[3];
      err = vht->scase->solution(vht->scase,x[i],uu,duu,pp,dpp,ee,dee);dCHK(err);
      volume += jw[i];
      shift  += (pp[0] - p[i]) * jw[i]; // The computed pressure sum is zero, but the continuous integral may not be
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  *pressureshift = shift/volume;
  dFunctionReturn(0);
}

static dErr VHTErrorNorms(VHT vht,Vec X,dReal N0u[3],dReal N1u[3],dReal N0p[3],dReal N1p[3],dReal N0e[3],dReal N1e[3])
{
  dErr             err;
  Vec              Coords,Xu,Xp,Xe;
  PetscScalar      pressureshift;
  dRulesetIterator iter;

  dFunctionBegin;
  err = dNormsStart(N0u,N1u);dCHK(err);
  err = dNormsStart(N0p,N1p);dCHK(err);
  err = dNormsStart(N0e,N1e);dCHK(err);
  err = VHTExtractGlobalSplit(vht,X,&Xu,&Xp,&Xe);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = VHTGetPressureShift(vht,Xp,&pressureshift);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,NULL, Xp,dFS_INHOMOGENEOUS,NULL, Xe,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw;
    const dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*p)[1],(*dp)[3],(*e)[1],(*de)[3];
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,NULL, &p,&dp,NULL,NULL, &e,&de,NULL,NULL);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[3],duu[9],pp[1],dpp[3],ee[1],dee[3];
      err = vht->scase->solution(vht->scase,x[i],uu,duu,pp,dpp,ee,dee);dCHK(err);
      pp[0] -= pressureshift;
      err = dNormsUpdate(N0u,N1u,jw[i],3,uu,u[i],duu,du[i]);dCHK(err);
      err = dNormsUpdate(N0p,N1p,jw[i],1,pp,p[i],dpp,dp[i]);dCHK(err);
      err = dNormsUpdate(N0e,N1e,jw[i],1,ee,e[i],dee,de[i]);dCHK(err);
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = dNormsFinish(N0u,N1u);dCHK(err);
  err = dNormsFinish(N0p,N1p);dCHK(err);
  err = dNormsFinish(N0e,N1e);dCHK(err);
  dFunctionReturn(0);
}

// This function cannot runs separately for each field because the nodal basis may be different for each field
static dErr VHTGetSolutionField_All(VHT vht,dFS fs,dInt fieldnumber,Vec *insoln)
{
  Vec      Sol,Xc,Cvecg,Cvec;
  dScalar *x;
  const dScalar *coords;
  dInt     n,bs;
  dErr     err;

  dFunctionBegin;
  *insoln = 0;
  err = dFSCreateGlobalVector(fs,&Sol);dCHK(err);
  err = VecDohpGetClosure(Sol,&Xc);dCHK(err);
  err = dFSGetNodalCoordinatesGlobal(fs,&Cvecg);dCHK(err);
  err = VecDohpGetClosure(Cvecg,&Cvec);dCHK(err);
  err = VecGetLocalSize(Xc,&n);dCHK(err);
  err = VecGetBlockSize(Xc,&bs);dCHK(err);
  {
    dInt nc;
    err = VecGetLocalSize(Cvec,&nc);dCHK(err);
    if (nc*bs != n*3) dERROR(PETSC_COMM_SELF,1,"Coordinate vector has inconsistent size");
  }
  err = VecGetArray(Xc,&x);dCHK(err);
  err = VecGetArrayRead(Cvec,&coords);dCHK(err);
  for (dInt i=0; i<n/bs; i++) {
    dScalar u_unused[3],p_unused[1],du_unused[3*3],dp_unused[3],e_unused[1],de_unused[3];
    switch (fieldnumber) {
    case 0:
      err = vht->scase->solution(vht->scase,&coords[3*i],&x[i*bs],du_unused,p_unused,dp_unused,e_unused,de_unused);dCHK(err);
      break;
    case 1:
      err = vht->scase->solution(vht->scase,&coords[3*i],u_unused,du_unused,&x[i*bs],dp_unused,e_unused,de_unused);dCHK(err);
      break;
    case 2:
      err = vht->scase->solution(vht->scase,&coords[3*i],u_unused,du_unused,p_unused,dp_unused,&x[i*bs],de_unused);dCHK(err);
      break;
    default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Requested field number %D",fieldnumber);
    }
  }
  err = VecRestoreArray(Xc,&x);dCHK(err);
  err = VecRestoreArrayRead(Cvec,&coords);dCHK(err);
  err = VecDohpRestoreClosure(Cvecg,&Cvec);dCHK(err);
  err = dFSInhomogeneousDirichletCommit(fs,Xc);dCHK(err);
  err = VecDohpRestoreClosure(Sol,&Xc);dCHK(err);
  *insoln = Sol;
  dFunctionReturn(0);
}

/** Creates a solution vector, commits the closure to each FS, returns packed solution vector */
static dErr VHTGetSolutionVector(VHT vht,Vec *insoln)
{
  dErr err;
  Vec Xu,Xp,Xe,spacked;

  dFunctionBegin;
  *insoln = 0;
  err = VHTGetSolutionField_All(vht,vht->fsu,0,&Xu);dCHK(err);
  err = VHTGetSolutionField_All(vht,vht->fsp,1,&Xp);dCHK(err);
  err = VHTGetSolutionField_All(vht,vht->fse,2,&Xe);dCHK(err);
  err = VecDuplicate(vht->gpacked,&spacked);dCHK(err);
  err = VecScatterBegin(vht->all.extractVelocity,Xu,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractVelocity,Xu,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(vht->all.extractPressure,Xp,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractPressure,Xp,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(vht->all.extractEnthalpy,Xe,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractEnthalpy,Xe,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecDestroy(&Xu);dCHK(err);
  err = VecDestroy(&Xp);dCHK(err);
  err = VecDestroy(&Xe);dCHK(err);
  *insoln = spacked;
  dFunctionReturn(0);
}

static dErr VHTGetNullSpace(VHT vht,MatNullSpace *matnull)
{
  dErr err;
  Vec r;

  dFunctionBegin;
  err = VecDuplicate(vht->gpacked,&r);dCHK(err);
  err = VecZeroEntries(r);dCHK(err);
  err = VecSet(vht->gpressure,1);dCHK(err);
  err = VecScatterBegin(vht->all.extractPressure,vht->gpressure,r,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractPressure,vht->gpressure,r,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecNormalize(r,PETSC_NULL);dCHK(err);
  err = MatNullSpaceCreate(vht->comm,dFALSE,1,&r,matnull);dCHK(err);
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
    err = SNESDefaultComputeJacobian(snes,U,&expmat_fd,&expmat_fd,&mstruct,NULL);dCHK(err);
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
  VHT vht;
  MPI_Comm comm;
  Mat J,B;
  Vec R,X,Xsoln = NULL;
  SNES snes;
  dBool check_error,check_null,compute_explicit,use_jblock,viewdhm;
  dErr err;

  err = dInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  err = PetscLogEventRegister("VHTShellMult",MAT_CLASSID,&LOG_VHTShellMult);dCHK(err);

  err = VHTCaseRegisterAll();dCHK(err);
  err = VHTCreate(comm,&vht);dCHK(err);
  err = VHTSetFromOptions(vht);dCHK(err);

  err = VecDuplicate(vht->gpacked,&R);dCHK(err);
  err = VecDuplicate(R,&X);dCHK(err);

  err = PetscOptionsBegin(vht->comm,NULL,"VHT solver options",__FILE__);dCHK(err); {
    check_error = vht->scase->reality ? dFALSE : dTRUE;
    err = PetscOptionsBool("-check_error","Compute errors","",check_error,&check_error,NULL);dCHK(err);
    err = PetscOptionsBool("-use_jblock","Use blocks to apply Jacobian instead of unified (more efficient) version","",use_jblock=dFALSE,&use_jblock,NULL);dCHK(err);
    err = PetscOptionsBool("-viewdhm","View the solution","",viewdhm=dFALSE,&viewdhm,NULL);dCHK(err);
    err = PetscOptionsBool("-check_null","Check that constant pressure really is in the null space","",check_null=dFALSE,&check_null,NULL);dCHK(err);
    if (check_null) {
      err = PetscOptionsBool("-compute_explicit","Compute explicit Jacobian (only very small sizes)","",compute_explicit=dFALSE,&compute_explicit,NULL);dCHK(err);
    }
  } err = PetscOptionsEnd();dCHK(err);
  err = VHTGetMatrices(vht,use_jblock,&J,&B);dCHK(err);
  err = SNESCreate(comm,&snes);dCHK(err);
  err = SNESSetFunction(snes,R,VHTFunction,vht);dCHK(err);
  err = SNESSetJacobian(snes,J,B,VHTJacobian,vht);dCHK(err);
  err = SNESSetFromOptions(snes);dCHK(err);
  {
    KSP    ksp;
    PC     pc;

    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPGetPC(ksp,&pc);dCHK(err);
    err = PCFieldSplitSetIS(pc,"u",vht->all.ublock);dCHK(err);
    err = PCFieldSplitSetIS(pc,"p",vht->all.pblock);dCHK(err);
    err = PCFieldSplitSetIS(pc,"e",vht->all.eblock);dCHK(err);
  }
  err = VHTGetSolutionVector(vht,&Xsoln);dCHK(err);
  if (!vht->scase->reality) {
    dReal nrm;
    MatStructure mstruct;
    Vec b;
    err = VecDuplicate(X,&b);dCHK(err);
    err = VecZeroEntries(X);dCHK(err);
    err = SNESComputeFunction(snes,X,b);dCHK(err); /* -f */
    err = SNESComputeFunction(snes,Xsoln,R);dCHK(err);
    err = VecNorm(R,NORM_2,&nrm);dCHK(err);
    err = dPrintf(comm,"Norm of discrete residual for exact solution %g\n",nrm);dCHK(err);
    err = SNESComputeJacobian(snes,Xsoln,&J,&B,&mstruct);dCHK(err);
    err = MatMult(J,Xsoln,R);dCHK(err);
    err = VecAXPY(R,1,b);dCHK(err); /* Jx - f */
    err = VecNorm(R,NORM_2,&nrm);dCHK(err);
    err = dPrintf(comm,"Norm of discrete linear residual at exact solution %g\n",nrm);dCHK(err);
    err = VecDestroy(&b);dCHK(err);
  }

  if (vht->alldirichlet) {                             /* Set null space */
    KSP ksp;
    MatNullSpace matnull;
    err = VHTGetNullSpace(vht,&matnull);dCHK(err);
    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPSetNullSpace(ksp,matnull);dCHK(err);
    if (Xsoln) {err = MatNullSpaceRemove(matnull,Xsoln,NULL);dCHK(err);}
    err = MatNullSpaceDestroy(&matnull);dCHK(err);
  }
  if (check_null) {
    err = CheckNullSpace(snes,R,compute_explicit);dCHK(err);
  }
  err = VecZeroEntries(R);dCHK(err);
  err = VecZeroEntries(X);dCHK(err);
  err = SNESSolve(snes,NULL,X);dCHK(err); /* ###  SOLVE  ### */
  if (vht->alldirichlet) {
    MatNullSpace matnull;
    KSP ksp;
    err = SNESGetKSP(snes,&ksp);dCHK(err);
    err = KSPGetNullSpace(ksp,&matnull);dCHK(err); /* does not reference */
    err = MatNullSpaceRemove(matnull,X,NULL);dCHK(err);
  }
  if (check_error) {
    dReal NAu[3],NIu[3],N0u[3],N1u[3],N0p[3],N1p[3],N0e[3],N1e[3];
    err = VHTErrorNorms(vht,X,N0u,N1u,N0p,N1p,N0e,N1e);dCHK(err);
    err = dNormsAlgebraicScaled(NAu,R);dCHK(err);
    err = VecWAXPY(R,-1,Xsoln,X);dCHK(err);
    err = dNormsAlgebraicScaled(NIu,R);dCHK(err);
    err = dPrintf(comm,"Algebraic residual        |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",NAu[0],NAu[1],NAu[2]);dCHK(err);
    err = dPrintf(comm,"Interpolation residual    |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",NIu[0],NIu[1],NIu[2]);dCHK(err);
    err = dPrintf(comm,"Integral velocity error 0 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N0u[0],N0u[1],N0u[2]);dCHK(err);
    err = dPrintf(comm,"Integral velocity error 1 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N1u[0],N1u[1],N1u[2]);dCHK(err);
    err = dPrintf(comm,"Integral pressure error 0 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N0p[0],N0p[1],N0p[2]);dCHK(err);
    err = dPrintf(comm,"Integral pressure error 1 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N1p[0],N1p[1],N1p[2]);dCHK(err);
    err = dPrintf(comm,"Integral enthalpy error 0 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N0e[0],N0e[1],N0e[2]);dCHK(err);
    err = dPrintf(comm,"Integral enthalpy error 1 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N1e[0],N1e[1],N1e[2]);dCHK(err);
  }
  if (viewdhm) {
    Vec Xu,Xp,Xe;
    dViewer view;
    err = PetscViewerCreate(comm,&view);dCHK(err);
    err = PetscViewerSetType(view,PETSCVIEWERDHM);dCHK(err);
    err = PetscViewerFileSetName(view,"vht.dhm");dCHK(err);
    err = PetscViewerFileSetMode(view,FILE_MODE_WRITE);dCHK(err);
    err = VHTExtractGlobalSplit(vht,X,&Xu,&Xp,&Xe);dCHK(err);
    err = dFSDirichletProject(vht->fsu,Xu,dFS_INHOMOGENEOUS);dCHK(err);
    err = dFSDirichletProject(vht->fsp,Xp,dFS_INHOMOGENEOUS);dCHK(err);
    err = dFSDirichletProject(vht->fse,Xe,dFS_INHOMOGENEOUS);dCHK(err);
    err = VecView(Xu,view);dCHK(err);
    err = VecView(Xp,view);dCHK(err);
    err = VecView(Xe,view);dCHK(err);
    err = PetscViewerDestroy(&view);dCHK(err);
  }

  err = VecDestroy(&R);dCHK(err);
  err = VecDestroy(&X);dCHK(err);
  err = VecDestroy(&Xsoln);dCHK(err);
  err = SNESDestroy(&snes);dCHK(err);
  if (J != B) {err = MatDestroy(&J);dCHK(err);}
  err = MatDestroy(&B);dCHK(err);
  err = VHTDestroy(&vht);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
