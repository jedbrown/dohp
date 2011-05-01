#include "stokesimpl.h"

// Trivial gravity model
static void StokesCaseSolution_Gravity(StokesCase dUNUSED scase,const dReal dUNUSED x[3],dScalar u[],dScalar du[],dScalar *p,dScalar dp[])
{                               /* Defines inhomogeneous Dirichlet boundary conditions */
  u[0] = u[1] = u[2] = 0;
  for (dInt i=0; i<9; i++) du[i] = 0;
  *p = 0;
  for (dInt i=0; i<3; i++) dp[i] = 0;
}
static void StokesCaseForcing_Gravity(StokesCase scase,const dReal dUNUSED x[3],dScalar fu[],dScalar *fp)
{
  fu[0] = 0;
  fu[1] = 0;
  fu[2] = scase->gravity;
  fp[0] = 0;
}
static dErr StokesCaseCreate_Gravity(StokesCase scase)
{
  dFunctionBegin;
  scase->reality = dTRUE;
  scase->solution = StokesCaseSolution_Gravity;
  scase->forcing  = StokesCaseForcing_Gravity;
  dFunctionReturn(0);
}

// A real implementation
typedef struct {
  int placeholder;
} StokesCase_Jako;

static dErr StokesCaseSolution_Jako(StokesCase scase,const dReal x[3],dScalar u[],dScalar du[],dScalar *p,dScalar dp[])
{                               /* Defines inhomogeneous Dirichlet boundary conditions */
  StokesCase_Jako *jako = scase->data;

  dFunctionBegin;
  if (x[0] > 580000. || jako->placeholder) {
    u[0] = -200;
    u[1] = -100;
    u[2] = 0;
  } else {
    u[0] = -400;
    u[1] = -300;
    u[2] = 0;
  }

  for (dInt i=0; i<9; i++) du[i] = 0;
  *p = 0;
  for (dInt i=0; i<3; i++) dp[i] = 0;
  dFunctionReturn(0);
}
static void StokesCaseSolution_Jako_Void(StokesCase scase,const dReal x[3],dScalar u[],dScalar du[],dScalar *p,dScalar dp[])
{
  dErr err;
  err = StokesCaseSolution_Jako(scase,x,u,du,p,dp);CHKERRV(err);
}
static dErr StokesCaseSetFromOptions_Jako(StokesCase scase)
{
  char fname[256] = "unknown";
  dBool flg;
  dErr err;
  dFunctionBegin;
  err = PetscOptionsHead("StokesCase_Jako options");dCHK(err); {
    err = PetscOptionsString("-jako_surface_velocity","File to read surface velocity from (assume same projection as model, e.g. UTM)","",fname,fname,sizeof(fname),&flg);dCHK(err);
    if (!flg) dERROR(scase->comm,PETSC_ERR_USER,"User must provide surface velocity file with -jako_surface_velocity FILENAME");
  } err = PetscOptionsTail();dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseDestroy_Jako(StokesCase scase)
{
  StokesCase_Jako *jako = scase->data;
  dErr err;

  dFunctionBegin;
  jako->placeholder = 0;
  err = dFree(scase->data);dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseCreate_Jako(StokesCase scase)
{
  StokesCase_Jako *jako;
  dErr err;

  dFunctionBegin;
  scase->reality = dTRUE;
  scase->solution = StokesCaseSolution_Jako_Void;
  scase->forcing  = StokesCaseForcing_Gravity;
  scase->setfromoptions = StokesCaseSetFromOptions_Jako;
  scase->destroy = StokesCaseDestroy_Jako;

  err = dNew(StokesCase_Jako,&jako);dCHK(err);
  scase->data = jako;
  dFunctionReturn(0);
}

dErr StokesCaseRegisterAll_Jako(void)
{
  dErr err;

  dFunctionBegin;
  err = StokesCaseRegisterAll_Exact();dCHK(err);
  err = StokesCaseRegister("gravity",StokesCaseCreate_Gravity);dCHK(err);
  err = StokesCaseRegister("jako",StokesCaseCreate_Jako);dCHK(err);
  dFunctionReturn(0);
}
