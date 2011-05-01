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


dErr StokesCaseRegisterAll_Jako(void)
{
  dErr err;

  dFunctionBegin;
  err = StokesCaseRegisterAll_Exact();dCHK(err);
  err = StokesCaseRegister("gravity",StokesCaseCreate_Gravity);dCHK(err);
  dFunctionReturn(0);
}
