static const char help[] = "Tests the dJacobi object.";

#include "dohpjacobi.h"

int main(int argc,char *argv[])
{
  dJacobi jac;
  MPI_Comm comm;
  PetscViewer viewer;
  dErr err;

  dFunctionBegin;
  err = PetscInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;

  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetDegrees(jac,8,4);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);
  err = dJacobiView(jac,viewer);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);
  err = dJacobiDestroy(jac);dCHK(err);
  err = PetscFinalize();dCHK(err);
  dFunctionReturn(0);
}
