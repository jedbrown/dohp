static const char help[] = "Test the construction of dQuotient objects.\n";

#include <dohp.h>
#include <dohpjacobi.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char *argv[])
{
  dQuotient quot;
  dJacobi jac;
  dMesh mesh;
  MPI_Comm comm;
  PetscViewer viewer;
  dErr err;

  dFunctionBegin;
  err = dInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = dMeshCreate(comm,&mesh);dCHK(err);
  err = dMeshLoad(mesh,"dblock.h5m","");dCHK(err);
  err = dQuotientCreate(mesh,0,0,&quot);dCHK(err);
  err = dQuotientSetFromOptions(quot);dCHK(err);
  err = dQuotientSetUp(quot);dCHK(err);
  err = dQuotientView(quot,viewer);dCHK(err);
  err = dQuotientDestroy(quot);dCHK(err);

  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiView(jac,viewer);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);
  err = dJacobiDestroy(jac);dCHK(err);
  err = dFinalize();dCHK(err);
  dFunctionReturn(0);
}
