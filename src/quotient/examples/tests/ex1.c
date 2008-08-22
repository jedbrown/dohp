static const char help[] = "Test the construction of DohpQuotient objects.\n";

#include "dohp.h"
#include "dohpjacobi.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char *argv[])
{
  DohpQuotient quot;
  DohpJacobi jac;
  DohpMesh mesh;
  MPI_Comm comm;
  PetscViewer viewer;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscInitialize(&argc,&argv,0,help);CHKERRQ(ierr);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  ierr = DohpMeshCreate(comm,&mesh);CHKERRQ(ierr);
  ierr = DohpMeshLoad(mesh,"dblock.h5m","");CHKERRQ(ierr);
  ierr = DohpQuotientCreate(mesh,0,0,&quot);CHKERRQ(ierr);
  ierr = DohpQuotientSetFromOptions(quot);CHKERRQ(ierr);
  ierr = DohpQuotientSetUp(quot);CHKERRQ(ierr);
  ierr = DohpQuotientView(quot,viewer);CHKERRQ(ierr);
  ierr = DohpQuotientDestroy(quot);CHKERRQ(ierr);

  ierr = DohpJacobiCreate(comm,&jac);CHKERRQ(ierr);
  ierr = DohpJacobiSetFromOptions(jac);CHKERRQ(ierr);
  ierr = DohpJacobiView(jac,viewer);CHKERRQ(ierr);
  ierr = DohpJacobiSetUp(jac);CHKERRQ(ierr);
  ierr = DohpJacobiDestroy(jac);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
