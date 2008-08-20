static const char help[] = "Test the construction of DohpQuotient objects.\n";

#include "dohp.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char *argv[])
{
  DohpQuotient quot;
  DohpMesh mesh;
  MPI_Comm comm = PETSC_COMM_WORLD;

  PetscFunctionBegin;
  ierr = PetscInitialize(argc,argv,0,help);CHKERRQ(ierr);
  ierr = DohpMeshCreate(comm,&mesh);CHKERRQ(ierr);
  ierr = DohpMeshLoad(mesh,"dblock.h5m","");CHKERRQ(ierr);
  ierr = DohpMeshView(mesh,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
