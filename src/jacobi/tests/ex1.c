static const char help[] = "Tests the dJacobi object.";

#include "dohpjacobi.h"
#include "private/fsimpl.h"

dErr checkRulesAndEFS(dJacobi);
  
int main(int argc,char *argv[])
{
  dJacobi jac;
  dEFS efs;
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
  /* err = dJacobiView(jac,viewer);dCHK(err); */
  err = dJacobiSetUp(jac);dCHK(err);
  err = checkRulesAndEFS(jac);dCHK(err);  
  err = dJacobiDestroy(jac);dCHK(err);
  err = PetscFinalize();dCHK(err);
  dFunctionReturn(0);
}

dErr checkRulesAndEFS(dJacobi jac)
{
  const dInt N = 10;
  const dInt rsize[3] = {4,5,6};
  const dInt bsize[3] = {3,4,5};
  const dTopology topo = iMesh_QUADRILATERAL;
  PetscViewer viewer = PETSC_VIEWER_STDOUT_WORLD;
  MPI_Comm comm = ((PetscObject)jac)->comm;
  dRule *rule;
  dEFS *efs;
  dInt index;
  void **rbase,**ebase;
  dErr err;

  dFunctionBegin;
  err = dMalloc(N*sizeof(dRule),&rule);dCHK(err);

  index = 0; rbase = NULL;
  do {
    if (index) {
      err = dMalloc(index*sizeof(*rbase),&rbase);dCHK(err);
      index = 0;
    }
    for (dInt i=0; i<N; i++) {
      err = dJacobiGetRule(jac,topo,rsize,&rule[i],rbase,&index);dCHK(err);
    }
  } while (!rbase);
  
  for (dInt i=0; i<N; i++) {
    err = dPrintf(comm,"Rule for element %d\n",i);dCHK(err);
    err = dRuleView(&rule[i],viewer);dCHK(err);
  }
  dFunctionReturn(0);
}
