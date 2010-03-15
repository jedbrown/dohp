static const char help[] = "Tests modal bases for dJacobi.";

#include <dohpjacobi.h>
#include <dohpsys.h>
#include <dohp.h>

static dErr TestModalBases(dJacobi jac,PetscViewer viewer)
{
  dEntTopology topo[] = {dTOPO_HEX,dTOPO_HEX,dTOPO_HEX};
  dPolynomialOrder rdegree[] = {dPolynomialOrderCreate(6,0,0,0),
                                dPolynomialOrderCreate(6,0,0,0),
                                dPolynomialOrderCreate(6,0,0,0)};
  dPolynomialOrder bdegree[] = {dPolynomialOrderCreate(0,0,0,0),
                                dPolynomialOrderCreate(1,0,0,0),
                                dPolynomialOrderCreate(2,0,0,0)};
  dRule *rules;
  dEFS *efs;
  dQuadrature quad;
  dErr err;

  dFunctionBegin;
  err = dJacobiGetQuadrature(jac,dQUADRATURE_METHOD_FAST,&quad);dCHK(err);
  err = dQuadratureGetRules(quad,3,topo,rdegree,&rules);dCHK(err);
  err = dJacobiGetEFS(jac,3,topo,bdegree,rules,&efs);dCHK(err);

  for (dInt i=0; i<3; i++) {
    const dScalar primes[10] = {2,3,5,7,11,13,17,19,23,29};
    dScalar *modes,*values,*derivs;
    dReal *coord,*weight;
    dInt m,n;
    err = dRuleGetSize(rules[i],NULL,&n);dCHK(err);
    err = dEFSGetSizes(efs[i],NULL,NULL,&m);dCHK(err);
    err = dMallocA5(m,&modes,3*n,&coord,n,&weight,n,&values,n*3,&derivs);dCHK(err);
    err = dMemcpy(modes,primes,m*sizeof(*modes));dCHK(err);
    err = dRuleGetNodeWeight(rules[i],coord,weight);dCHK(err);
    err = dEFSApply(efs[i],NULL,1,modes,values,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Values of %d modes evaluated at quadrature points\n",m);dCHK(err);
    err = PetscScalarView(n,values,viewer);dCHK(err);
    for (dInt j=0; j<n; j++) values[j] *= weight[j];
    err = dEFSApply(efs[i],NULL,1,values,modes,dAPPLY_INTERP_TRANSPOSE,INSERT_VALUES);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Modal values after application of mass matrix\n");dCHK(err);
    err = PetscScalarView(m,modes,viewer);dCHK(err);
    err = dFree5(modes,coord,weight,values,derivs);dCHK(err);
  }

  err = dFree(rules);dCHK(err);
  err = dFree(efs);dCHK(err);
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  dJacobi jac;
  MPI_Comm comm;
  PetscViewer viewer;
  dErr err;

  err = dInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;

  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = TestModalBases(jac,viewer);dCHK(err);
  err = dJacobiDestroy(jac);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
