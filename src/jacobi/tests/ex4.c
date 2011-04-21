static const char help[] = "Tests dRuleGetPatches for dJacobi.";

#include <dohpjacobi.h>
#include <dohpsys.h>
#include <dohp.h>

static dErr TestPatches(dJacobi jac,PetscViewer viewer)
{
  dEntTopology topo[] = {dTOPO_HEX,dTOPO_HEX,dTOPO_HEX,dTOPO_HEX};
  dPolynomialOrder rdegree[4],bdegree[4];
  dRule *rules;
  dEFS *efs;
  dQuadrature quad;
  dQuadratureMethod method = dQUADRATURE_METHOD_FAST;
  dInt bp = 5;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"TestPatches Options",NULL);dCHK(err);
  err = PetscOptionsInt("-basis_degree","Degree of tensor-product basis polynomial",NULL,bp,&bp,NULL);dCHK(err);
  err = PetscOptionsEnum("-quadrature_method","Method to use for sample quadrature",NULL,dQuadratureMethods,(PetscEnum)method,(PetscEnum*)&method,NULL);dCHK(err);
  err = PetscOptionsEnd();dCHK(err);
  for (dInt i=0; i<4; i++) {
    rdegree[i] = dPolynomialOrderCreate(0,(i+1)*2-2,(i+1)*2-2,(i+1)*2-2);
    bdegree[i] = dPolynomialOrderCreate(0,bp,bp,bp);
  }
  err = dJacobiGetQuadrature(jac,method,&quad);dCHK(err);
  err = dQuadratureGetRules(quad,4,topo,rdegree,&rules);dCHK(err);
  err = dJacobiGetEFS(jac,4,topo,bdegree,rules,&efs);dCHK(err);

  for (dInt i=0; i<4; i++) {
    dInt npatches,patchsize,dim,tsize[3];
    const dInt *ind;
    const dReal *weight,*tcoord[3],*tweight[3];
    dRule rule;
    err = dEFSGetRule(efs[i],&rule);dCHK(err);
    err = dRuleGetTensorNodeWeight(rule,&dim,tsize,tcoord,tweight);dCHK(err);
    err = dRuleGetPatches(rule,&npatches,&patchsize,&ind,&weight);dCHK(err);
    if (npatches*patchsize != tsize[0]*tsize[1]*tsize[2]) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Inconsistent patch sizes");
    err = PetscViewerASCIIPrintf(viewer,"Element %D: npatches %D  patchsize %D\n",i,npatches,patchsize);dCHK(err);
    err = dIntTableView(npatches,patchsize,ind,viewer,"indices");dCHK(err);
    err = dRealTableView(npatches,patchsize,weight,viewer,"weights");dCHK(err);
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
  err = TestPatches(jac,viewer);dCHK(err);
  err = dJacobiDestroy(&jac);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
