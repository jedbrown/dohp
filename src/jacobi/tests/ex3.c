static const char help[] = "Tests explicit bases for dJacobi.";

#include <dohpjacobi.h>
#include <dohpsys.h>
#include <dohp.h>

static dErr MatRealView(dInt m,dInt n,const dReal *rvalues,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<m; i++) {
    for (dInt j=0; j<n; j++) {
      err = PetscViewerASCIIPrintf(viewer,"%+10.6f%s",rvalues[i*n+j],j+1==n?"\n":" ");dCHK(err);
    }
  }
  dFunctionReturn(0);
}

static dErr GetEFS(dJacobi jac,dRule **rules,dEFS **efs)
{
  dErr err;
  const dJacobiType jtype;
  enum {TENSOR,MODAL} type;
  dQuadratureMethod qmethod = dQUADRATURE_METHOD_FAST;
  dEntTopology topo[] = {dTOPO_HEX,dTOPO_HEX,dTOPO_HEX,dTOPO_HEX};
  dPolynomialOrder rdegree[4],bdegree[4];
  dQuadrature quad;
  dInt rp = 6;

  dFunctionBegin;
  err = dJacobiGetType(jac,&jtype);dCHK(err);
  if (!strcmp(jtype,dJACOBI_TENSOR)) type = TENSOR;
  else if (!strcmp(jtype,dJACOBI_MODAL)) type = MODAL;
  else dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unexpected Jacobi type '%s'",jtype);

  err = PetscOptionsGetEnum(NULL,"-qmethod",dQuadratureMethods,(PetscEnum*)&qmethod,NULL);dCHK(err);
  err = PetscOptionsGetInt(NULL,"-rule_degree",&rp,NULL);dCHK(err);
  for (dInt i=0; i<4; i++) {
    switch (type) {
      case TENSOR:
        rdegree[i] = dPolynomialOrderCreate(0,rp,rp,rp);
        bdegree[i] = dPolynomialOrderCreate(0,i+1,i+1,i+1);
        break;
      case MODAL:
        rdegree[i] = dPolynomialOrderCreate(rp,0,0,0);
        bdegree[i] = dPolynomialOrderCreate(i,0,0,0);
    }
  }
  err = dJacobiGetQuadrature(jac,qmethod,&quad);dCHK(err);
  err = dQuadratureGetRules(quad,4,topo,rdegree,rules);dCHK(err);
  err = dJacobiGetEFS(jac,4,topo,bdegree,*rules,efs);dCHK(err);
  dFunctionReturn(0);
}

static dErr TestExplicitBases(dJacobi jac,PetscViewer viewer)
{
  dRule *rules;
  dEFS *efs;
  dBool just_view = dFALSE;
  dErr err;

  dFunctionBegin;
  err = GetEFS(jac,&rules,&efs);dCHK(err);
  err = PetscOptionsGetBool(NULL,"-just_view",&just_view,NULL);dCHK(err);

  for (dInt i=0; i<4; i++) {
    dScalar *values,*derivs;
    dReal *coord,*weight;
    const dReal *interp,*deriv;
    dInt m,n,Q,P;
    if (just_view) {
      err = dEFSView(efs[i],PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
      continue;
    }
    err = dRuleGetSize(rules[i],NULL,&n);dCHK(err);
    err = dEFSGetSizes(efs[i],NULL,NULL,&m);dCHK(err);
    err = dMallocA4(3*n,&coord,n,&weight,n,&values,n*3,&derivs);dCHK(err);
    err = dRuleGetNodeWeight(rules[i],coord,weight);dCHK(err);
    err = dEFSGetExplicit(efs[i],NULL,&Q,&P,&interp,&deriv);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Interpolation matrix of size (%d,%d)\n",Q,P);dCHK(err);
    err = MatRealView(Q,P,interp,viewer);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Gradient matrix of size (%d,%d*3)\n",Q,P);dCHK(err);
    err = MatRealView(Q,P*3,deriv,viewer);dCHK(err);
    {
      dReal mass[P][P];
      err = dMemzero(mass,sizeof(mass));dCHK(err);
      for (dInt j=0; j<Q; j++) {
        for (dInt k=0; k<P; k++) {
          for (dInt l=0; l<P; l++) {
            mass[k][l] += interp[j*P+k] * weight[j] * interp[j*P+l];
          }
        }
      }
      err = PetscViewerASCIIPrintf(viewer,"Mass matrix of size (%d,%d)\n",P,P);dCHK(err);
      err = MatRealView(P,P,&mass[0][0],viewer);dCHK(err);
    }
    err = dEFSRestoreExplicit(efs[i],NULL,&Q,&P,&interp,&deriv);dCHK(err);
    err = dFree4(coord,weight,values,derivs);dCHK(err);
  }

  err = dFree(rules);dCHK(err);
  err = dFree(efs);dCHK(err);
  dFunctionReturn(0);
}

static dErr TestExplicitSparseBases(dJacobi jac,PetscViewer viewer)
{
  dErr err;
  dRule *rules;
  dEFS *efs;

  dFunctionBegin;
  err = GetEFS(jac,&rules,&efs);dCHK(err);
  for (dInt i=0; i<3; i++) {
    dRule rule;
    dInt npatches,Q,P;
    const dInt *qidx;
    const dReal *jdet;
    dInt *bidx;
    dReal *interp,*deriv;
    err = dEFSGetRule(efs[i],&rule);dCHK(err);
    err = dEFSGetSizes(efs[i],NULL,NULL,&P);dCHK(err); /* Patch support may be smaller, P will be redefined below */
    err = dRuleGetPatches(rule,&npatches,&Q,&qidx,&jdet);dCHK(err);
    err = dMallocA3(npatches*P,&bidx,npatches*Q*P,&interp,npatches*Q*P*3,&deriv);dCHK(err);
    err = dEFSGetExplicitSparse(efs[i],npatches,Q,qidx,NULL,0,&P,bidx,interp,deriv);dCHK(err);
    err = dIntTableView(npatches,Q,qidx,viewer,"elem %D: qidx(patch,q)",i);dCHK(err);
    err = dIntTableView(npatches,P,bidx,viewer,"elem %D: bidx(patch,b)",i);dCHK(err);
    for (dInt j=0; j<npatches; j++) {
      err = dRealTableView(Q,P,interp+Q*P*j,viewer,"elem %D  patch %D  interp",i,j);dCHK(err);
      err = dRealTableView(Q,P*3,deriv+Q*P*3*j,viewer,"elem %D  patch %D  deriv",i,j);dCHK(err);
    }
    err = dFree3(bidx,interp,deriv);dCHK(err);
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
  dBool sparse = dFALSE;

  err = dInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;

  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = PetscOptionsGetBool(NULL,"-sparse",&sparse,NULL);dCHK(err);
  if (sparse) {
    err = TestExplicitSparseBases(jac,viewer);dCHK(err);
  } else {
    err = TestExplicitBases(jac,viewer);dCHK(err);
  }
  err = dJacobiDestroy(jac);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
