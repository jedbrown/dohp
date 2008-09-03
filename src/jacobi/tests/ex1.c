static const char help[] = "Tests the dJacobi object.";

#include "dohpjacobi.h"
#include "private/fsimpl.h"
#include "petscvec.h"
#include <stdlib.h>

#define JACOBI_VIEW 0
#define RULE_VIEW 0
#define EFS_VIEW 0

dInt productInt(dInt,const dInt[]);
dErr checkRulesAndEFS(dJacobi);
dErr createFS(MPI_Comm comm,dInt dim,dInt N,dEFS *efs,Vec *U);
dErr checkInterp(dInt N,dEFS *efs,Vec u);
dErr exact_0(dInt dim,dReal x[],dScalar val[]);
  
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
  err = dJacobiSetDegrees(jac,15,4);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);
#if JACOBI_VIEW 
  err = dJacobiView(jac,viewer);dCHK(err);
#endif
  err = checkRulesAndEFS(jac);dCHK(err);  
  err = dJacobiDestroy(jac);dCHK(err);
  err = PetscFinalize();dCHK(err);
  dFunctionReturn(0);
}

dErr checkRulesAndEFS(dJacobi jac)
{
  const dInt N = 10;
  const dInt rsize[3] = {12,12,12};
  const dInt bsize[3] = {6,6,6};
  const dInt dim = 3;
  const dTopology topo = iMesh_HEXAHEDRON;
  PetscViewer viewer = PETSC_VIEWER_STDOUT_WORLD;
  MPI_Comm comm = ((PetscObject)jac)->comm;
  Vec u;
  dRule *rule;
  dEFS *efs;
  dInt index;
  void **rbase,**ebase;
  dErr err;

  dFunctionBegin;
  err = dMalloc(N*sizeof(dRule),&rule);dCHK(err);

  /* find the size of all rules on the elements, the allocate buffer and fill with rules */
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

#if RULE_VIEW
  for (dInt i=0; i<N; i++) {
    err = dPrintf(comm,"Rule for element %d\n",i);dCHK(err);
    err = dRuleView(rule[i],viewer);dCHK(err);
  }
#endif

  /* Same as above but for dEFS.  Note that these are all the same order, but that could be changed by changing 'bsize'
  * and different topology can be handled by changing 'topo' */
  err = dMalloc(N*sizeof(dEFS),&efs);dCHK(err);
  index = 0; ebase = NULL;
  do {
    if (index) {
      err = dMalloc(index*sizeof(*ebase),&ebase);dCHK(err);
      index = 0;
    }
    for (dInt i=0; i<N; i++) {
      err = dJacobiGetEFS(jac,topo,bsize,rule[i],&efs[i],ebase,&index);dCHK(err);
    }
  } while (!ebase);

#if EFS_VIEW
  for (dInt i=0; i<N; i++) {
    err = dPrintf(comm,"EFS for element %d\n",i);dCHK(err);
    err = dEFSView(efs[i],viewer);dCHK(err);
  }
#endif

  err = createFS(comm,N,dim,efs,&u);dCHK(err);
  err = checkInterp(N,efs,u);dCHK(err);

  err = dFree(rule);dCHK(err);
  err = dFree(efs);dCHK(err);
  err = dFree(rbase);dCHK(err);
  err = dFree(ebase);dCHK(err);
  err = VecDestroy(u);dCHK(err);
  dFunctionReturn(0);
}

dErr createFS(MPI_Comm comm,dInt dim,dInt N,dEFS *efs,Vec *U)
{
  Vec u;
  dInt m,n;
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,4);
  dValidPointer(U,5);
  *U = 0;
  n = 0;
  for (dInt i=0; i<N; i++) {
    err = dEFSGetSizes(efs[i],NULL,NULL,&m);dCHK(err);
    n += m;
  }
  err = VecCreate(comm,&u);dCHK(err);
  err = VecSetSizes(u,n*dim,PETSC_DECIDE);dCHK(err);
  err = VecSetFromOptions(u);dCHK(err);
  *U = u;
  dFunctionReturn(0);
}

dErr checkInterp(dInt N,dEFS efs[],Vec u)
{
  dRule rule;
  dReal *x[3],y[3],w,z;
  const dReal *qx[3],*qw[3];
  dScalar *f,*g=NULL,*work=NULL,h[3];
  dInt P[3],Q[3],ind,size,wsize=0,dim,needed,qind,gsize=0;
  dErr err;

  dFunctionBegin;
  {
    dInt m,M;
    err = VecGetLocalSize(u,&m);dCHK(err);
    err = VecGetSize(u,&M);dCHK(err);
    err = dPrintf(PETSC_COMM_WORLD,"vec sizes local=%d global=%d\n",m,M);dCHK(err);
  }
  err = VecGetArray(u,&f);dCHK(err);
  ind = 0;
  for (dInt i=0; i<N; i++) {
    err = dEFSGetTensorNodes(efs[i],&dim,P,x);dCHK(err);
    switch (dim) {
      case 1: P[1] = 1; P[2] = 1; break;
      case 2: P[2] = 1; break;
      case 3: break;
      default: dERROR(1,"dim %d out of range",dim);
    }
    for (dInt j=0; j<P[0]; j++) {
      for (dInt k=0; k<P[1]; k++) {
        for (dInt l=0; l<P[2]; l++) {
          y[0] = x[0][j];
          y[1] = (dim > 1) ? x[1][k] : 0.0;
          y[2] = (dim > 2) ? x[2][l] : 0.0;
          err = exact_0(dim,y,&f[ind]);dCHK(err);
          ind += dim;
        }
      }
    }
  }
  err = VecRestoreArray(u,&f);dCHK(err);

  ind = 0;
  err = VecGetArray(u,&f);dCHK(err);
  for (dInt i=0; i<N; i++) {
    err = dEFSGetSizes(efs[i],&dim,NULL,&size);dCHK(err);
    err = dEFSGetRule(efs[i],&rule);dCHK(err);
    err = dRuleGetTensorNodeWeight(rule,&dim,Q,qx,qw);dCHK(err);
    needed = dim * productInt(dim,Q);
    if (needed > gsize) {
      err = dFree(g);dCHK(err);
      gsize = needed * 2;
      err = dMalloc(gsize*sizeof g[0],&g);dCHK(err);
    }

    err = dMemzero(g,gsize*sizeof g[0]);dCHK(err);
    err = dEFSApply(efs[i],dim,&wsize,&work,&f[ind],g,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);

    /* compare to exact solution */
    switch (dim) {
      case 1: Q[1] = 1; Q[2] = 1; break;
      case 2: Q[2] = 1; break;
      case 3: break;
      default: dERROR(1,"dim %d out of range",dim);
    }
    qind = 0; z = 0.0;
    for (dInt j=0; j<Q[0]; j++) {
      for (dInt k=0; k<Q[1]; k++) {
        for (dInt l=0; l<Q[2]; l++) {
          y[0] = qx[0][j];
          y[1] = (dim > 1) ? qx[1][k] : 0.0;
          y[2] = (dim > 2) ? qx[2][l] : 0.0;
          err = exact_0(dim,y,h);dCHK(err);
          w = qw[0][j] * (dim>1 ? qw[1][k] : 1.0) * (dim>2 ? qw[2][l] : 1.0);
          for (dInt d=0; d<dim; d++) {
            z += dSqr(g[qind+d] - h[d]) * w;
          }
          qind += dim;
        }
      }
    }
    printf("element %d\n",i);
    printf("L2 error for element %3d = %g\n",i,z); 
    ind += dim*size;
  }
  err = dFree(g);dCHK(err);
  err = dFree(work);dCHK(err);
  dFunctionReturn(0);
}

dErr exact_0(dInt dim,dReal x[],dScalar f[])
{

  dFunctionBegin;
  switch (dim) {
    case 1:
      f[0] = sin(x[0]);
      break;
    case 2:
      f[0] = sin(x[0]) * cos(x[1]);
      f[1] = cos(x[0]) * sin(x[1]);
      break;
    case 3:
      f[0] = sin(x[0]) * cos(x[1]) * tanh(x[2]);
      f[1] = cos(x[0]) * sin(x[1]) * tanh(x[2]);
      f[2] = sinh(x[0]) * tanh(x[1]) * cosh(x[2]);
      break;
    default:
      dERROR(1,"dimension %d not in range",dim);
  }
  dFunctionReturn(0);
}

dInt productInt(dInt N,const dInt a[])
{
  dInt z=1;
  for (dInt i=0; i<N; i++) z *= a[i];
  return z;
}
