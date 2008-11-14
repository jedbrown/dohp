static const char help[] = "Tests the dJacobi object.";

#include "dohpjacobi.h"
#include "private/jacimpl.h"
#include "petscvec.h"
#include <stdlib.h>

#define JACOBI_VIEW 0
#define RULE_VIEW 0
#define EFS_VIEW 0

static struct {
  dErr (*function)(dInt,const dReal[],dScalar[]);
  dErr (*deriv)(dInt,const dReal[],dScalar[]);
} exact;

static dErr exact_0(dInt dim,const dReal x[],dScalar f[])
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
static dErr dexact_0(dInt dim,const dReal x[],dScalar f[])
{

  dFunctionBegin;
#define dsin(a) cos(a)
#define dcos(a) (-1.0 * sin(a))
#define dsinh(a) cosh(a)
#define dcosh(a) sinh(a)
#define dtanh(a) (1.0 - dSqr(tanh(a)))
  switch (dim) {
    case 1:
      f[0] = dsin(x[0]);
      break;
    case 2:
      f[0] = dsin(x[0]) *  cos(x[1]);
      f[1] =  sin(x[0]) * dcos(x[1]);
      f[2] = dcos(x[0]) *  sin(x[1]);
      f[3] =  cos(x[0]) * dsin(x[1]);
      break;
    case 3:
      f[0] = dsin(x[0]) *  cos(x[1]) *  tanh(x[2]);
      f[1] =  sin(x[0]) * dcos(x[1]) *  tanh(x[2]);
      f[2] =  sin(x[0]) *  cos(x[1]) * dtanh(x[2]);
      f[3] = dcos(x[0]) *  sin(x[1]) *  tanh(x[2]);
      f[4] =  cos(x[0]) * dsin(x[1]) *  tanh(x[2]);
      f[5] =  cos(x[0]) *  sin(x[1]) * dtanh(x[2]);
      f[6] = dsinh(x[0]) *  tanh(x[1]) *  cosh(x[2]);
      f[7] =  sinh(x[0]) * dtanh(x[1]) *  cosh(x[2]);
      f[8] =  sinh(x[0]) *  tanh(x[1]) * dcosh(x[2]);
      break;
  }
#undef dsin
#undef dcos
#undef dsinh
#undef dcosh
#undef dtanh
  dFunctionReturn(0);
}

static dErr exact_1(dInt D,const dReal x[],dScalar f[])
{

  dFunctionBegin;
  switch (D) {
    case 3:
      f[0] = 0.5*dSqr(x[0]);
      f[1] = x[1];
      f[2] = x[2];
      break;
    default: dERROR(1,"not implemented");
  }
  dFunctionReturn(0);
}

static dErr dexact_1(dInt D,const dReal x[],dScalar f[])
{

  dFunctionBegin;
  switch (D) {
    case 3:
      f[0] = x[0];
      f[1] = 0;
      f[2] = 0;
      f[3] = 0;
      f[4] = 1;
      f[5] = 0;
      f[6] = 0;
      f[7] = 0;
      f[8] = 1;
      break;
    default: dERROR(1,"not implemented");
  }
  dFunctionReturn(0);
}


static dInt productInt(dInt N,const dInt a[])
{
  dInt z=1;
  for (dInt i=0; i<N; i++) z *= a[i];
  return z;
}

static dErr checkInterp(dInt N,s_dEFS efs[],Vec u)
{
  dRule rule;
  dReal *x[3],y[3],w,z,dz;
  const dReal *qx[3],*qw[3];
  dScalar *f,*g=NULL,*dg=NULL,*work=NULL,h[3],dh[9];
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
    err = dEFSGetTensorNodes(&efs[i],&dim,P,x);dCHK(err);
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
          {
            dInt vsize;
            err = VecGetLocalSize(u,&vsize);dCHK(err);
            if (vsize < ind + dim)
              dERROR(1,"Will overwrite Vec size %d: at %d;%d,%d,%d ind=%d D=%d\n",vsize,i,j,k,l,ind,dim);
          }
          err = exact.function(dim,y,&f[ind]);dCHK(err);
          ind += dim;
        }
      }
    }
  }

  ind = 0;
  for (dInt i=0; i<N; i++) {
    err = dEFSGetSizes(&efs[i],&dim,NULL,&size);dCHK(err);
    err = dEFSGetRule(&efs[i],&rule);dCHK(err);
    err = dRuleGetTensorNodeWeight(rule,&dim,Q,qx,qw);dCHK(err);
    needed = dim * productInt(dim,Q);
    if (needed > gsize) {
      err = dFree2(g,dg);dCHK(err);
      gsize = needed * 2;
      err = dMallocA2(gsize,&g,dim*gsize,&dg);dCHK(err);
    }

    err = dMemzero(g,gsize*sizeof g[0]);dCHK(err);
    for (dInt j=0; j<gsize; j++) { g[j] = NAN; }
    for (dInt j=0; j<gsize*dim; j++) { dg[j] = NAN; }
    err = dEFSApply(&efs[i],dim,&wsize,&work,f+ind,g,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSApply(&efs[i],dim,&wsize,&work,f+ind,dg,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);

    /* compare to exact solution */
    switch (dim) {
      case 1: Q[1] = 1; Q[2] = 1; break;
      case 2: Q[2] = 1; break;
      case 3: break;
      default: dERROR(1,"dim %d out of range",dim);
    }
    qind = 0; z = dz = 0.0;
    for (dInt j=0; j<Q[0]; j++) {
      for (dInt k=0; k<Q[1]; k++) {
        for (dInt l=0; l<Q[2]; l++) {
          y[0] = qx[0][j];
          y[1] = (dim > 1) ? qx[1][k] : 0.0;
          y[2] = (dim > 2) ? qx[2][l] : 0.0;
          err = exact.function(dim,y,h);dCHK(err);
          err = exact.deriv(dim,y,dh);dCHK(err);
          w = qw[0][j] * (dim>1 ? qw[1][k] : 1.0) * (dim>2 ? qw[2][l] : 1.0);
          for (dInt d=0; d<dim; d++) {
            z += dSqr(g[qind+d] - h[d]) * w;
            for (dInt e=0; e<dim; e++) {
              /* The gradient function gives us chunks with each column of the Jacobian */
              const dInt chunk = dim*Q[0]*Q[1]*Q[2],qii = qind + d;
              dReal l2error;
              if (e*chunk+qii >= gsize*dim) {
                dERROR(1,"About to make invalid read");
              }
              l2error = dSqr(dg[e*chunk+qii] - dh[d*dim+e]);
              dz += l2error * w;
            }
          }
          qind += dim;
        }
      }
    }
    {
      dInt D;
      err = dEFSGetTensorNodes(&efs[i],&D,P,NULL);dCHK(err);
      err = dPrintf(PETSC_COMM_WORLD,"L2 error for element %3d (%2d,%2d,%2d) with rule (%2d,%2d,%2d) interp %8.1e, derivative %8.1e\n",i,P[0],P[1],P[2],Q[0],Q[1],Q[2],z,dz);dCHK(err);
    }
    ind += dim*size;
  }
  err = VecRestoreArray(u,&f);dCHK(err);
  err = dFree2(g,dg);dCHK(err);
  err = dFree(work);dCHK(err);
  dFunctionReturn(0);
}

static dErr createFS(MPI_Comm comm,dInt dim,dInt N,dEFS efs,Vec *U)
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
    err = dEFSGetSizes(&efs[i],NULL,NULL,&m);dCHK(err);
    n += m * dim;
  }
  err = VecCreate(comm,&u);dCHK(err);
  err = VecSetSizes(u,n,PETSC_DECIDE);dCHK(err);
  err = VecSetFromOptions(u);dCHK(err);
  *U = u;
  dFunctionReturn(0);
}

static dErr checkRulesAndEFS(dJacobi jac)
{
  const dInt dim = 3;
  PetscViewer viewer = PETSC_VIEWER_STDOUT_WORLD;
  MPI_Comm comm = ((PetscObject)jac)->comm;
  Vec u;
  dInt N,minrdeg,maxrdeg,minbdeg,maxbdeg,*rdeg,*bdeg;
  dEntTopology *topo;
  dTruth showrules,showefs;
  s_dRule *rule;
  s_dEFS *efs;
  dErr err;

  dFunctionBegin;
  minrdeg = 10;
  maxrdeg = 10;
  minbdeg = 2;
  maxbdeg = 6;
  showrules = dFALSE;
  showefs = dFALSE;
  err = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Options for testing Jacobi",NULL);dCHK(err);
  err = PetscOptionsInt("-min_rdeg","Minimum rule degree",NULL,minrdeg,&minrdeg,NULL);dCHK(err);
  err = PetscOptionsInt("-max_rdeg","Maximum rule degree",NULL,maxrdeg,&maxrdeg,NULL);dCHK(err);
  err = PetscOptionsInt("-min_bdeg","Minimum basis degree",NULL,minbdeg,&minbdeg,NULL);dCHK(err);
  err = PetscOptionsInt("-max_bdeg","Maximum basis degree",NULL,maxbdeg,&maxbdeg,NULL);dCHK(err);
  err = PetscOptionsTruth("-show_rules","Show rules",NULL,showrules,&showrules,NULL);dCHK(err);
  err = PetscOptionsTruth("-show_efs","Show EFS",NULL,showefs,&showefs,NULL);dCHK(err);
  err = PetscOptionsEnd();dCHK(err);
  N = (maxrdeg-minrdeg+1) * (maxbdeg-minbdeg+1);
  err = dMallocA5(N,&topo,3*N,&rdeg,3*N,&bdeg,N,&rule,N,&efs);dCHK(err);
  for (dInt i=minrdeg; i<=maxrdeg; i++) {
    for (dInt j=minbdeg; j<=maxbdeg; j++) {
      const dInt e = (i-minrdeg)*(maxbdeg-minbdeg+1) + j-minbdeg;
      rdeg[3*e+0] = rdeg[3*e+1] = rdeg[3*e+2] = i;
      bdeg[3*e+0] = bdeg[3*e+1] = bdeg[3*e+2] = j;
      topo[e] = dTOPO_HEX;
    }
  }

  err = dJacobiGetRule(jac,N,topo,rdeg,rule);dCHK(err);
  if (showrules) {
    for (dInt i=0; i<N; i++) {
      err = dPrintf(comm,"Rule for element %d\n",i);dCHK(err);
      err = dRuleView(&rule[i],viewer);dCHK(err);
    }
  }

  err = dJacobiGetEFS(jac,N,topo,bdeg,rule,efs);dCHK(err);
  if (showefs) {
    for (dInt i=0; i<N; i++) {
      err = dPrintf(comm,"EFS for element %d\n",i);dCHK(err);
      err = dEFSView(&efs[i],viewer);dCHK(err);
    }
  }

  err = createFS(comm,dim,N,efs,&u);dCHK(err);
  err = checkInterp(N,efs,u);dCHK(err);

  err = dFree5(topo,rdeg,bdeg,rule,efs);dCHK(err);
  err = VecDestroy(u);dCHK(err);
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  dJacobi jac;
  MPI_Comm comm;
  PetscViewer viewer;
  dInt ex;
  dErr err;

  dFunctionBegin;
  err = PetscInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;

  ex = 1;
  err = PetscOptionsBegin(comm,NULL,"Jacobi ex1 test driver options",NULL);dCHK(err);
  err = PetscOptionsInt("-exact","exact solution number",NULL,ex,&ex,NULL);dCHK(err);
  err = PetscOptionsEnd();dCHK(err);
  switch (ex) {
    case 0:
      exact.function = exact_0;
      exact.deriv = dexact_0;
      break;
    case 1:
      exact.function = exact_1;
      exact.deriv = dexact_1;
      break;
    default: dERROR(1,"exact solution %d not implemented",ex);
  }
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
