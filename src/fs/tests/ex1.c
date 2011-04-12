static const char help[] = "Test the construction of dFS objects and anisotropic propogation.\n";

#include <petscsnes.h>
#include <dohpfs.h>
#include <dohpmesh.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>

#define ALEN(a) (dInt)(sizeof(a)/sizeof((a)[0]))

struct Options {
  dInt constBDeg;
  dBool  showsoln;
  dInt cycles;
  dInt  proj_version;
  dQuadratureMethod proj_qmethod;
  dQuadratureMethod jac_qmethod;
  dReal q1scale;
  dReal frequency[3];
  dReal normRequirePtwise[3];
  dReal normRequireGrad[3];
} gopt;

/* A global variable for the exact solution */
struct {
  dErr (*function)(const dReal[3],dScalar[1]);
  dErr (*gradient)(const dReal[3],dScalar[1][3]);
} exact;

static dErr exact_0_function(const dReal x[3],dScalar f[1])
{
  const dReal a = gopt.frequency[0],b = gopt.frequency[1],c = gopt.frequency[2];
  f[0] = sin(x[0])*cosh(b*x[1]) + cos(a*x[0])*exp(x[1]) + sinh(x[1])*tanh(c*x[2]);
  return 0;
}
static dErr exact_0_gradient(const dReal x[3],dScalar df[1][3])
{
  const dReal a = gopt.frequency[0],b = gopt.frequency[1],c = gopt.frequency[2];
  df[0][0] = cos(x[0])*cosh(b*x[1]) - a*sin(a*x[0])*exp(x[1]);
  df[0][1] = sin(x[0])*b*sinh(b*x[1]) + cos(a*x[0])*exp(x[1]) + cosh(x[1])*tanh(c*x[2]);
  df[0][2] = sinh(x[1]) * c*(1.0 - dSqr(tanh(c*x[2])));
  return 0;
}

static dErr exact_1_function(const dReal x[3],dScalar f[1])
{ f[0] = 10*x[0] + 1*x[1] + 0.1*x[2]; return 0; }
static dErr exact_1_gradient(const dUNUSED dReal x[3],dScalar df[1][3])
{ df[0][0] = 10; df[0][1] = 1; df[0][2] = 0.1; return 0; }

static dErr createHexMesh(iMesh_Instance mi)
{
  static const int vtxlen = 12*3;
  static const double vtx_affine[12*3] = {0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                          0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                          2,0,0, 2,1,0, 2,1,1, 2,0,1};
  static const double vtx_parametric[12*3] = {0.1,0.2,0.3, 1.3,0.2,0.1, 1.1,1.2,0.3, 0.3,1.2,0.1,
                                              0.1,0.2,1.3, 1.3,0.2,1.1, 1.1,1.2,1.3, 0.3,1.2,1.1,
                                              2.1,0.2,0.3, 2.3,1.2,0.1, 2.1,1.2,1.3, 2.3,0.2,1.1};
  const double *vtx;
  dBool  iaffine = dFALSE;
  int rconn[16] = {0,1,2,3,4,5,6,7, 1,2,6,5,8,9,10,11};
  int fconn[11*4] = {0,1,2,3, 1,2,6,5, 2,3,7,6, 0,3,7,4, 0,1,5,4, 5,6,7,4,
                     8,11,5,1, 2,1,8,9, 8,9,10,11, 11,10,6,5, 9,10,6,2};
  int econn[20*2] = {0,1,4,5,7,6,3,2, 1,8,5,11,6,10,2,9,
                     0,3,3,7,7,4,4,0, 1,2,2,6,6,5,5,1, 8,11,11,10,10,9,9,8};
  iBase_EntityHandle work[100];
  MeshListEH v=MLZ,e=MLZ,f=MLZ,r=MLZ,tv=MLZ;
  MeshListInt stat=MLZ,off=MLZ;
  dIInt ierr;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Hex mesh options",NULL);dCHK(err);
  err = PetscOptionsBool("-affine","Use an affine coordinate map",NULL,iaffine,&iaffine,NULL);dCHK(err);
  err = PetscOptionsEnd();dCHK(err);
  vtx = iaffine ? vtx_affine : vtx_parametric;
  iMesh_createVtxArr(mi,vtxlen/3,iBase_INTERLEAVED,vtx,vtxlen,MLREF(v),&ierr);dICHK(mi,ierr);
  for (int i=0; i<ALEN(rconn); i++) work[i] = v.v[rconn[i]];
  iMesh_createEntArr(mi,iMesh_HEXAHEDRON,work,ALEN(rconn),MLREF(r),MLREF(stat),&ierr);dICHK(mi,ierr);
  {                             /* Check to see if any orientations changed */
    iMesh_getEntArrAdj(mi,r.v,r.s,iMesh_POINT,MLREF(tv),MLREF(off),&ierr);dICHK(mi,ierr);
    if (tv.s != v.s + 4) dERROR(PETSC_COMM_SELF,1,"wrong number of vertices returned"); /* interface verts counted twice */
    for (int i=0; i<tv.s; i++) {
      if (v.v[rconn[i]] != tv.v[i]) dERROR(PETSC_COMM_SELF,1,"unexpected vertex ordering");
    }
    MeshListFree(tv);
    MeshListFree(off);
  }
  MeshListFree(r); MeshListFree(stat);
  for (int i=0; i<ALEN(fconn); i++) work[i] = v.v[fconn[i]];
  iMesh_createEntArr(mi,iMesh_QUADRILATERAL,work,ALEN(fconn),MLREF(f),MLREF(stat),&ierr);dICHK(mi,ierr);
  MeshListFree(f); MeshListFree(stat);
  for (int i=0; i<ALEN(econn); i++) work[i] = v.v[econn[i]];
  iMesh_createEntArr(mi,iMesh_LINE_SEGMENT,work,ALEN(econn),MLREF(e),MLREF(stat),&ierr);dICHK(mi,ierr);
  MeshListFree(e); MeshListFree(stat);
  MeshListFree(v);
  dFunctionReturn(0);
}

/* Verify the adjacencies against these hand-verified values, should prevent regressions */
static dErr verifyAdjacencies(dMesh mesh)
{
  const dInt nents = 45,toff[5] = {0,12,32,43,45};
  const dInt adjoff[46] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,2,4, 6,8,10,12,14, 16,18,20,22,24, 26,28,30,32,34,
                           36,38,40,44,48, 52,56,60,64,68, 72,76,80,84,90, 96};
  const dInt adjind[96] = {0,1,4,5,7, 6,3,2,1,8, 5,11,6,10,2, 9,0,3,3,7,
                           7,4,4,0,1, 2,2,6,6,5, 5,1,8,11,11, 10,10,9,9,8,
                           12,24,15,20,24, 25,26,27,15,21, 14,25,20,21,22, 23,12,27,13,23,
                           26,14,22,13,28, 17,27,16,24,16, 31,19,31,30,29, 28,29,18,26,17,
                           30,18,25,19,36, 33,34,35,32,37, 39,42,41,38,33, 40};
  const dInt adjperm[96] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
                            0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
                            0,0,1,1,0, 0,0,0,1,0, 0,1,0,0,0, 0,0,1,1,0,
                            1,1,0,0,0, 1,0,0,1,0, 1,1,1,1,1, 1,0,1,0,0,
                            1,1,1,0,0, 0,0,5,4,1, 5,7,2,2,4, 0};
  dMeshAdjacency ma;
  dInt i;
  dErr err;

  dFunctionBegin;
  err = dMeshGetAdjacency(mesh,0,&ma);dCHK(err);
  dASSERT(ma->nents == nents);
  for (i=0; i<ALEN(toff); i++) dASSERT(ma->toff[i] == toff[i]);
  for (i=0; i<ALEN(adjoff); i++) dASSERT(ma->adjoff[i] == adjoff[i]);
  for (i=0; i<ALEN(adjind); i++) dASSERT(ma->adjind[i] == adjind[i]);
  for (i=0; i<ALEN(adjperm); i++) dASSERT(ma->adjperm[i] == adjperm[i]);
  err = dMeshRestoreAdjacency(mesh,0,&ma);dCHK(err);
  dFunctionReturn(0);
}

static dErr tagHexes(dMesh mesh,dMeshTag *intag)
{
  dMeshTag tag;
  dMeshEH ents[100],hex[2];
  dInt nhex,nents;
  dPolynomialOrder adata[ALEN(ents)],hexDegree[ALEN(hex)] = {dPolynomialOrderCreate(0,3,4,5),dPolynomialOrderCreate(0,6,7,8)};
  dErr err;

  dFunctionBegin;
  *intag = 0;
  if (gopt.constBDeg) {
    for (dInt i=0; i<ALEN(hexDegree); i++) hexDegree[i] = dPolynomialOrderCreate(0,gopt.constBDeg,gopt.constBDeg,gopt.constBDeg);
  }
  err = dMeshTagCreate(mesh,"ex1_anisotropic",1,dDATA_INT,&tag);dCHK(err);
  /* tag edges and faces with high values, currently needed to propagate degrees (will overwrite high values) */
  for (dEntType type=dTYPE_VERTEX; type<=dTYPE_FACE; type++) {
    err = dMeshGetEnts(mesh,0,type,dTOPO_ALL,ents,ALEN(ents),&nents);
    for (dInt i=0; i<nents; i++) adata[i] = dPolynomialOrderCreate(0,30,30,30);
    err = dMeshTagSetData(mesh,tag,ents,nents,adata,nents,dDATA_INT);dCHK(err);
  }
  /* tag the hexes with meaningful values */
  err = dMeshGetEnts(mesh,0,dTYPE_ALL,dTOPO_HEX,hex,ALEN(hex),&nhex);dCHK(err);
  if (nhex != 2) dERROR(PETSC_COMM_SELF,1,"wrong number of hexes");
  err = dMeshTagSetData(mesh,tag,hex,nhex,hexDegree,nhex,dDATA_INT);dCHK(err);
  *intag = tag;
  dFunctionReturn(0);
}

static dErr examine(dMesh mesh,dMeshTag tag)
{
  dMeshEH ents[100];
  dEntTopology topo[100];
  MeshListEH conn=MLZ;
  MeshListInt connoff=MLZ;
  dInt data[3*ALEN(ents)],nents,nv;
  dIInt ierr;
  iMesh_Instance mi;
  dErr err;

  dFunctionBegin;
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  for (dEntType type=dTYPE_EDGE; type<dTYPE_ALL; type++) {
    err = dMeshGetEnts(mesh,0,type,dTOPO_ALL,ents,ALEN(ents),&nents);dCHK(err);
    err = dMeshGetTopo(mesh,nents,ents,topo);dCHK(err);
    iMesh_getEntArrAdj(mi,ents,nents,iBase_VERTEX,MLREF(conn),MLREF(connoff),&ierr);dICHK(mi,ierr);
    err = dMeshTagGetData(mesh,tag,ents,nents,data,ALEN(data),dDATA_INT);dCHK(err);
    for (dInt i=0; i<nents; i++) {
      switch (topo[i]) {
        case dTOPO_HEX: nv = 8; break;
        case dTOPO_QUAD: nv = 4; break;
        case dTOPO_LINE: nv = 2; break;
        default: dERROR(PETSC_COMM_SELF,1,"unsupported topology");
      }
      err = dPrintf(PETSC_COMM_WORLD,"%30s %d %d %d   [",dMeshEntTopologyName(topo[i]),data[3*i],data[3*i+1],data[3*i+2]);dCHK(err);
      for (dInt j=0; j<nv; j++) { err = dPrintf(PETSC_COMM_WORLD," %p",conn.v[connoff.v[i]+j]);dCHK(err); }
      err = dPrintf(PETSC_COMM_WORLD," ]\n");dCHK(err);
    }
    MeshListFree(conn); MeshListFree(connoff);
  }
  dFunctionReturn(0);
}

static dErr useFS(dFS fs)
{
  Vec x,y,g;
  dScalar *xx;
  dInt nx=-1;
  dErr err;

  dFunctionBegin;
  err = dFSCreateGlobalVector(fs,&g);dCHK(err);
  err = dFSCreateExpandedVector(fs,&x);dCHK(err);
  err = VecGetSize(x,&nx);dCHK(err);
  {
    dInt n,rstart,rend;
    err = VecGetSize(g,&n);dCHK(err);
    err = VecGetOwnershipRange(g,&rstart,&rend);dCHK(err);
    err = PetscPrintf(PETSC_COMM_WORLD,"expanded size %d, global size %d, range %d .. %d \n",nx,n,rstart,rend);dCHK(err);
  }
  err = VecDuplicate(x,&y);dCHK(err);
  err = VecGetArray(x,&xx);dCHK(err);
  for (dInt i=0; i<nx; i++) {
    xx[i] = 1000 + 1.0 * i;
  }
  err = VecRestoreArray(x,&xx);dCHK(err);
  err = VecSet(x,1.0);dCHK(err);
  err = dFSExpandedToGlobal(fs,x,g,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = dFSGlobalToExpanded(fs,g,y,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  {
    dScalar xsum,ysum,gsum;
    err = VecSum(x,&xsum);dCHK(err);
    err = VecSum(g,&gsum);dCHK(err);
    err = VecSum(y,&ysum);dCHK(err);
    err = PetscPrintf(PETSC_COMM_WORLD,"|x| = %f, |g| = %f, |y| = %f\n",xsum,gsum,ysum);dCHK(err);
    if (dAbs(xsum-gsum) > 1e-14) dERROR(PETSC_COMM_SELF,1,"Expanded sum does not match global sum");
    /* There are 16 points on the interface between the elements, these points get double-counted in both elements, so
    * the expanded sum is larger than the global sum by 32. */
    if (gopt.constBDeg) {
      if (dAbs(ysum-gsum-2.0*dSqr(gopt.constBDeg+1)) > 1e-14) dERROR(PETSC_COMM_SELF,1,"Unexpected expanded sum %f != %f + 2*%D^2",ysum,gsum,gopt.constBDeg+1);
    } else {
      dERROR(PETSC_COMM_SELF,1,"Don't know how to check for non-const Basis Degree");
    }
  }
  err = VecDestroy(x);dCHK(err);
  err = VecDestroy(y);dCHK(err);
  err = VecDestroy(g);dCHK(err);
  dFunctionReturn(0);
}

static dErr confirmSufficientResiduals(const char name[],const dReal resnorm[3],const dReal required[3])
{
  static const char *normName[3] = {"1","2","infty"};
  dFunctionBegin;
  for (dInt i=0; i<3; i++) {
    if (required[i] > 0 && resnorm[i] > required[i])
      dERROR(PETSC_COMM_SELF,1,"Norm requirement not met: ||%s||_%s = %g > %g",name,normName[i],resnorm[i],required[i]);
  }
  dFunctionReturn(0);
}

struct ProjContext {
  dFS fs;
  Vec x,y;
  dRuleset ruleset[dQUADRATURE_METHOD_INVALID];
  dRulesetIterator regioniter[dQUADRATURE_METHOD_INVALID];
};

static dErr ProjGetRuleset(struct ProjContext *proj,dQuadratureMethod qmethod,dRuleset *ruleset)
{
  dErr err;

  dFunctionBegin;
  if (!proj->ruleset[qmethod]) {
    dMeshESH domain;
    err = dFSGetDomain(proj->fs,&domain);dCHK(err);
    err = dFSGetPreferredQuadratureRuleSet(proj->fs,domain,dTYPE_REGION,dTOPO_ALL,qmethod,&proj->ruleset[qmethod]);dCHK(err);
  }
  *ruleset = proj->ruleset[qmethod];
  dFunctionReturn(0);
}

static dErr ProjGetRegionIterator(struct ProjContext *proj,dQuadratureMethod qmethod,dRulesetIterator *riter)
{
  dErr err;

  dFunctionBegin;
  if (!proj->regioniter[qmethod]) {
    dRulesetIterator iter;
    dRuleset ruleset;
    dFS cfs;
    err = ProjGetRuleset(proj,qmethod,&ruleset);dCHK(err);
    err = dFSGetCoordinateFS(proj->fs,&cfs);dCHK(err);
    err = dRulesetCreateIterator(ruleset,cfs,&iter);dCHK(err);
    err = dRulesetIteratorAddFS(iter,proj->fs);dCHK(err);
    proj->regioniter[qmethod] = iter;
  }
  *riter = proj->regioniter[qmethod];
  dFunctionReturn(0);
}

static dErr ProjResidual(dUNUSED SNES snes,Vec gx,Vec gy,void *ctx)
{
  dErr err;
  struct ProjContext *proj = ctx;
  dFS fs = proj->fs;
  dInt nelem,n;
  const dEFS *efs,*cefs;
  dScalar *y,*u,*v,*q,*cjac,*cjinv,*jw;
  const dScalar *x,*coords;
  dRuleset ruleset;
  dFS cfs;
  Vec Coords;

  dFunctionBegin;
  err = dFSGlobalToExpanded(fs,gx,proj->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArrayRead(proj->x,&x);dCHK(err);
  err = VecZeroEntries(proj->y);dCHK(err);
  err = VecGetArray(proj->y,&y);dCHK(err);
  err = ProjGetRuleset(proj,gopt.proj_qmethod,&ruleset);dCHK(err);
  err = dFSGetEFS(fs,ruleset,&nelem,&efs);dCHK(err);

  err = dFSGetCoordinateFS(fs,&cfs);dCHK(err);
  err = dFSGetGeometryVectorExpanded(fs,&Coords);dCHK(err);
  err = VecGetArrayRead(Coords,&coords);dCHK(err);
  err = dFSGetEFS(cfs,ruleset,&n,&cefs);dCHK(err);
  dASSERT(n == nelem);

  err = dRulesetGetWorkspace(ruleset,&q,&cjac,&cjinv,&jw,1,&u,&v,NULL,NULL,0);dCHK(err);
  for (dInt e=0,off=0,coff=0; e<nelem; e++) {
    dInt P,Q;
    dRule rule;
    err = dEFSApply(cefs[e],NULL,3,coords+coff,q,dAPPLY_INTERP,INSERT_VALUES);dCHK(err); /* Location of quadrature nodes */
    err = dEFSApply(cefs[e],NULL,3,coords+coff,cjac,dAPPLY_GRAD,INSERT_VALUES);dCHK(err); /* Coordinate Jacobian */
    err = dEFSGetRule(cefs[e],&rule);dCHK(err);
    err = dRuleComputePhysical(rule,cjac,cjinv,jw);dCHK(err);
    err = dEFSApply(efs[e],cjinv,1,x+off,u,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dRuleGetSize(rule,0,&Q);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar f[1];             /* Scalar problem */
      err = exact.function(&q[i*3],f);dCHK(err);
      v[i*1+0] = (u[i*1+0] - f[0]) * jw[i];
    }
    err = dEFSApply(efs[e],cjinv,1,v,y+off,dAPPLY_INTERP_TRANSPOSE,ADD_VALUES);dCHK(err);
    err = dEFSGetSizes(efs[e],NULL,NULL,&P);dCHK(err);
    off += P;
    err = dEFSGetSizes(cefs[e],NULL,NULL,&P);dCHK(err);
    coff += P*3;
  }
  err = dFSRestoreEFS(fs,ruleset,&n,&efs);dCHK(err);
  err = dFSRestoreEFS(cfs,ruleset,&n,&cefs);dCHK(err);
  err = dRulesetRestoreWorkspace(ruleset,&q,&cjac,&cjinv,&jw,1,&u,&v,NULL,NULL,0);dCHK(err);
  err = VecRestoreArrayRead(Coords,&coords);dCHK(err);
  err = VecRestoreArrayRead(proj->x,&x);dCHK(err);
  err = VecRestoreArray(proj->y,&y);dCHK(err);
  err = dFSExpandedToGlobal(fs,proj->y,gy,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  dFunctionReturn(0);
}

static dErr ProjResidual2(dUNUSED SNES snes,Vec gx,Vec gy,void *ctx)
{
  dErr err;
  struct ProjContext *proj = ctx;
  dFS fs = proj->fs;
  Vec Coords;
  dRulesetIterator iter;

  dFunctionBegin;
  err = VecZeroEntries(gy);dCHK(err);
  err = dFSGetGeometryVectorExpanded(fs,&Coords);dCHK(err);
  err = ProjGetRegionIterator(proj,gopt.proj_qmethod,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,gy,dFS_INHOMOGENEOUS);dCHK(err);

  while (dRulesetIteratorHasPatch(iter)) {
    dRule rule;
    dEFS cefs,efs;
    dScalar *ex,*eu,*ev;
    dScalar *cjinv,*jw,*x,*dx,*u,*du,*v,*dv;
    dInt Q;
    err = dRulesetIteratorGetElement(iter,&rule,&cefs,&ex,NULL,&efs,&eu,&ev);dCHK(err);
    err = dRulesetIteratorGetPatchSpace(iter,&cjinv,&jw,&x,&dx,NULL,NULL,&u,&du,&v,&dv);dCHK(err);
    err = dEFSApply(cefs,NULL,3,ex,x,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSApply(cefs,NULL,3,ex,dx,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
    err = dRuleComputePhysical(rule,dx,cjinv,jw);dCHK(err);
    err = dEFSApply(efs,cjinv,1,eu,u,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
    err = dEFSApply(efs,cjinv,1,eu,du,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
    err = dRuleGetSize(rule,0,&Q);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar f[1];             /* Scalar problem */
      err = exact.function(&x[i*3],f);dCHK(err);
      v[i*1+0] = (u[i*1+0] - f[0]) * jw[i];
      dv[i*3+0] = 0;
      dv[i*3+1] = 0;
      dv[i*3+2] = 0;
    }
    err = dEFSApply(efs,cjinv,1,v,ev,dAPPLY_INTERP_TRANSPOSE,INSERT_VALUES);dCHK(err);
    err = dEFSApply(efs,cjinv,1,dv,ev,dAPPLY_GRAD_TRANSPOSE,ADD_VALUES);dCHK(err);
    err = dRulesetIteratorRestorePatchSpace(iter,&cjinv,&jw,&x,&dx,NULL,NULL,&u,&du,&v,&dv);dCHK(err);
    err = dRulesetIteratorCommitPatch(iter,NULL,ev);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }

  err = dRulesetIteratorFinish(iter);dCHK(err);
  dFunctionReturn(0);
}

static dErr ProjResidual3(dUNUSED SNES snes,Vec gx,Vec gy,void *ctx)
{
  dErr err;
  struct ProjContext *proj = ctx;
  Vec Coords;
  dRulesetIterator iter;

  dFunctionBegin;
  err = VecZeroEntries(gy);dCHK(err);
  err = dFSGetGeometryVectorExpanded(proj->fs,&Coords);dCHK(err);
  err = ProjGetRegionIterator(proj,gopt.proj_qmethod,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,gy,dFS_INHOMOGENEOUS);dCHK(err);

  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar *x,*dx,*u,*du,*v,*dv;
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,&v,&dv);dCHK(err);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar f[1];             /* Scalar problem */
      err = exact.function(&x[i*3],f);dCHK(err);
      v[i*1+0] = (u[i*1+0] - f[0]) * jw[i];
      dv[i*3+0] = 0;
      dv[i*3+1] = 0;
      dv[i*3+2] = 0;
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES,NULL,NULL,v,dv);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }

  err = dRulesetIteratorFinish(iter);dCHK(err);
  dFunctionReturn(0);
}

static dErr ProjJacobian1(SNES dUNUSED snes,Vec gx,Mat *J,Mat *Jp,MatStructure *structure,void *ctx)
{
  struct ProjContext *proj = ctx;
  dRulesetIterator iter;
  dScalar *Kflat;
  dErr err;
  Vec Coords;

  dFunctionBegin;
  err = MatZeroEntries(*Jp);dCHK(err);
  err = dFSGetGeometryVectorExpanded(proj->fs,&Coords);dCHK(err);
  err = ProjGetRegionIterator(proj,gopt.jac_qmethod,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSplit(iter,NULL,NULL,NULL,&Kflat);dCHK(err);

  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw,*interp_flat,*deriv;
    const dInt *rowcol;
    dScalar *x,*dx,*u,*du;
    dInt Q,P;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,NULL,NULL);dCHK(err);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv);dCHK(err);
    {                           /* Scope so that we can declare new VLA pointers for convenient assembly */
      const dReal (*interp)[P] = (const dReal(*)[P])interp_flat;
      dScalar (*K)[P] = (dScalar(*)[P])Kflat;
      err = PetscMemzero(K,P*P*sizeof(K[0][0]));dCHK(err);
      for (dInt i=0; i<P; i++) {
        for (dInt q=0; q<Q; q++) {
          for (dInt j=0; j<P; j++) {
            K[i][j] += interp[q][i] * jw[q] * interp[q][j];
          }
        }
      }
      err = dFSMatSetValuesBlockedExpanded(proj->fs,*Jp,P,rowcol,P,rowcol,&K[0][0],ADD_VALUES);dCHK(err);
    }
    err = dRulesetIteratorRestorePatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }

  err = dRulesetIteratorFinish(iter);dCHK(err);

  err = MatAssemblyBegin(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  /* Dummy assembly, necessary to reset base vector in -snes_mf_operator */
  if (J != Jp) {
    err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  }
  *structure = SAME_NONZERO_PATTERN;
  dFunctionReturn(0);
}

#if 0
static dErr ProjJacobian(SNES dUNUSED snes,Vec gx,Mat dUNUSED *J,Mat *Jp,MatStructure *structure,void *ctx)
{
  struct ProjContext *proj = ctx;
  dRule *rule;
  dEFS *efs;
  dReal (*nx)[3];
  dScalar *x;
  dFS fs = proj->fs;
  dInt n,*off,*geomoff;
  dReal (*geom)[3];
  dErr err;

  dFunctionBegin;
  err = MatZeroEntries(*Jp);dCHK(err);
  err = dFSGlobalToExpanded(fs,gx,proj->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = VecGetArray(proj->x,&x);dCHK(err);
  err = dFSGetElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = dFSGetWorkspace(fs,__func__,&nx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err);
  for (dInt e=0; e<n; e++) {
    dInt three,P[3];
    err = dEFSGetGlobalCoordinates(efs[e],(const dReal(*)[3])(geom+geomoff[e]),&three,P,nx);dCHK(err);
    if (three != 3) dERROR(PETSC_COMM_SELF,1,"Dimension not equal to 3");
    for (dInt i=0; i<P[0]-1; i++) { /* P-1 = number of sub-elements in each direction */
      for (dInt j=0; j<P[1]-1; j++) {
        for (dInt k=0; k<P[2]-1; k++) {
          dQ1CORNER_CONST_DECLARE(c,rowcol,corners,off[e],nx,P,i,j,k);
          //const dScalar (*u)[1] = (const dScalar(*)[1])x+off[e]; /* Scalar valued function, can be indexed at corners as u[c[#]][0] */
          const dReal (*qx)[3],*jw,*flatBasis,*flatDeriv;
          dInt qn;
          dScalar K[8][8];
          err = dMemzero(K,sizeof(K));dCHK(err);
          err = dQ1HexComputeQuadrature(corners,&qn,&qx,&jw,&flatBasis,&flatDeriv);dCHK(err);
          {                     /* Scoping for VLA-pointer access */
            const dReal (*basis)[8] = (const dReal(*)[8])flatBasis;
            // const dReal (*deriv)[qn][8] = (const dReal(*)[qn][8])flatDeriv; /* UNUSED */
            if (qn != 8) dERROR(PETSC_COMM_SELF,1,"Unexpected number of quadrature points %d, but it *should* work, disable this error",qn);
            for (dInt lq=0; lq<qn; lq++) {           /* Loop over quadrature points */
              for (dInt ltest=0; ltest<8; ltest++) { /* Loop over test basis functions (corners) */
                for (dInt lp=0; lp<8; lp++) {        /* Loop over trial basis functions (corners) */
                  K[ltest][lp] += gopt.q1scale * (basis[lq][ltest] * jw[lq] * basis[lq][lp]);
                }
              }
            }
          }
          err = dFSMatSetValuesBlockedExpanded(fs,*Jp,8,rowcol,8,rowcol,&K[0][0],ADD_VALUES);dCHK(err);
        }
      }
    }
  }
  err = dFSRestoreWorkspace(fs,__func__,&nx,NULL,NULL,NULL,NULL,NULL,NULL);dCHK(err);
  err = dFSRestoreElements(fs,&n,&off,&rule,&efs,&geomoff,&geom);dCHK(err);
  err = VecRestoreArray(proj->x,&x);dCHK(err);
  err = MatAssemblyBegin(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*Jp,MAT_FINAL_ASSEMBLY);dCHK(err);
  /* Dummy assembly, somehow important for -snes_mf_operator */
  err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  *structure = DIFFERENT_NONZERO_PATTERN;
  dFunctionReturn(0);
}
#endif

static dErr ProjResidualNorms(struct ProjContext *proj,Vec gx,dReal residualNorms[static 3],dReal gresidualNorms[static 3])
{
  Vec Coords;
  dRulesetIterator iter;
  dErr err;

  dFunctionBegin;
  err = dMemzero(residualNorms,3*sizeof(residualNorms));dCHK(err);
  err = dMemzero(gresidualNorms,3*sizeof(gresidualNorms));dCHK(err);
  err = dFSGetGeometryVectorExpanded(proj->fs,&Coords);dCHK(err);
  err = ProjGetRegionIterator(proj,gopt.proj_qmethod,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter,Coords,dFS_INHOMOGENEOUS,NULL,gx,dFS_INHOMOGENEOUS,NULL);dCHK(err);

  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar *x,*dx,*u,*du;
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, &x,&dx,NULL,NULL, &u,&du,NULL,NULL);dCHK(err);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar f[1],df[1][3],r[1],gr[1][3],grsum;             /* Scalar problem */
      err = exact.function(&x[i*3],f);dCHK(err);
      err = exact.gradient(&x[i*3],df);dCHK(err);
      r[0] = u[i*1+0] - f[0];
      gr[0][0] = du[i*3+0] - df[0][0];
      gr[0][1] = du[i*3+1] - df[0][1];
      gr[0][2] = du[i*3+2] - df[0][2];
      grsum = dSqr(gr[0][0]) + dSqr(gr[0][1]) + dSqr(gr[0][2]);
      residualNorms[0] += dAbs(r[0]) * jw[i];               /* 1-norm */
      residualNorms[1] += dSqr(r[0]) * jw[i];               /* 2-norm */
      residualNorms[2] = dMax(residualNorms[2],dAbs(r[0])); /* Sup-norm */
      gresidualNorms[0] += grsum * jw[i];
      gresidualNorms[1] += dSqr(grsum) * jw[i];
      gresidualNorms[2] = dMax(gresidualNorms[2],grsum);
#if 0
      printf("pointwise stats %8g %8g %8g %8g\n",jw[i],r[0],dSqr(r[0]),residualNorms[1]);
      printf("pointwise grads %8g %8g %8g (%8g)\n",gr[0][0],gr[0][1],gr[0][2],grsum);
      printf("jinv[%2d][%3d]   %+3.1f %+3.1f %+3.1f    %+3.1f %+3.1f %+3.1f    %+3.1f %+3.1f %+3.1f\n",e,i,
             jinv[i][0][0],jinv[i][0][1],jinv[i][0][2],
             jinv[i][1][0],jinv[i][1][1],jinv[i][1][2],
             jinv[i][2][0],jinv[i][2][1],jinv[i][2][2]);
#endif
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  residualNorms[1] = dSqrt(residualNorms[1]);
  gresidualNorms[1] = dSqrt(gresidualNorms[1]);
  dFunctionReturn(0);
}

static dErr doProjection(dFS fs)
{
  struct ProjContext proj;
  MPI_Comm comm;
  SNES snes;
  Vec r,x;
  Mat J,Jp;
  dErr err;

  dFunctionBegin;
  err = dMemzero(&proj,sizeof(proj));dCHK(err);
  err = dObjectGetComm((dObject)fs,&comm);dCHK(err);
  proj.fs = fs;
  err = dFSCreateExpandedVector(fs,&proj.x);dCHK(err);
  err = VecDuplicate(proj.x,&proj.y);dCHK(err);

  err = dFSCreateGlobalVector(fs,&r);dCHK(err);
  err = dFSGetMatrix(fs,MATSEQAIJ,&Jp);dCHK(err);
  err = MatSetOptionsPrefix(Jp,"q1");dCHK(err);
  err = MatSeqAIJSetPreallocation(Jp,27,NULL);dCHK(err);
  J = Jp;                       /* Use -snes_mf_operator to apply J matrix-free instead of actually using Jp as the Krylov matrix */
  err = SNESCreate(comm,&snes);dCHK(err);
  switch (gopt.proj_version) {
  case 1:
    err = SNESSetFunction(snes,r,ProjResidual,(void*)&proj);dCHK(err);
    break;
  case 2:
    err = SNESSetFunction(snes,r,ProjResidual2,(void*)&proj);dCHK(err);
    break;
  case 3:
    err = SNESSetFunction(snes,r,ProjResidual3,(void*)&proj);dCHK(err);
    break;
  default: dERROR(PETSC_COMM_SELF,PETSC_ERR_USER,"Invalid residual version (-proj_version)");
  }
  err = SNESSetJacobian(snes,J,Jp,ProjJacobian1,(void*)&proj);dCHK(err);
  err = SNESSetTolerances(snes,PETSC_DEFAULT,1e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);dCHK(err);
  err = SNESSetFromOptions(snes);dCHK(err);
  err = VecDuplicate(r,&x);dCHK(err);
  for (dInt i=0; i<gopt.cycles; i++) {
    err = VecZeroEntries(x);dCHK(err);
    err = SNESSolve(snes,NULL,x);dCHK(err);
  }
  if (gopt.showsoln) {err = VecView(x,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);}
  {
    dReal norm[2],norminf,resNorms[3],gresNorms[3];
    err = VecNorm(r,NORM_1_AND_2,norm);dCHK(err);
    err = VecNorm(r,NORM_INFINITY,&norminf);dCHK(err);
    err = dPrintf(PETSC_COMM_WORLD,"Algebraic projection residual  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",norm[0],norm[1],norminf);dCHK(err);
    err = ProjResidualNorms(&proj,x,resNorms,gresNorms);dCHK(err);
    err = dPrintf(PETSC_COMM_WORLD,"Pointwise projection residual  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",resNorms[0],resNorms[1],resNorms[2]);dCHK(err);
    err = dPrintf(PETSC_COMM_WORLD,"Gradient ptwise proj residual  |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",gresNorms[0],gresNorms[1],gresNorms[2]);dCHK(err);
    err = confirmSufficientResiduals("Pointwise residuals",resNorms,gopt.normRequirePtwise);dCHK(err);
    err = confirmSufficientResiduals("Pointwise gradients",gresNorms,gopt.normRequireGrad);dCHK(err);
  }
  err = SNESDestroy(snes);dCHK(err);
  err = MatDestroy(Jp);dCHK(err);
  err = VecDestroy(r);dCHK(err);
  err = VecDestroy(x);dCHK(err);
  err = VecDestroy(proj.x);dCHK(err);
  err = VecDestroy(proj.y);dCHK(err);
  for (dInt i=0; i<ALEN(proj.ruleset); i++) {
    if (proj.ruleset[i]) {err = dRulesetDestroy(proj.ruleset[i]);dCHK(err);}
  }
  for (dInt i=0; i<ALEN(proj.regioniter); i++) {
    if (proj.regioniter[i]) {err = dRulesetIteratorDestroy(proj.regioniter[i]);dCHK(err);}
  }
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  /* const char pTagName[] = "dohp_partition"; */
  iMesh_Instance mi;
  dJacobi jac;
  dFS fs;
  dMesh mesh;
  dMeshESH domain;
  dMeshTag dtag;
  MPI_Comm comm;
  PetscViewer viewer;
  dBool  showconn,showmesh,flg;
  dInt exactChoice,nset;
  dErr err;

  err = dInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = PetscOptionsBegin(comm,NULL,"Test options","ex1");dCHK(err); {
    gopt.constBDeg = 4; gopt.showsoln = dFALSE; gopt.cycles = 1; gopt.proj_version = 1;
    gopt.proj_qmethod = dQUADRATURE_METHOD_FAST;
    gopt.jac_qmethod = dQUADRATURE_METHOD_SPARSE; gopt.q1scale = 1.0;
    gopt.frequency[0] = 1; gopt.frequency[1] = 1; gopt.frequency[2] = 1;
    exactChoice = 0; showconn = dFALSE; showmesh = dFALSE;
    err = PetscOptionsInt("-const_bdeg","Use constant isotropic degree on all elements",NULL,gopt.constBDeg,&gopt.constBDeg,NULL);dCHK(err);
    err = PetscOptionsBool("-show_soln","Show solution vector immediately after solving",NULL,gopt.showsoln,&gopt.showsoln,NULL);dCHK(err);
    err = PetscOptionsInt("-exact","Exact solution choice (0=transcendental,1=x coord)",NULL,exactChoice,&exactChoice,NULL);dCHK(err);
    err = PetscOptionsInt("-cycles","Number of times to solve the equation, useful for profiling",NULL,gopt.cycles,&gopt.cycles,NULL);dCHK(err);
    err = PetscOptionsInt("-proj_version","Residual evaluation version",NULL,gopt.proj_version,&gopt.proj_version,NULL);dCHK(err);
    err = PetscOptionsEnum("-proj_qmethod","Quadrature method for residual evaluation",NULL,dQuadratureMethods,(PetscEnum)gopt.proj_qmethod,(PetscEnum*)&gopt.proj_qmethod,NULL);dCHK(err);
    err = PetscOptionsEnum("-jac_qmethod","Quadrature method for Jacobian evaluation",NULL,dQuadratureMethods,(PetscEnum)gopt.jac_qmethod,(PetscEnum*)&gopt.jac_qmethod,NULL);dCHK(err);
    err = PetscOptionsReal("-q1scale","Scale matrix entries of Q1 preconditioning matrix",NULL,gopt.q1scale,&gopt.q1scale,NULL);dCHK(err);
    nset = 3;
    err = PetscOptionsRealArray("-frequency","Frequency of oscillation in each cartesion direction",NULL,gopt.frequency,&nset,&flg);dCHK(err);
    if (flg && nset > 3) dERROR(PETSC_COMM_SELF,1,"frequency may be at most 3 values");
    err = PetscOptionsBool("-show_conn","Show connectivity",NULL,showconn,&showconn,NULL);dCHK(err);
    err = PetscOptionsBool("-show_mesh","Show mesh immediately after createHexMesh()",NULL,showmesh,&showmesh,NULL);dCHK(err);
    nset = 3;
    err = PetscOptionsRealArray("-require_ptwise","<L^1,L^2,L^infty> Error if pointwise norms exceed given values, negative to disable",NULL,gopt.normRequirePtwise,&nset,&flg);dCHK(err);
    if (flg && nset != 3) dERROR(PETSC_COMM_SELF,1,"You must set 3 values for -require_ptwise, %d set",nset);
    nset = 3;
    err = PetscOptionsRealArray("-require_grad","<L^1,L^2,L^infty> Error if pointwise gradient norms exceed given values, negative to disable",NULL,gopt.normRequireGrad,&nset,&flg);dCHK(err);
    if (flg && nset != 3) dERROR(PETSC_COMM_SELF,1,"You must set 3 values for -require_grad");
  } err = PetscOptionsEnd();dCHK(err);
  err = dMeshCreate(comm,&mesh);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  err = createHexMesh(mi);dCHK(err);
  if (showmesh) {err = dMeshView(mesh,viewer);dCHK(err);}
  err = verifyAdjacencies(mesh);dCHK(err);

  /* The root set cannot be tagged (at least in MOAB) */
  err = dMeshGetRoot(mesh,&domain);dCHK(err);
  err = dMeshSetDuplicateEntsOnly(mesh,domain,&domain);dCHK(err);

  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);

  err = tagHexes(mesh,&dtag);dCHK(err);
  if (showconn) {err = examine(mesh,dtag);dCHK(err);}

  err = dFSCreate(comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,domain);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  err = dFSSetFromOptions(fs);dCHK(err);

  err = useFS(fs);dCHK(err);

  switch (exactChoice) {
    case 0:
      exact.function = exact_0_function;
      exact.gradient = exact_0_gradient;
      break;
    case 1:
      exact.function = exact_1_function;
      exact.gradient = exact_1_gradient;
      break;
    default: dERROR(PETSC_COMM_SELF,1,"Exact solution choice %d invalid",exactChoice);
  }
  err = doProjection(fs);dCHK(err);

  err = dFSDestroy(fs);dCHK(err);
  err = dJacobiDestroy(jac);dCHK(err);
  err = dMeshDestroy(mesh);
  err = dFinalize();dCHK(err);
  return 0;
}
