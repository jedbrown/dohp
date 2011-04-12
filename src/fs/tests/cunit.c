static const char help[] = "Unit test for coordinate evaluation\n\n";

#include <dohpfs.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>
#include <petscpf.h>

#define ALEN(a) ((dInt)(sizeof(a)/sizeof((a)[0])))
typedef enum {DEFORM_IDENTITY,DEFORM_AFFINE} DeformType;
static const char *const DeformTypes[] = {"IDENTITY","AFFINE","DeformTypes","DEFORM_",0};

typedef struct CUnitCtx *CU;
struct CUnitCtx {
  MPI_Comm comm;
  Vec      x,gx;
  dJacobi  jac;
  dMesh    mesh;
  dFS      fs;
  dInt     constBDeg;
  DeformType deform_type;
  dBool    expanded_view,cexp_view,nce_view,nce_compare;
};

static dBool FuzzyEquals(dScalar a,dScalar b) {return (dBool)(dAbs(a-b)/(1.e-50+dAbs(a)+dAbs(b)) < 1e-12);}
static dBool FuzzyEquals3(const dScalar *a,const dScalar *b)
{return FuzzyEquals(a[0],b[0]) && FuzzyEquals(a[1],b[1]) && FuzzyEquals(a[2],b[2]);}

typedef struct {
  dScalar A[3][3];
  dScalar b[3];
} Deform_Affine;
static dErr DeformApply_Affine(void *ctx,dInt n,const dScalar *x,dScalar *y)
{
  Deform_Affine *aff = ctx;
  dInt i=0;

  dFunctionBegin;
  for (i=0; i<n; i++) {
    y[3*i+0] = aff->A[0][0]*x[3*i+0] + aff->A[0][1]*x[3*i+1] + aff->A[0][2]*x[3*i+2] + aff->b[0];
    y[3*i+1] = aff->A[1][0]*x[3*i+0] + aff->A[1][1]*x[3*i+1] + aff->A[1][2]*x[3*i+2] + aff->b[1];
    y[3*i+2] = aff->A[2][0]*x[3*i+0] + aff->A[2][1]*x[3*i+1] + aff->A[2][2]*x[3*i+2] + aff->b[2];
  }
  dFunctionReturn(0);
}

static dErr DeformFree1(void *p)
{
  dErr err;

  dFunctionBegin;
  err = dFree(p);dCHK(err);
  dFunctionReturn(0);
}

static dErr CUCreate(MPI_Comm comm,CU *cunit)
{
  dErr err;
  CU cu;

  dFunctionBegin;
  *cunit = 0;
  err = dNew(struct CUnitCtx,&cu);dCHK(err);

  cu->comm        = comm;
  cu->deform_type = DEFORM_IDENTITY;
  cu->constBDeg   = 2;          /* With only one element, this is the minimum possible with Dirichlet boundary conditions */

  *cunit = cu;
  dFunctionReturn(0);
}

/* Create a mesh with two hexagons */
static dErr CUCreateHexMesh(CU cu)
{
  static const int vtxlen = 12*3;
  static const double vtx_base[12*3] = {0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                        0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                        2,0,0, 2,1,0, 2,1,1, 2,0,1};
  double vtx_mapped[12*3];
  int rconn[16] = {0,1,2,3,4,5,6,7, 1,2,6,5,8,9,10,11};
  int fconn[11*4] = {0,1,2,3, 1,2,6,5, 2,3,7,6, 0,3,7,4, 0,1,5,4, 5,6,7,4,
                     8,11,5,1, 2,1,8,9, 8,9,10,11, 11,10,6,5, 9,10,6,2};
  int econn[20*2] = {0,1,4,5,7,6,3,2, 1,8,5,11,6,10,2,9,
                     0,3,3,7,7,4,4,0, 1,2,2,6,6,5,5,1, 8,11,11,10,10,9,9,8};
  iBase_EntityHandle work[100];
  iMesh_Instance mi;
  MeshListEH v=MLZ,e=MLZ,f=MLZ,r=MLZ,tv=MLZ;
  MeshListInt stat=MLZ,off=MLZ;
  dIInt ierr;
  PF pf;
  dErr err;

  dFunctionBegin;
  err = PFCreate(cu->comm,3,3,&pf);dCHK(err);
  switch (cu->deform_type) {
    case DEFORM_IDENTITY:
      err = PFSetType(pf,PFIDENTITY,0);dCHK(err);
      break;
    case DEFORM_AFFINE: {
      Deform_Affine *ctx,tmpctx = {.A = {{1,2,3},{4,5,6},{7,8,9}},.b = {10,11,12}};
      err = PetscMalloc(sizeof(Deform_Affine),&ctx);dCHK(err);
      *ctx = tmpctx;
      err = PFSet(pf,DeformApply_Affine,0,0,DeformFree1,ctx);dCHK(err);
    } break;
    default: dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Deformation number %D not recognized",(PetscInt)cu->deform_type);
  }
  err = PFApply(pf,vtxlen/3,vtx_base,vtx_mapped);dCHK(err);

  err = dMeshCreate(cu->comm,&cu->mesh);dCHK(err);
  err = dMeshSetType(cu->mesh,dMESHSERIAL);dCHK(err);
  err = dMeshGetInstance(cu->mesh,&mi);dCHK(err);
  iMesh_createVtxArr(mi,vtxlen/3,iBase_INTERLEAVED,vtx_mapped,vtxlen,MLREF(v),&ierr);dICHK(mi,ierr);
  for (dInt i=0; i<ALEN(rconn); i++) work[i] = v.v[rconn[i]];
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
  err = PFDestroy(pf);dCHK(err);
  dFunctionReturn(0);
}

static dErr CUSetFromOptions(CU cu)
{
  dFS fs;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag dtag;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Hex mesh options",NULL);dCHK(err);
  err = PetscOptionsInt("-const_Bdeg","Constant interpolation degree",NULL,cu->constBDeg,&cu->constBDeg,NULL);dCHK(err);
  err = PetscOptionsEnum("-deform_type","Deformation to apply",NULL,DeformTypes,cu->deform_type,(PetscEnum*)&cu->deform_type,NULL);dCHK(err);
  err = PetscOptionsBool("-expanded_view","Show expanded vector which shows connectivity",NULL,cu->expanded_view,&cu->expanded_view,NULL);dCHK(err);
  err = PetscOptionsBool("-cexp_view","Show expanded coordinates (in coordinate basis)",NULL,cu->cexp_view,&cu->cexp_view,NULL);dCHK(err);
  err = PetscOptionsBool("-nce_view","Show coordinates evaluated at nodes of solution basis",NULL,cu->nce_view,&cu->nce_view,NULL);dCHK(err);
  err = PetscOptionsBool("-nce_compare","Compare global nodal coordinates re-expanded (consistency check)",NULL,cu->nce_compare,&cu->nce_compare,NULL);dCHK(err);
  err = PetscOptionsEnd();dCHK(err);

  /* Create mesh */
  err = CUCreateHexMesh(cu);dCHK(err);
  err = dMeshGetRoot(cu->mesh,&domain);dCHK(err);
  err = dMeshSetDuplicateEntsOnly(cu->mesh,domain,&domain);dCHK(err);

  err = dJacobiCreate(cu->comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  cu->jac = jac;

  err = dMeshCreateRuleTagIsotropic(cu->mesh,domain,"cu_efs_degree",cu->constBDeg,&dtag);dCHK(err);

  err = dFSCreate(cu->comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,cu->mesh,domain);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  cu->fs = fs;

  err = dFSSetFromOptions(fs);dCHK(err);

  err = dFSCreateGlobalVector(fs,&cu->gx);dCHK(err);
  err = dFSCreateExpandedVector(fs,&cu->x);dCHK(err);
  dFunctionReturn(0);
}

static dErr CUDestroy(CU cu)
{
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(cu->fs);dCHK(err);
  err = dJacobiDestroy(cu->jac);dCHK(err);
  err = dMeshDestroy(cu->mesh);dCHK(err);
  err = VecDestroy(cu->x);dCHK(err);
  err = VecDestroy(cu->gx);dCHK(err);
  err = dFree(cu);dCHK(err);
  dFunctionReturn(0);
}

static dErr CUView(CU cu,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = PetscViewerASCIIPrintf(viewer,"Boundary Unit test object (CU)\n");dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"constant basis degree %2d\n",cu->constBDeg);dCHK(err);
  err = dFSView(cu->fs,viewer);dCHK(err);
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr VecBlockCompare(Vec X,Vec Y,dViewer viewer)
{
  dErr err;
  dInt i,m,n,bs;
  const dScalar *x,*y;

  dFunctionBegin;
  dValidHeader(X,VEC_CLASSID,1);
  dValidHeader(Y,VEC_CLASSID,2);
  dValidHeader(viewer,PETSC_VIEWER_CLASSID,3);
  err = VecGetLocalSize(X,&m);dCHK(err);
  err = VecGetLocalSize(Y,&n);dCHK(err);
  if (m != n) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Size of X %D not compatible with Y %D",m,n);
  err = VecGetBlockSize(X,&bs);dCHK(err);
  if (bs != 3) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"X block size %D, but only implemented for bs=3",bs);
  err = VecGetBlockSize(Y,&bs);dCHK(err);
  if (bs != 3) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"Y block size %D, but only implemented for bs=3",bs);

  err = VecGetArrayRead(X,&x);dCHK(err);
  err = VecGetArrayRead(Y,&y);dCHK(err);
  for (i=0; i<m; i+=bs) {
    if (!FuzzyEquals3(&x[i],&y[i])) {
      err = PetscViewerASCIISynchronizedPrintf(viewer,"Node %3D: X {%10G %10G %10G} != Y {%10G %10G %10G}\n",i/bs,PetscRealPart(x[i]),PetscRealPart(x[i+1]),PetscRealPart(x[i+2]),
                                               PetscRealPart(y[i]),PetscRealPart(y[i+1]),PetscRealPart(y[i+2]));dCHK(err);
    }
  }
  err = PetscViewerFlush(viewer);dCHK(err);
  err = VecRestoreArrayRead(X,&x);dCHK(err);
  err = VecRestoreArrayRead(Y,&y);dCHK(err);
  dFunctionReturn(0);
}

static dErr CUTestEvaluate(CU cu,dViewer viewer)
{
  dErr err;
  dFS cfs,fs3;
  Vec ncglobal,ncexpanded,expanded2;

  dFunctionBegin;
  err = dFSGetCoordinateFS(cu->fs,&cfs);dCHK(err);
  err = dFSGetNodalCoordinateFS(cu->fs,&fs3);dCHK(err);
  {
    err = dFSGetNodalCoordinatesGlobal(cu->fs,&ncglobal);dCHK(err);
    err = dFSGetNodalCoordinatesExpanded(cu->fs,&ncexpanded);dCHK(err);
    if (cu->nce_view) {
      err = PetscViewerASCIIPrintf(viewer,"Nodal coordinates in expanded space (checks evaluation; -nce_view)\n");dCHK(err);
      err = VecBlockView(ncexpanded,viewer);dCHK(err);
    }
    err = VecDuplicate(ncexpanded,&expanded2);dCHK(err);
    err = dFSGlobalToExpanded(fs3,ncglobal,expanded2,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err);
    if (cu->nce_compare) {
      err = PetscViewerASCIIPrintf(viewer,"Comparison of re-expanded global coordinates (-nce_compare)\n");dCHK(err);
      err = VecBlockCompare(ncexpanded,expanded2,viewer);dCHK(err);
    }
    err = VecDestroy(expanded2);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr CUTest(CU cu,dViewer viewer)
{
  dErr     err;
  Vec      gc,coords;
  dScalar *g;
  dInt     low,hi,N,bs;

  dFunctionBegin;
  /* Check connectivity */
  err = VecDohpGetClosure(cu->gx,&gc);dCHK(err);
  err = VecGetOwnershipRange(gc,&low,&hi);dCHK(err);
  err = VecGetArray(gc,&g);dCHK(err);
  for (dInt i=0; i<hi-low; i++) g[i] = 1000 + low + i;
  err = VecRestoreArray(gc,&g);dCHK(err);
  err = VecDohpRestoreClosure(cu->gx,&gc);dCHK(err);

  err = dFSGlobalToExpanded(cu->fs,cu->gx,cu->x,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err);
  if (cu->expanded_view) {
    err = PetscViewerASCIIPrintf(viewer,"Expanded vector (-expanded_view)\n");dCHK(err);
    err = VecView(cu->x,viewer);dCHK(err);
  }

  /* View coordinates in Q1 basis derived from the mesh */
  err = dFSGetGeometryVectorExpanded(cu->fs,&coords);dCHK(err);
  err = VecGetSize(coords,&N);dCHK(err);
  err = VecGetBlockSize(coords,&bs);dCHK(err);
  if (N%3 || bs !=3) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Vector size unexpected");
  if (cu->cexp_view) {
    err = PetscViewerASCIIPrintf(viewer,"Expanded coordinate vector (3 dofs per expanded node, %D nodes; -cexp_view)\n",N/3);dCHK(err);
    err = VecBlockView(coords,viewer);dCHK(err);
  }

  err = CUTestEvaluate(cu,viewer);dCHK(err);
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  CU          cu;
  MPI_Comm    comm;
  PetscViewer viewer;
  dErr        err;

  err = dInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = CUCreate(comm,&cu);dCHK(err);
  err = CUSetFromOptions(cu);dCHK(err);
  err = CUView(cu,viewer);dCHK(err);
  err = CUTest(cu,viewer);dCHK(err);
  err = CUDestroy(cu);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
