static const char help[] = "Unit test for boundary condition manipulation\n\n";

#include <dohpfs.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>

typedef struct BUnitCtx *BU;
struct BUnitCtx {
  MPI_Comm comm;
  Vec      x,y,gx;
  dJacobi  jac;
  dMesh    mesh;
  dFS      fs;
  dInt     constBDeg;
};

static dErr BUCreate(MPI_Comm comm,BU *bunit)
{
  dErr err;
  BU bu;

  dFunctionBegin;
  *bunit = 0;
  err = dNew(struct BUnitCtx,&bu);dCHK(err);

  bu->comm        = comm;
  bu->constBDeg   = 2;          /* With only one element, this is the minimum possible with Dirichlet boundary conditions */

  *bunit = bu;
  dFunctionReturn(0);
}

static dErr BUSetBoundaries(BU bu)
{
  dErr err;

  dFunctionBegin;
  err = dFSRegisterBoundary(bu->fs,100,dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  dFunctionReturn(0);
}

static dErr BUSetFromOptions(BU bu)
{
  dMesh mesh;
  dFS fs;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag dtag;
  dErr err;

  dFunctionBegin;
  /* no options yet */

  /* Create mesh */
  err = dMeshCreate(bu->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);dCHK(err);
  bu->mesh = mesh;
  err = dMeshGetRoot(mesh,&domain);dCHK(err);
  err = dMeshSetDuplicateEntsOnly(mesh,domain,&domain);dCHK(err);

  err = dJacobiCreate(bu->comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  bu->jac = jac;

  err = dMeshCreateRuleTagIsotropic(mesh,domain,"bu_efs_degree",bu->constBDeg,&dtag);dCHK(err);

  err = dFSCreate(bu->comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,domain);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  bu->fs = fs;

  err = BUSetBoundaries(bu);dCHK(err);
  err = dFSSetFromOptions(fs);dCHK(err);

  err = dFSCreateGlobalVector(fs,&bu->gx);dCHK(err);
  err = dFSCreateExpandedVector(fs,&bu->x);dCHK(err);
  err = VecDuplicate(bu->x,&bu->y);dCHK(err);
  dFunctionReturn(0);
}

static dErr BUDestroy(BU bu)
{
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(bu->fs);dCHK(err);
  err = dJacobiDestroy(bu->jac);dCHK(err);
  err = dMeshDestroy(bu->mesh);dCHK(err);
  err = VecDestroy(bu->x);dCHK(err);
  err = VecDestroy(bu->y);dCHK(err);
  err = VecDestroy(bu->gx);dCHK(err);
  err = dFree(bu);dCHK(err);
  dFunctionReturn(0);
}

static dErr BUView(BU bu,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = PetscViewerASCIIPrintf(viewer,"Boundary Unit test object (BU)\n");dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"constant basis degree %2d\n",bu->constBDeg);dCHK(err);
  err = dFSView(bu->fs,viewer);dCHK(err);
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr BUAssemble(BU bu,Mat P)
{
  dErr     err;
  Vec      gc,coords;
  dScalar *g;
  dViewer  viewer = PETSC_VIEWER_STDOUT_WORLD;
  dInt     low,hi;

  dFunctionBegin;
  err = VecDohpGetClosure(bu->gx,&gc);dCHK(err);
  err = VecSet(gc,10);dCHK(err);
  err = dFSInhomogeneousDirichletCommit(bu->fs,gc);dCHK(err);
  err = VecGetOwnershipRange(gc,&low,&hi);dCHK(err);
  err = VecGetArray(gc,&g);dCHK(err);
  for (dInt i=0; i<hi-low; i++) g[i] = 1000 + low + i;
  err = VecRestoreArray(gc,&g);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Global closure vector, should be a sequence starting with 1000\n");dCHK(err);
  err = VecView(gc,viewer);dCHK(err);
  err = VecDohpRestoreClosure(bu->gx,&gc);dCHK(err);
  err = dFSGlobalToExpanded(bu->fs,bu->gx,bu->x,dFS_HOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Expanded vector projected into homogeneous space\n");dCHK(err);
  err = VecView(bu->x,viewer);dCHK(err);
  err = dFSGlobalToExpanded(bu->fs,bu->gx,bu->x,dFS_INHOMOGENEOUS,INSERT_VALUES);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Expanded vector projected into inhomogeneous space\n");dCHK(err);
  err = VecView(bu->x,viewer);dCHK(err);
  err = dFSGetNodalCoordinatesGlobal(bu->fs,&coords);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Global coordinate vector (3 dofs per closure node)\n");dCHK(err);
  err = VecView(coords,viewer);dCHK(err);
  err = MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd  (P,MAT_FINAL_ASSEMBLY);dCHK(err);
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  BU          bu;
  MPI_Comm    comm;
  PetscViewer viewer;
  Mat         Jp;
  dErr        err;

  err = dInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;

  err = BUCreate(comm,&bu);dCHK(err);
  err = BUSetFromOptions(bu);dCHK(err);

  err = BUView(bu,viewer);dCHK(err);

  err = dFSGetMatrix(bu->fs,MATSEQAIJ,&Jp);dCHK(err);
  err = BUAssemble(bu,Jp);dCHK(err);
  err = MatView(Jp,viewer);dCHK(err);

  err = MatDestroy(Jp);dCHK(err);
  err = BUDestroy(bu);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
