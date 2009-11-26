static const char help[] = "Test viewer\n";

#include <dohp.h>
#include <dohpfs.h>
#include <dohpviewer.h>

#define ALEN(a) (dInt)(sizeof(a)/sizeof((a)[0]))

static dErr FSEx4CreateMesh(MPI_Comm comm,dMesh *inmesh)
{
  static const double vtx[14*3] = {0,0,0, 1,0,0, 0.5,1.5,0, -0.3,0.8,0, -1.5,0,0, -0.3,-0.8,0, 0.5,-1.5,0,
                                   0,0,1, 1,0,1, 0.5,1.5,1, -0.3,0.8,1, -1.5,0,1, -0.3,-0.8,1, 0.5,-1.5,1};
  static const dIInt
    econn[2*(7+2*9)] = {0,7, 1,8, 2,9, 3,10, 4,11, 5,12, 6,13,
                        0,1, 1,2, 2,3, 0,3, 3,4, 4,5, 0,5, 5,6, 6,1,
                        7,8, 8,9, 9,10, 7,10, 10,11, 11,12, 7,12, 12,13, 13,8}
  ,fconn[4*(2*3+9)] = {0,1,2,3, 0,3,4,5, 0,5,6,1,
                       7,8,9,10, 7,10,11,12, 7,12,13,8,
                       0,1,8,7, 1,2,9,8, 2,3,10,9,  0,3,10,7, 3,4,11,10, 4,5,12,11,  0,5,12,7, 5,6,13,12, 6,1,8,13}
  ,rconn[24] = {0,1,2,3,7,8,9,10, 0,3,4,5,7,10,11,12, 0,5,6,1,7,12,13,8};

  MeshListEH verts=MLZ,edges=MLZ,faces=MLZ,regions=MLZ;
  MeshListInt stat=MLZ;
  iBase_EntityHandle work[ALEN(econn)+ALEN(fconn)+ALEN(rconn)];
  dIInt ierr;
  dMesh mesh;
  iMesh_Instance mi;
  dErr err;

  dFunctionBegin;
  err = dMeshCreate(comm,&mesh);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  iMesh_createVtxArr(mi,ALEN(vtx)/3,iBase_INTERLEAVED,vtx,ALEN(vtx),MLREF(verts),&ierr);dICHK(mi,ierr);
  for (dInt i=0; i<ALEN(econn); i++) work[i] = verts.v[econn[i]];
  iMesh_createEntArr(mi,iMesh_LINE_SEGMENT,work,ALEN(econn),MLREF(edges),MLREF(stat),&ierr);dICHK(mi,ierr);
  MeshListFree(edges); MeshListFree(stat);
  for (dInt i=0; i<ALEN(fconn); i++) work[i] = verts.v[fconn[i]];
  iMesh_createEntArr(mi,iMesh_QUADRILATERAL,work,ALEN(fconn),MLREF(faces),MLREF(stat),&ierr);dICHK(mi,ierr);
  MeshListFree(faces); MeshListFree(stat);
  for (dInt i=0; i<ALEN(rconn); i++) work[i] = verts.v[rconn[i]];
  iMesh_createEntArr(mi,iMesh_HEXAHEDRON,work,ALEN(rconn),MLREF(regions),MLREF(stat),&ierr);dICHK(mi,ierr);
  MeshListFree(regions); MeshListFree(stat);
  *inmesh = mesh;
  dFunctionReturn(0);
}

typedef struct _n_FSEx4 {
  dInt rdeg,bdeg;
} *FSEx4;

static dErr FSEx4CheckSubMesh(FSEx4 ex4,dFS fs)
{
  dEntTopology *topo;
  PetscViewer  viewer;
  dInt         p,nelems,nverts,nconn;
  dInt         *off,*conn;
  dScalar      *coords;
  Vec          X;
  dErr         err;

  dFunctionBegin;
  err = dFSGetSubElementMeshSize(fs,&nelems,&nverts,&nconn);dCHK(err);
  p = ex4->bdeg;
  dASSERT(nelems == 3*(p-1)*(p-1)*(p-1));
  dASSERT(nverts == 3*p*p*p - 3*p*p + p);
  dASSERT(nconn  == 8*nelems);
  err = dMallocA3(nelems,&topo,nelems+1,&off,nconn,&conn);dCHK(err);
  err = dFSGetSubElementMesh(fs,nelems,nconn,topo,off,conn);dCHK(err);

  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = PetscViewerASCIIPrintf(viewer,"Offsets for %d subelements\n",nelems);dCHK(err);
  err = PetscIntView(nelems+1,off,viewer);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Connectivity\n");dCHK(err);
  err = PetscIntView(nconn,conn,viewer);dCHK(err);
  err = dFree3(topo,off,conn);dCHK(err);

  err = PetscViewerASCIIPrintf(viewer,"Coordinates\n");dCHK(err);
  err = dFSGetCoordinates(fs,&X);dCHK(err);
  err = VecView(X,viewer);dCHK(err);
  err = VecGetArray(X,&coords);dCHK(err);
  err = VecRestoreArray(X,&coords);dCHK(err);
  err = VecDestroy(X);dCHK(err);
  dFunctionReturn(0);
}


int main(int argc,char *argv[])
{
  MPI_Comm comm;
  FSEx4    ex4;
  dFS      fs;
  dErr     err;
  dMesh    mesh;
  dViewer  viewer;
  dJacobi  jac;
  dMeshTag rtag,dtag;
  dTruth   read;
  Vec      X;

  err = dInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  err = dCalloc(sizeof(struct _n_FSEx4),&ex4);dCHK(err);
  ex4->rdeg = 4;
  ex4->bdeg = 4;
  read      = dFALSE;
  err = PetscOptionsBegin(comm,NULL,"FS-Ex4 Options","");dCHK(err);
  {
    err = PetscOptionsInt("-bdeg","Number of nodes per element in each Cartesian direction","",ex4->bdeg,&ex4->bdeg,NULL);dCHK(err);
    err = PetscOptionsTruth("-read_back","Read the mesh back in","",read,&read,NULL);dCHK(err);
  }
  err = PetscOptionsEnd();dCHK(err);
  err = FSEx4CreateMesh(comm,&mesh);dCHK(err);
  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetDegrees(jac,9,2);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);

  err = dMeshCreateRuleTagIsotropic(mesh,0,jac,"fsex4_rule_degree",ex4->rdeg,&rtag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,0,jac,"fsex4_efs_degree",ex4->bdeg,&dtag);dCHK(err);
  err = dFSCreate(comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,0);dCHK(err);
  err = dFSSetRuleTag(fs,jac,rtag);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  err = dFSSetFromOptions(fs);dCHK(err);
  err = dFSCreateGlobalVector(fs,&X);dCHK(err);
  err = PetscObjectSetName((PetscObject)X,"my_vec");dCHK(err);

  err = PetscViewerCreate(comm,&viewer);dCHK(err);
  err = PetscViewerSetType(viewer,PETSC_VIEWER_DHM);dCHK(err);
  err = PetscViewerFileSetName(viewer,"fsex4.dhm");dCHK(err);
  err = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);dCHK(err);
  err = dViewerDHMSetTimeUnits(viewer,"hour",PETSC_PI*1e7/3600);dCHK(err);
  err = dViewerDHMSetTime(viewer,0.1);dCHK(err);
  err = VecView(X,viewer);dCHK(err);
  err = VecDestroy(X);dCHK(err);
  err = PetscViewerDestroy(viewer);dCHK(err);

  err = FSEx4CheckSubMesh(ex4,fs);dCHK(err);

  err = dJacobiDestroy(jac);dCHK(err);
  err = dFSDestroy(fs);dCHK(err);
  err = dMeshDestroy(mesh);dCHK(err);

  if (read) {
    PetscMPIInt rank;
    dInt nsteps;
    dReal *steptimes;
    err = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);dCHK(err);
    err = PetscViewerCreate(PETSC_COMM_SELF,&viewer);dCHK(err);
    err = PetscViewerSetType(viewer,PETSC_VIEWER_DHM);dCHK(err);
    err = PetscViewerFileSetName(viewer,"fsex4.dhm");dCHK(err);
    err = PetscViewerFileSetMode(viewer,FILE_MODE_READ);dCHK(err);
    err = dViewerDHMGetSteps(viewer,&nsteps,&steptimes);dCHK(err);
    err = dPrintf(PETSC_COMM_SELF,"[%d] DHM has %d steps available\n",rank,nsteps);dCHK(err);
    for (dInt i=0; i<nsteps; i++) {
      err = dPrintf(PETSC_COMM_SELF,"[%d] step %d  time %g\n",rank,i,steptimes[i]);dCHK(err);
    }
    err = dViewerDHMSetTimeStep(viewer,0);dCHK(err);
    err = dFSCreate(PETSC_COMM_SELF,&fs);dCHK(err);
    err = dFSSetType(fs,dFSCONT);dCHK(err);
    err = dFSLoadIntoFS(viewer,"my_vec",fs);dCHK(err);

    err = FSEx4CheckSubMesh(ex4,fs);dCHK(err);

    err = dViewerDHMRestoreSteps(viewer,&nsteps,&steptimes);dCHK(err);
    err = dFSDestroy(fs);dCHK(err);
    err = PetscViewerDestroy(viewer);dCHK(err);
  }

  err = dFree(ex4);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
