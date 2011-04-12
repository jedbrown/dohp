static const char help[] = "Read a DHM file and dump its contents\n";

#include <dohpfs.h>
#include <dohpviewer.h>
#include <dohp.h>
#include <dohpsys.h>

static dErr SubElementMeshView(dFS fs,dViewer view)
{
  dErr err;
  dInt nelem,nverts,nconn,*suboff,*subind,n,bs;
  dEntTopology *subtopo;
  Vec X;
  const dReal *x;

  dFunctionBegin;
  err = dFSGetSubElementMeshSize(fs,&nelem,&nverts,&nconn);dCHK(err);
  dASSERT(nconn == 8*nelem);
  err = dMallocA3(nelem,&subtopo,nelem+1,&suboff,nconn,&subind);dCHK(err);
  err = dFSGetSubElementMesh(fs,nelem,nconn,subtopo,suboff,subind);dCHK(err);
  err = PetscViewerASCIIPrintf(view,"SubElementMesh\n");dCHK(err);
  err = dIntTableView(nelem,8,subind,view,"subconn");dCHK(err);
  err = dFree3(subtopo,suboff,subind);dCHK(err);

  err = dFSGetNodalCoordinatesGlobal(fs,&X);dCHK(err);
  err = VecGetLocalSize(X,&n);dCHK(err);
  err = VecGetBlockSize(X,&bs);dCHK(err);
  err = VecGetArrayRead(X,&x);dCHK(err);
  err = dRealTableView(n/bs,bs,x,view,"coords");dCHK(err);
  err = VecRestoreArrayRead(X,&x);dCHK(err);
  dFunctionReturn(0);
}

int main(int argc, char *argv[])
{
  char filename[PETSC_MAX_PATH_LEN];
  dViewer load,view;
  dBool flg;
  dErr err;
  dInt nsteps;
  dReal *times;
  MPI_Comm comm;

  dInitialize(&argc,&argv,0,help);
  comm = PETSC_COMM_WORLD;
  view = PETSC_VIEWER_STDOUT_WORLD;
  err = PetscOptionsBegin(comm,NULL,"DHM Dump options",__FILE__);dCHK(err); {
    err = PetscOptionsString("-f","DHM file to read","",filename,filename,sizeof filename,&flg);dCHK(err);
    if (!flg) dERROR(comm,PETSC_ERR_USER,"Must specify a file to load with -f FILENAME");
  } err = PetscOptionsEnd();dCHK(err);
  err = PetscViewerCreate(comm,&load);dCHK(err);
  err = PetscViewerSetType(load,PETSCVIEWERDHM);dCHK(err);
  err = PetscViewerFileSetName(load,filename);dCHK(err);
  err = PetscViewerFileSetMode(load,FILE_MODE_READ);dCHK(err);

  err = dViewerDHMGetSteps(load,&nsteps,&times);dCHK(err);
  err = PetscViewerASCIIPrintf(view,"number of time steps: %D\n",nsteps);dCHK(err);
  err = PetscViewerASCIIPushTab(view);dCHK(err);
  for (dInt step=0; step<nsteps; step++) {
    const struct dViewerDHMSummaryFS *fspaces;
    const struct dViewerDHMSummaryField *fields;
    dInt nfspaces,nfields;
    err = dViewerDHMSetTimeStep(load,step);dCHK(err);
    err = dViewerDHMGetStepSummary(load,&nfspaces,&fspaces,&nfields,&fields);dCHK(err);
    err = PetscViewerASCIIPrintf(view,"Step %D (%G): %D fspaces, %D fields\n",step,times[step],nfspaces,nfields);dCHK(err);
    err = PetscViewerASCIIPushTab(view);dCHK(err);
    for (dInt i=0; i<nfspaces; i++) {
      dFS fs;
      err = PetscViewerASCIIPrintf(view,"FS %D: name=%s, nblocks=%D\n",i,fspaces[i].name,fspaces[i].nblocks);dCHK(err);
      err = PetscViewerASCIIPushTab(view);dCHK(err);
      err = dFSCreate(comm,&fs);dCHK(err);
      err = dFSSetType(fs,dFSCONT);dCHK(err);
      err = dFSSetOrderingType(fs,MATORDERINGNATURAL);dCHK(err);
      err = dFSLoadIntoFS(load,fspaces[i].name,fs);dCHK(err);
      err = dFSView(fs,view);dCHK(err);
      err = SubElementMeshView(fs,view);dCHK(err);
      err = dFSDestroy(fs);dCHK(err);
      err = PetscViewerASCIIPopTab(view);dCHK(err);
    }
    err = PetscViewerASCIIPopTab(view);dCHK(err);
    err = dViewerDHMRestoreStepSummary(load,&nfspaces,&fspaces,&nfields,&fields);dCHK(err);
  }
  err = PetscViewerASCIIPopTab(view);dCHK(err);
  err = dViewerDHMRestoreSteps(load,&nsteps,&times);dCHK(err);
  err = PetscViewerDestroy(load);dCHK(err);

  dFinalize();
  return 0;
}
