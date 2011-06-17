static const char help[] = "Read a DHM file and dump its contents\n";

#include <dohpfs.h>
#include <dohpviewer.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>

int main(int argc, char *argv[])
{
  char filename[PETSC_MAX_PATH_LEN];
  dViewer load,view;
  dBool flg;
  dErr err;
  dInt nsteps;
  dReal *times;
  MPI_Comm comm;
  dBool view_submesh,view_fields;

  dInitialize(&argc,&argv,0,help);
  comm = PETSC_COMM_WORLD;
  view = PETSC_VIEWER_STDOUT_WORLD;
  err = PetscOptionsBegin(comm,NULL,"DHM Dump options",__FILE__);dCHK(err); {
    err = PetscOptionsString("-f","DHM file to read","",filename,filename,sizeof filename,&flg);dCHK(err);
    if (!flg) dERROR(comm,PETSC_ERR_USER,"Must specify a file to load with -f FILENAME");
    err = PetscOptionsBool("-submesh","View the sub-element mesh (big)","",view_submesh=dFALSE,&view_submesh,NULL);dCHK(err);
    err = PetscOptionsBool("-fields","View all the fields defined on the mesh (big)","",view_fields=dFALSE,&view_fields,NULL);dCHK(err);
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
      err = PetscViewerASCIIPrintf(view,"FS %D: %-20s nblocks=%D\n",i,fspaces[i].name,fspaces[i].nblocks);dCHK(err);
      err = PetscViewerASCIIPushTab(view);dCHK(err);
      err = dFSCreate(comm,&fs);dCHK(err);
      err = dFSSetType(fs,dFSCONT);dCHK(err);
      err = dFSSetOrderingType(fs,MATORDERINGNATURAL);dCHK(err);
      err = dFSLoadIntoFS(load,fspaces[i].name,fs);dCHK(err);
      err = dFSView(fs,view);dCHK(err);
      if (view_submesh) {err = dFSSubElementMeshView(fs,view);dCHK(err);}
      err = dFSDestroy(&fs);dCHK(err);
      err = PetscViewerASCIIPopTab(view);dCHK(err);
    }
    for (dInt i=0; i<nfields; i++) {
      const char *vecname = fields[i].name;
      Vec X;
      dFS fs;
      dReal norms[3];
      err = dFSCreate(comm,&fs);dCHK(err);
      err = dFSSetType(fs,dFSCONT);dCHK(err);
      err = dFSSetOrderingType(fs,MATORDERINGNATURAL);dCHK(err);
      err = dFSLoadIntoFS(load,vecname,fs);dCHK(err);
      err = dFSCreateGlobalVector(fs,&X);dCHK(err);
      err = VecDohpLoadIntoVector(load,vecname,X);dCHK(err);
      err = dFSDestroy(&fs);dCHK(err);
      err = VecNorm(X,NORM_1_AND_2,norms);dCHK(err);
      err = VecNorm(X,NORM_INFINITY,&norms[2]);dCHK(err);
      err = PetscViewerASCIIPrintf(view,"Vec %D: %-20s |x|_1 % 12.6e  |x|_2 % 12.6e  |x|_max % 12.6e\n",i,vecname,norms[0],norms[1],norms[2]);dCHK(err);
      if (view_fields) {err = VecView(X,view);dCHK(err);}
      err = VecDestroy(&X);dCHK(err);
    }
    err = PetscViewerASCIIPopTab(view);dCHK(err);
    err = dViewerDHMRestoreStepSummary(load,&nfspaces,&fspaces,&nfields,&fields);dCHK(err);
  }
  err = PetscViewerASCIIPopTab(view);dCHK(err);
  err = dViewerDHMRestoreSteps(load,&nsteps,&times);dCHK(err);
  err = PetscViewerDestroy(&load);dCHK(err);

  dFinalize();
  return 0;
}
