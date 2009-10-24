#include "dhm.h"
#include <../src/sys/viewer/viewerimpl.h>

PetscErrorCode PetscViewerCreate_DHM(PetscViewer);

static dErr DHMSetUp(PetscViewer viewer)
{
  dViewer_DHM *dhm = viewer->data;
  dErr         err;
  herr_t       herr;
  hid_t        plist_id,fid;

  dFunctionBegin;
  if (dhm->file_id) dFunctionReturn(0);
  /* Set attributes for opening the file */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  herr = H5Pset_fapl_mpio(plist_id,((PetscObject)viewer)->comm,MPI_INFO_NULL);dHCHK(herr);
  /* Create or open the file collectively */
  switch (dhm->btype) {
    case FILE_MODE_READ:
      fid = H5Fopen(dhm->filename,H5F_ACC_RDONLY,plist_id);
      if (fid < 0) dERROR(PETSC_ERR_LIB,"H5Fopen(\"%s\",H5F_ACC_RDONLY,...) failed",dhm->filename);
      break;
    case FILE_MODE_WRITE:
      fid = H5Fcreate(dhm->filename,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
      if (fid < 0) dERROR(PETSC_ERR_LIB,"H5Fcreate(\"%s\",H5F_ACC_TRUNC,...) failed",dhm->filename);
      break;
    case FILE_MODE_APPEND:
      fid = H5Fopen(dhm->filename,H5F_ACC_RDWR,plist_id);
      if (fid < 0) dERROR(PETSC_ERR_LIB,"H5Fopen(\"%s\",H5F_ACC_RDWR,...) failed",dhm->filename);
      break;
    default:
      dERROR(PETSC_ERR_ORDER,"Must call PetscViewerFileSetMode() before PetscViewerFileSetName()");
  }
  dhm->file_id = fid;
  herr = H5Pclose(plist_id);dHCHK(herr);

  /* \todo set up the groups */
  dFunctionReturn(0);
}

static dErr PetscViewerDestroy_DHM(PetscViewer v)
{
  dErr err;
  herr_t herr;
  dViewer_DHM *dhm = v->data;

  dFunctionBegin;
  err = dFree(dhm->filename);dCHK(err);
  if (dhm->file_id) {herr = H5Fclose(dhm->file_id);dHCHK(herr);}
  if (dhm->fs)      {err = dFSDestroy(dhm->fs);dCHK(err);}
  err = dFree(dhm);dCHK(err);
  dFunctionReturn(0);
}

static dErr PetscViewerFileSetName_DHM(PetscViewer v,const char *name)
{
  dViewer_DHM *dhm = v->data;
  dErr err;

  dFunctionBegin;
  err = dFree(dhm->filename);dCHK(err);
  err = PetscStrallocpy(name,&dhm->filename);dCHK(err);
  dFunctionReturn(0);
}

static dErr PetscViewerFileSetMode_DHM(PetscViewer v,PetscFileMode btype)
{
  dViewer_DHM *dhm = v->data;

  dFunctionBegin;
  dhm->btype = btype;
  dFunctionReturn(0);
}

static dErr dViewerDHMSetTime_DHM(PetscViewer v,dReal time)
{
  dViewer_DHM *dhm = v->data;

  dFunctionBegin;
  dhm->time = time;
  dFunctionReturn(0);
}

static dErr dViewerDHMSetFS_DHM(PetscViewer v,dFS fs)
{
  dViewer_DHM *dhm = v->data;
  dErr err;

  dFunctionBegin;
  err = PetscObjectReference((dObject)fs);dCHK(err);
  if (dhm->fs) {err = dFSDestroy(dhm->fs);dCHK(err);}
  dhm->fs = fs;
  dFunctionReturn(0);
}

PetscErrorCode PetscViewerCreate_DHM(PetscViewer v)
{
  dErr err;
  dViewer_DHM *dhm;

  dFunctionBegin;
  err = PetscNewLog(v,dViewer_DHM,&dhm);dCHK(err);
  v->data         = (void *) dhm;
  v->ops->destroy = PetscViewerDestroy_DHM;
  v->ops->flush   = 0;
  v->iformat      = 0;
  dhm->btype     = (PetscFileMode) -1;
  dhm->filename  = 0;

  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"PetscViewerFileSetName_C","PetscViewerFileSetName_DHM",
                                           PetscViewerFileSetName_DHM);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"PetscViewerFileSetMode_C","PetscViewerFileSetMode_DHM",
                                           PetscViewerFileSetMode_DHM);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"dViewerDHMSetFS_C","dViewerDHMSetFS_DHM",dViewerDHMSetFS_DHM);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"dViewerDHMSetTime_C","dViewerDHMSetTime_DHM",dViewerDHMSetTime_DHM);dCHK(err);
  dFunctionReturn(0);
}
