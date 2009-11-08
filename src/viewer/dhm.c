#include "dhm.h"

PetscErrorCode PetscViewerCreate_DHM(PetscViewer);

dErr dViewerDHMSetUp(PetscViewer viewer)
{
  dViewer_DHM *dhm = viewer->data;
  herr_t       herr;
  hid_t        plist_id,fid;

  dFunctionBegin;
  if (dhm->file) dFunctionReturn(0);
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
  dhm->file = fid;
  herr = H5Pclose(plist_id);dHCHK(herr);

  if (dhm->btype == FILE_MODE_READ) dFunctionReturn(0);

  /* set up the groups */
  dhm->dohproot = H5Gcreate(dhm->file,"dohp",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);     if (dhm->dohproot < 0) dERROR(PETSC_ERR_ORDER,"H5Gcreate");
  dhm->meshroot = H5Gcreate(dhm->dohproot,"mesh",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); if (dhm->meshroot < 0) dERROR(PETSC_ERR_ORDER,"H5Gcreate");
  dhm->fsroot   = H5Gcreate(dhm->dohproot,"fs",  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); if (dhm->fsroot   < 0) dERROR(PETSC_ERR_ORDER,"H5Gcreate");
  dhm->steproot = H5Gcreate(dhm->dohproot,"step",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); if (dhm->steproot < 0) dERROR(PETSC_ERR_ORDER,"H5Gcreate");
  dFunctionReturn(0);
}

static dErr PetscViewerDestroy_DHM(PetscViewer v)
{
  dErr err;
  herr_t herr;
  dViewer_DHM *dhm = v->data;

  dFunctionBegin;
  err = dFree(dhm->filename);dCHK(err);
  if (dhm->meshroot) {herr = H5Gclose(dhm->meshroot); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");}
  if (dhm->fsroot)   {herr = H5Gclose(dhm->fsroot);   if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");}
  if (dhm->steproot) {herr = H5Gclose(dhm->steproot); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");}
  if (dhm->dohproot) {herr = H5Gclose(dhm->dohproot); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");}
  if (dhm->file)     {herr = H5Fclose(dhm->file);     if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Fclose");}
  if (dhm->fs) {err = dFSDestroy(dhm->fs);dCHK(err);}
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
