#include "dhm.h"

extern PetscErrorCode PetscViewerCreate_DHM(PetscViewer);

dErr dViewerDHMGetStringTypes(PetscViewer viewer,hid_t *fstring,hid_t *mstring,hid_t *sspace)
{
  dViewer_DHM *dhm = viewer->data;
  herr_t       herr;

  dFunctionBegin;
  if (dhm->h5t_fstring < 0) {
    dhm->h5t_fstring = H5Tcopy(H5T_FORTRAN_S1);dH5CHK(dhm->h5t_fstring,H5Tcopy);
    herr = H5Tset_size(dhm->h5t_fstring,H5T_VARIABLE);dH5CHK(herr,H5Tset_size);
    herr = H5Tcommit(dhm->typeroot,"string",dhm->h5t_fstring,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(herr,H5Tcommit);
  }
  if (dhm->h5t_mstring < 0) {
    dhm->h5t_mstring = H5Tcopy(H5T_C_S1);dH5CHK(dhm->h5t_mstring,H5Tcopy);
    herr = H5Tset_size(dhm->h5t_mstring,H5T_VARIABLE);dH5CHK(herr,H5Tset_size);
  }
  if (fstring) *fstring = dhm->h5t_fstring;
  if (mstring) *mstring = dhm->h5t_mstring;
  if (sspace)  *sspace  = dhm->h5s_scalar;
  dFunctionReturn(0);
}

dErr dViewerDHMAttributeStringWrite(PetscViewer viewer,hid_t grp,const char *attname,const char *str)
{
  hid_t fstring,mstring,scalar,att;
  herr_t herr;
  dErr err;

  dFunctionBegin;
  err = dViewerDHMGetStringTypes(viewer,&fstring,&mstring,&scalar);dCHK(err);
  att = H5Acreate(grp,attname,fstring,scalar,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(att,H5Acreate);
  herr = H5Awrite(att,mstring,&str);dH5CHK(herr,H5Awrite);
  herr = H5Aclose(att);dH5CHK(herr,H5Aclose);
  dFunctionReturn(0);
}

dErr dViewerDHMWriteDimensions(PetscViewer viewer,hid_t grp,const char *name,const char *units,dReal scale)
{
  herr_t herr;
  hid_t dimgrp,satt,sspace;
  dErr err;

  dFunctionBegin;
  dimgrp = H5Gcreate(grp,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);                   if (dimgrp < 0) dERROR(PETSC_ERR_LIB,"H5Gcreate");
  sspace = H5Screate(H5S_SCALAR); if (sspace < 0) dERROR(PETSC_ERR_LIB,"H5Screate");  if (sspace < 0) dERROR(PETSC_ERR_LIB,"H5Screate");
  satt = H5Acreate(dimgrp,"scale",dH5T_REAL,sspace,H5P_DEFAULT,H5P_DEFAULT);if (satt < 0) dERROR(PETSC_ERR_LIB,"H5Acreate");
  herr = H5Awrite(satt,dH5T_REAL,&scale);  if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Awrite");
  err = dViewerDHMAttributeStringWrite(viewer,dimgrp,"dimensions",units);dCHK(err);
  herr = H5Aclose(satt);   if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Aclose");
  herr = H5Sclose(sspace); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Sclose");
  herr = H5Gclose(dimgrp); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");
  dFunctionReturn(0);
}

dErr dViewerDHMSetUp(PetscViewer viewer)
{
  dViewer_DHM *dhm = viewer->data;
  herr_t       herr;
  hid_t        plist_id,fid;
  dErr         err;

  dFunctionBegin;
  if (dhm->file > 0) dFunctionReturn(0);
  /* Set attributes for opening the file */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  herr = H5Pset_fapl_mpio(plist_id,((PetscObject)viewer)->comm,MPI_INFO_NULL);dH5CHK(herr,H5Pset_fapl_mpio);
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
  herr = H5Pclose(plist_id);dH5CHK(herr,H5Pclose);

  if (dhm->btype == FILE_MODE_READ) dFunctionReturn(0);

  /* set up the groups */
  dhm->dohproot = H5Gcreate(dhm->file,"dohp",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->dohproot,H5Gcreate);
  dhm->meshroot = H5Gcreate(dhm->dohproot,"mesh",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->meshroot,H5Gcreate);
  dhm->fsroot   = H5Gcreate(dhm->dohproot,"fs",  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->fsroot,H5Gcreate);
  dhm->steproot = H5Gcreate(dhm->dohproot,"step",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->steproot,H5Gcreate);
  dhm->typeroot = H5Gcreate(dhm->dohproot,"type",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->typeroot,H5Gcreate);

  dhm->h5s_scalar = H5Screate(H5S_SCALAR);dH5CHK(dhm->h5s_scalar,H5Screate);

  err = dViewerDHMWriteDimensions(viewer,dhm->steproot,"time",dhm->timeunits,dhm->timescale);dCHK(err);
  dFunctionReturn(0);
}

static dErr PetscViewerDestroy_DHM(PetscViewer v)
{
  dErr err;
  herr_t herr;
  dViewer_DHM *dhm = v->data;

  dFunctionBegin;
  err = dFree(dhm->filename);dCHK(err);
  if (dhm->h5t_fstring > 0) {herr = H5Tclose(dhm->h5t_fstring);dH5CHK(herr,H5Tclose);}
  if (dhm->h5t_mstring > 0) {herr = H5Tclose(dhm->h5t_mstring);dH5CHK(herr,H5Tclose);}
  if (dhm->h5s_scalar  > 0) {herr = H5Sclose(dhm->h5s_scalar);dH5CHK(herr,H5Sclose);}
  if (dhm->curstep     > 0) {herr = H5Gclose(dhm->curstep);dH5CHK(herr,H5Gclose);}
  if (dhm->meshroot    > 0) {herr = H5Gclose(dhm->meshroot);dH5CHK(herr,H5Gclose);}
  if (dhm->fsroot      > 0) {herr = H5Gclose(dhm->fsroot);dH5CHK(herr,H5Gclose);}
  if (dhm->steproot    > 0) {herr = H5Gclose(dhm->steproot);dH5CHK(herr,H5Gclose);}
  if (dhm->typeroot    > 0) {herr = H5Gclose(dhm->typeroot);dH5CHK(herr,H5Gclose);}
  if (dhm->dohproot    > 0) {herr = H5Gclose(dhm->dohproot);dH5CHK(herr,H5Gclose);}
  if (dhm->file        > 0) {herr = H5Fclose(dhm->file);dH5CHK(herr,H5Fclose);}
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
  if (time != dhm->time) {
    dhm->stepnumber++;
    if (dhm->curstep > 0) {
      hid_t herr;
      herr = H5Gclose(dhm->curstep);dH5CHK(herr,H5Gclose);
    }
  }
  dhm->time = time;
  dFunctionReturn(0);
}

static dErr dViewerDHMSetTimeUnits_DHM(PetscViewer v,const char *units,dReal scale)
{
  dViewer_DHM *dhm = v->data;
  dErr         err;

  dFunctionBegin;
  if (dhm->timeunits) {err = dFree(dhm->timeunits);}
  err = PetscStrallocpy(units,&dhm->timeunits);dCHK(err);
  dhm->timescale = scale;
  dFunctionReturn(0);
}


dErr dViewerDHMGetStep(PetscViewer viewer,hid_t *step)
{
  dViewer_DHM *dhm = viewer->data;
  dErr err;

  dFunctionBegin;
  if (dhm->curstep < 0) {
    char stepname[16];
    herr_t herr;
    hid_t att;
    err = PetscSNPrintf(stepname,sizeof stepname,"%03d",dhm->stepnumber);dCHK(err);
    dhm->curstep = H5Gcreate(dhm->steproot,stepname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->curstep,H5Gopen);
    att = H5Acreate(dhm->curstep,"time",dH5T_REAL,dhm->h5s_scalar,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(att,H5Acreate);
    herr = H5Awrite(att,dH5T_REAL,&dhm->time);dH5CHK(herr,H5Awrite);
    herr = H5Aclose(att);dH5CHK(herr,H5Aclose);
  }
  *step = dhm->curstep;
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
  dhm->btype      = (PetscFileMode) -1;
  dhm->filename   = 0;

  dhm->file        = -1;
  dhm->dohproot    = -1;
  dhm->typeroot    = -1;
  dhm->steproot    = -1;
  dhm->meshroot    = -1;
  dhm->fsroot      = -1;
  dhm->curstep     = -1;
  dhm->h5t_fstring = -1;
  dhm->h5t_mstring = -1;
  dhm->h5s_scalar  = -1;

  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"PetscViewerFileSetName_C","PetscViewerFileSetName_DHM",
                                           PetscViewerFileSetName_DHM);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"PetscViewerFileSetMode_C","PetscViewerFileSetMode_DHM",
                                           PetscViewerFileSetMode_DHM);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"dViewerDHMSetTime_C","dViewerDHMSetTime_DHM",dViewerDHMSetTime_DHM);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"dViewerDHMSetTimeUnits_C","dViewerDHMSetTimeUnits_DHM",dViewerDHMSetTimeUnits_DHM);dCHK(err);
  dFunctionReturn(0);
}
