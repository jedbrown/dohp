#include "dhm.h"
#include <dohp.h>

extern PetscErrorCode PetscViewerCreate_DHM(PetscViewer);

static dErr dViewerDHM_H5Tcommit(PetscViewer viewer,const char *typename,hid_t type)
{
  dViewer_DHM *dhm = viewer->data;
  dErr        err;
  dTruth      match;
  herr_t      herr;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_DHM,&match);dCHK(err);
  if (!match) dERROR(PETSC_ERR_ARG_WRONG,"Viewer type must be DHM");
  if (dhm->btype == FILE_MODE_WRITE) {
    herr = H5Tcommit(dhm->typeroot,typename,type,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(herr,H5Tinsert);
  }
  dFunctionReturn(0);
}

dErr dViewerDHMGetStringTypes(PetscViewer viewer,hid_t *fstring,hid_t *mstring,hid_t *sspace)
{
  dViewer_DHM *dhm = viewer->data;
  herr_t       herr;
  dErr         err;

  dFunctionBegin;
  if (dhm->h5t_fstring < 0) {
    dhm->h5t_fstring = H5Tcopy(H5T_C_S1);dH5CHK(dhm->h5t_fstring,H5Tcopy);
    herr = H5Tset_size(dhm->h5t_fstring,H5T_VARIABLE);dH5CHK(herr,H5Tset_size);
    err = dViewerDHM_H5Tcommit(viewer,"string",dhm->h5t_fstring);dCHK(err);
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

static dErr dViewerDHMGetUnitsType(PetscViewer viewer,hid_t *unitstype)
{
  dViewer_DHM *dhm = viewer->data;
  dErr         err;

  dFunctionBegin;
  if (dhm->h5t_units < 0) {
    hid_t  type,strtype;
    herr_t herr;
    err  = dViewerDHMGetStringTypes(viewer,&strtype,NULL,NULL);dCHK(err);
    type = H5Tcreate(H5T_COMPOUND,sizeof(dht_Units));dH5CHK(type,H5Tcreate);
    herr = H5Tinsert(type,"dimensions",offsetof(dht_Units,dimensions),strtype);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(type,"scale",offsetof(dht_Units,scale),dH5T_SCALAR);dH5CHK(herr,H5Tinsert);
    err = dViewerDHM_H5Tcommit(viewer,"Units_type",type);dCHK(err);
    dhm->h5t_units = type;
  }
  *unitstype = dhm->h5t_units;
  dFunctionReturn(0);
}

static dErr dViewerDHMGetTimeType(PetscViewer viewer,hid_t *timetype)
{
  dViewer_DHM *dhm = viewer->data;
  dErr         err;

  dFunctionBegin;
  if (dhm->h5t_time < 0) {
    hid_t  type,unitstype;
    herr_t herr;
    err = dViewerDHMGetUnitsType(viewer,&unitstype);dCHK(err);
    type = H5Tcreate(H5T_COMPOUND,sizeof(dht_RealWithUnits));dH5CHK(type,H5Tcreate);
    herr = H5Tinsert(type,"value",offsetof(dht_RealWithUnits,value),dH5T_REAL);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(type,"units",offsetof(dht_RealWithUnits,units),unitstype);dH5CHK(herr,H5Tinsert);
    err = dViewerDHM_H5Tcommit(viewer,"Time_type",type);dCHK(err);
    dhm->h5t_time = type;
  }
  *timetype = dhm->h5t_time;
  dFunctionReturn(0);
}

dErr dViewerDHMGetFSType(PetscViewer viewer,hid_t *intype)
{
  dViewer_DHM *dhm = viewer->data;
  dErr err;

  dFunctionBegin;
  if (dhm->h5t_fs < 0) {
    hid_t strtype,fstype,unittype,fieldtype,fieldstype;
    herr_t herr;

    err = dViewerDHMGetStringTypes(viewer,&strtype,NULL,NULL);dCHK(err);
    err = dViewerDHMGetUnitsType(viewer,&unittype);dCHK(err);
    fieldtype = H5Tcreate(H5T_COMPOUND,sizeof(dht_Field));dH5CHK(fieldtype,H5Tcreate);
    herr = H5Tinsert(fieldtype,"name",offsetof(dht_Field,name),strtype);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(fieldtype,"units",offsetof(dht_Field,units),unittype);dH5CHK(herr,H5Tinsert);

    err = dViewerDHM_H5Tcommit(viewer,"Field_type",fieldtype);dCHK(err); /* skip this? */

    fieldstype = H5Tvlen_create(fieldtype);dH5CHK(fieldstype,H5Tvlen_create);

    fstype = H5Tcreate(H5T_COMPOUND,sizeof(dht_FS));dH5CHK(fstype,H5Tcreate);
    herr = H5Tinsert(fstype,"degree",offsetof(dht_FS,degree),strtype);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(fstype,"global_offset",offsetof(dht_FS,global_offset),strtype);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(fstype,"partition",offsetof(dht_FS,partition),strtype);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(fstype,"ordered_subdomain",offsetof(dht_FS,ordered_subdomain),strtype);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(fstype,"mesh",offsetof(dht_FS,mesh),H5T_STD_REF_OBJ);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(fstype,"time",offsetof(dht_FS,time),dH5T_REAL);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(fstype,"internal_state",offsetof(dht_FS,internal_state),dH5T_INT);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(fstype,"number_of_subdomains",offsetof(dht_FS,number_of_subdomains),dH5T_INT);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(fstype,"fields",offsetof(dht_FS,fields),fieldstype);dH5CHK(herr,H5Tinsert);

    err = dViewerDHM_H5Tcommit(viewer,"FS_type",fstype);dCHK(err);

    herr = H5Tclose(fieldstype);dH5CHK(herr,H5Tclose);
    herr = H5Tclose(fieldtype);dH5CHK(herr,H5Tclose);
    dhm->h5t_fs = fstype;
  }
  *intype = dhm->h5t_fs;
  dFunctionReturn(0);
}

dErr dViewerDHMGetVecType(PetscViewer viewer,hid_t *intype)
{
  dViewer_DHM *dhm = viewer->data;
  dErr        err;

  dFunctionBegin;
  if (dhm->h5t_vec < 0) {
    hid_t vectype;
    herr_t herr;
    vectype = H5Tcreate(H5T_COMPOUND,sizeof(dht_Vec));dH5CHK(vectype,H5Tcreate);
    herr = H5Tinsert(vectype,"fs",offsetof(dht_Vec,fs),H5T_STD_REF_DSETREG);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(vectype,"time",offsetof(dht_Vec,time),dH5T_REAL);dH5CHK(herr,H5Tinsert);
    herr = H5Tinsert(vectype,"internal_state",offsetof(dht_Vec,internal_state),dH5T_INT);dH5CHK(herr,H5Tinsert);
    err = dViewerDHM_H5Tcommit(viewer,"Vec_type",vectype);dCHK(err);
    dhm->h5t_vec = vectype;
  }
  *intype = dhm->h5t_vec;
  dFunctionReturn(0);
}

dErr dViewerDHMSetUp(PetscViewer viewer)
{
  dViewer_DHM *dhm = viewer->data;
  herr_t       herr;
  hid_t        plist_id,fid;

  dFunctionBegin;
  if (dhm->file >= 0) dFunctionReturn(0);
  /* Set attributes for opening the file */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
#if defined dUSE_PARALLEL_HDF5  /* Writing variable-length date unsupported in HDF-1.8.3 */
  herr = H5Pset_fapl_mpio(plist_id,((PetscObject)viewer)->comm,MPI_INFO_NULL);dH5CHK(herr,H5Pset_fapl_mpio);
#endif
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

  /* set up the groups */
  if (dhm->btype == FILE_MODE_READ) {
    dhm->dohproot = H5Gopen(dhm->file,"dohp",H5P_DEFAULT);dH5CHK(dhm->dohproot,H5Gopen);
    dhm->meshroot = H5Gopen(dhm->dohproot,"mesh",H5P_DEFAULT);dH5CHK(dhm->meshroot,H5Gopen);
    dhm->fsroot   = H5Gopen(dhm->dohproot,"fs",H5P_DEFAULT);dH5CHK(dhm->fsroot,H5Gopen);
    dhm->steproot = H5Gopen(dhm->dohproot,"step",H5P_DEFAULT);dH5CHK(dhm->steproot,H5Gopen);
    dhm->typeroot = H5Gopen(dhm->dohproot,"type",H5P_DEFAULT);dH5CHK(dhm->typeroot,H5Gopen);
  } else {
    dhm->dohproot = H5Gcreate(dhm->file,"dohp",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->dohproot,H5Gcreate);
    dhm->meshroot = H5Gcreate(dhm->dohproot,"mesh",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->meshroot,H5Gcreate);
    dhm->fsroot   = H5Gcreate(dhm->dohproot,"fs",  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->fsroot,H5Gcreate);
    dhm->steproot = H5Gcreate(dhm->dohproot,"step",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->steproot,H5Gcreate);
    dhm->typeroot = H5Gcreate(dhm->dohproot,"type",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->typeroot,H5Gcreate);
  }

  dhm->h5s_scalar = H5Screate(H5S_SCALAR);dH5CHK(dhm->h5s_scalar,H5Screate);
  dFunctionReturn(0);
}

static dErr PetscViewerDestroy_DHM(PetscViewer v)
{
  dErr err;
  herr_t herr;
  dViewer_DHM *dhm = v->data;

  dFunctionBegin;
  err = dFree(dhm->filename);dCHK(err);
  err = dFree(dhm->timeunits);dCHK(err);
  if (dhm->h5t_fstring >= 0) {herr = H5Tclose(dhm->h5t_fstring);dH5CHK(herr,H5Tclose);}
  if (dhm->h5t_mstring >= 0) {herr = H5Tclose(dhm->h5t_mstring);dH5CHK(herr,H5Tclose);}
  if (dhm->h5t_fs      >= 0) {herr = H5Tclose(dhm->h5t_fs);dH5CHK(herr,H5Tclose);}
  if (dhm->h5t_vec     >= 0) {herr = H5Tclose(dhm->h5t_vec);dH5CHK(herr,H5Tclose);}
  if (dhm->h5t_units   >= 0) {herr = H5Tclose(dhm->h5t_units);dH5CHK(herr,H5Tclose);}
  if (dhm->h5t_time    >= 0) {herr = H5Tclose(dhm->h5t_time);dH5CHK(herr,H5Tclose);}
  if (dhm->h5s_scalar  >= 0) {herr = H5Sclose(dhm->h5s_scalar);dH5CHK(herr,H5Sclose);}
  if (dhm->curstep     >= 0) {herr = H5Gclose(dhm->curstep);dH5CHK(herr,H5Gclose);}
  if (dhm->meshroot    >= 0) {herr = H5Gclose(dhm->meshroot);dH5CHK(herr,H5Gclose);}
  if (dhm->fsroot      >= 0) {herr = H5Gclose(dhm->fsroot);dH5CHK(herr,H5Gclose);}
  if (dhm->steproot    >= 0) {herr = H5Gclose(dhm->steproot);dH5CHK(herr,H5Gclose);}
  if (dhm->typeroot    >= 0) {herr = H5Gclose(dhm->typeroot);dH5CHK(herr,H5Gclose);}
  if (dhm->dohproot    >= 0) {herr = H5Gclose(dhm->dohproot);dH5CHK(herr,H5Gclose);}
  if (dhm->file        >= 0) {herr = H5Fclose(dhm->file);dH5CHK(herr,H5Fclose);}
  err = dFree(dhm);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"PetscViewerFileSetName_C","",NULL);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"PetscViewerFileSetMode_C","",NULL);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"dViewerDHMSetTime_C","",NULL);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"dViewerDHMSetTimeUnits_C","",NULL);dCHK(err);
  dFunctionReturn(0);
}

static dErr dViewerDHMInvalidateCurrentStep(PetscViewer viewer)
{
  dViewer_DHM *dhm = viewer->data;

  dFunctionBegin;
  if (dhm->curstep > 0) {
    hid_t herr;
    herr = H5Gclose(dhm->curstep);dH5CHK(herr,H5Gclose);
  }
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
  dErr        err;

  dFunctionBegin;
  if (time != dhm->time) {
    dhm->stepnumber++;
    dhm->totalsteps++;
    err = dViewerDHMInvalidateCurrentStep(v);dCHK(err);
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
    hid_t att,timetype;
    dht_RealWithUnits time;

    err = dViewerDHMGetTimeType(viewer,&timetype);dCHK(err);
    err = PetscSNPrintf(stepname,sizeof stepname,"%03d",dhm->stepnumber);dCHK(err);
    if (dhm->btype == FILE_MODE_READ) {
      dhm->curstep = H5Gopen(dhm->steproot,stepname,H5P_DEFAULT);dH5CHK(dhm->curstep,H5Gopen);
      att          = H5Aopen(dhm->curstep,"time",H5P_DEFAULT);dH5CHK(att,H5Aopen);
      herr         = H5Aread(att,timetype,&time);dCHK(err);
      dhm->time      = time.value;
      err = dViewerDHMSetTimeUnits(viewer,time.units.dimensions,time.units.scale);dCHK(err);
    } else {
      dhm->curstep = H5Gcreate(dhm->steproot,stepname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dhm->curstep,H5Gcreate);
      att = H5Acreate(dhm->curstep,"time",timetype,dhm->h5s_scalar,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(att,H5Acreate);
      time.value            = dhm->time;
      time.units.dimensions = dhm->timeunits;
      time.units.scale      = dhm->timescale;
      herr = H5Awrite(att,timetype,&time);dH5CHK(herr,H5Awrite);
    }
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
  dhm->h5t_fs      = -1;
  dhm->h5t_vec     = -1;
  dhm->h5t_units   = -1;
  dhm->h5t_time    = -1;
  dhm->h5s_scalar  = -1;

  dhm->stepnumber  = -1;

  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"PetscViewerFileSetName_C","PetscViewerFileSetName_DHM",PetscViewerFileSetName_DHM);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"PetscViewerFileSetMode_C","PetscViewerFileSetMode_DHM",PetscViewerFileSetMode_DHM);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"dViewerDHMSetTime_C","dViewerDHMSetTime_DHM",dViewerDHMSetTime_DHM);dCHK(err);
  err = PetscObjectComposeFunctionDynamic((PetscObject)v,"dViewerDHMSetTimeUnits_C","dViewerDHMSetTimeUnits_DHM",dViewerDHMSetTimeUnits_DHM);dCHK(err);
  dFunctionReturn(0);
}



/*
* These really need to fail if the action has not been done, thus we don't gain anything by dynamic calls.
*/
dErr dViewerDHMSetTimeStep(PetscViewer viewer,dInt step)
{
  dViewer_DHM *dhm = viewer->data;
  dErr        err;
  dTruth      match;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_DHM,&match);dCHK(err);
  if (!match) dERROR(PETSC_ERR_ARG_WRONG,"Viewer type must be DHM");
  if (step < 0 || dhm->totalsteps <= step) dERROR(PETSC_ERR_ARG_OUTOFRANGE,"Step %d out of range [0 .. %d]",step,dhm->totalsteps);
  if (step != dhm->stepnumber) {
    err = dViewerDHMInvalidateCurrentStep(viewer);dCHK(err);
  }
  dhm->stepnumber = step;
  dFunctionReturn(0);
}

typedef struct {
  dInt link;
  dInt nlinks;
  dReal *times;
} step_traversal_t;

static herr_t step_traverse_func(hid_t dUNUSED loc_id,const char *name,const H5L_info_t dUNUSED *info,void *data)
{
  step_traversal_t *ctx = data;

  printf("link[%d] = '%s'\n",ctx->link,name);
  ctx->link++;
  return 0;
}

dErr dViewerDHMGetSteps(PetscViewer viewer,dInt *nsteps,dReal **steptimes)
{
  dViewer_DHM      *dhm = viewer->data;
  dErr             err;
  dTruth           match;
  H5G_info_t       info;
  step_traversal_t ctx;
  hsize_t          idx;
  herr_t           herr;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_DHM,&match);dCHK(err);
  if (!match) dERROR(PETSC_ERR_ARG_WRONG,"Viewer type must be DHM");
  err = dViewerDHMSetUp(viewer);dCHK(err);

  herr = H5Gget_info(dhm->steproot,&info);dH5CHK(herr,H5Gget_info);
  err = dMallocA(info.nlinks,steptimes);dCHK(err);
  ctx.link   = 0;               /* Will count the number of interesting links */
  ctx.nlinks = (dInt)info.nlinks;
  ctx.times  = *steptimes;
  idx = 0;
  herr = H5Literate(dhm->steproot,H5_INDEX_NAME,H5_ITER_INC,&idx,step_traverse_func,&ctx);dCHK(err);
  dhm->totalsteps = ctx.link;
  *nsteps = ctx.link;
  dFunctionReturn(0);
}

dErr dViewerDHMRestoreSteps(PetscViewer viewer,dInt *nsteps,dReal **steptimes)
{
  dTruth match;
  dErr err;

  dFunctionBegin;
  dValidHeader(viewer,PETSC_VIEWER_COOKIE,1);
  dValidIntPointer(nsteps,2);
  dValidPointer(steptimes,3);
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_DHM,&match);dCHK(err);
  if (!match) dERROR(PETSC_ERR_ARG_WRONG,"Viewer type must be DHM");
  *nsteps = 0;
  err = dFree(*steptimes);dCHK(err);
  dFunctionReturn(0);
}
