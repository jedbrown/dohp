#include <dohpfs.h>
#include <dohpvec.h>
#include "cont.h"
#include "../../../viewer/dhm.h"

static dErr WriteTagNameAsAttribute(PetscViewer viewer,hid_t grp,const char *attname,dMesh mesh,dMeshTag tag)
{
  dErr err;
  char *tagname;

  dFunctionBegin;
  err = dMeshGetTagName(mesh,tag,&tagname);dCHK(err);
  err = dViewerDHMAttributeStringWrite(viewer,grp,attname,tagname);dCHK(err);
  err = dFree(tagname);dCHK(err);
  dFunctionReturn(0);
}

static dErr dFSGetDHMLink(dFS fs,dViewer viewer,hid_t *link)
{
  dViewer_DHM *dhm = viewer->data;
  char         fstatestr[16];
  const char  *fsname;
  dInt         fsstate;
  htri_t       hflg;
  hid_t        fsgrp;
  herr_t       herr;
  dErr         err;

  dFunctionBegin;
  err = PetscObjectGetName((dObject)fs,&fsname);dCHK(err);
  err = PetscObjectStateQuery((dObject)fs,&fsstate);dCHK(err);
  hflg = H5Lexists(dhm->fsroot,fsname,H5P_DEFAULT); if (hflg < 0) dERROR(PETSC_ERR_LIB,"H5Lexists");
  if (!hflg) {
    fsgrp = H5Gcreate(dhm->fsroot,fsname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); if (fsgrp < 0) dERROR(PETSC_ERR_LIB,"H5Gcreate");
  } else {
    fsgrp = H5Gopen(dhm->fsroot,fsname,H5P_DEFAULT); if (fsgrp < 0) dERROR(PETSC_ERR_LIB,"H5Gopen");
  }
  err = PetscSNPrintf(fstatestr,sizeof fstatestr,"%03d",fsstate);dCHK(err);
  hflg = H5Lexists(fsgrp,fstatestr,H5P_DEFAULT);dH5CHK(hflg,H5Lexists);
  if (!hflg) {
    *link = H5Gcreate(fsgrp,fstatestr,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(*link,H5Gcreate);
  } else {
    *link = H5Gopen(fsgrp,fstatestr,H5P_DEFAULT);dH5CHK(*link,H5Gopen);
  }
  herr = H5Gclose(fsgrp);dH5CHK(herr,H5Gclose);
  dFunctionReturn(0);
}

static dErr dVecGetDHMLink(Vec x,dViewer viewer,hid_t *link)
{
  char         xstatestr[16];
  const char  *xname;
  dInt         xstate;
  htri_t       hflg;
  hid_t        xgrp,curstep;
  herr_t       herr;
  dErr         err;

  dFunctionBegin;
  err = PetscObjectGetName((dObject)x,&xname);dCHK(err);
  err = PetscObjectStateQuery((dObject)x,&xstate);dCHK(err);
  err = dViewerDHMGetStep(viewer,&curstep);dCHK(err);
  hflg = H5Lexists(curstep,xname,H5P_DEFAULT);dH5CHK(hflg,H5Lexists);
  if (!hflg) {
    xgrp = H5Gcreate(curstep,xname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(xgrp,H5Gcreate);
  } else {
    xgrp = H5Gopen(curstep,xname,H5P_DEFAULT);dH5CHK(xgrp,H5Gopen);
  }
  err = PetscSNPrintf(xstatestr,sizeof xstatestr,"%03d",xstate);dCHK(err);
  *link = H5Gcreate(xgrp,xstatestr,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(*link,H5Gcreate);
  herr = H5Gclose(xgrp);dH5CHK(herr,H5Gclose);
  dFunctionReturn(0);
}

dErr dFSView_Cont_DHM(dFS fs,dViewer viewer)
{
  /* dFS_Cont *cont = fs->data; */
  dViewer_DHM *dhm = viewer->data;
  const char *meshname;
  char mstatestr[16];
  dInt meshstate;
  herr_t herr;
  htri_t hflg;
  hid_t fsgrp,meshgrp,att,dspace;
  dErr err;
  dIInt ierr;

  dFunctionBegin;
  err = dViewerDHMSetUp(viewer);dCHK(err);
  /* Check if current mesh has been written */
  err = PetscObjectGetName((dObject)fs->mesh,&meshname);dCHK(err);
  err = PetscObjectStateQuery((dObject)fs->mesh,&meshstate);dCHK(err);
  hflg = H5Lexists(dhm->meshroot,meshname,H5P_DEFAULT); if (hflg < 0) dERROR(PETSC_ERR_LIB,"H5Lexists");
  if (!hflg) {
    meshgrp = H5Gcreate(dhm->meshroot,meshname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(meshgrp,H5Gcreate);
  } else {
    meshgrp = H5Gopen(dhm->meshroot,meshname,H5P_DEFAULT);dH5CHK(meshgrp,H5Gopen);
  }
  err = PetscSNPrintf(mstatestr,sizeof mstatestr,"%03d",meshstate);dCHK(err);
  hflg = H5Lexists(meshgrp,mstatestr,H5P_DEFAULT); if (hflg < 0) dERROR(PETSC_ERR_LIB,"H5Lexists");
  if (!hflg) {                  /* Save mesh to external file and create link to it */
    char imeshpath[dNAME_LEN];
    iMesh_Instance mi;
    err = PetscSNPrintf(imeshpath,sizeof imeshpath,"imesh-%s-%03d.h5m",meshname,meshstate);dCHK(err);
    err = dMeshGetInstance(fs->mesh,&mi);dCHK(err);
    iMesh_save(mi,0,imeshpath,"",&ierr,(int)strlen(imeshpath),0);dICHK(mi,ierr);
    herr = H5Lcreate_external(imeshpath,"tstt",meshgrp,mstatestr,H5P_DEFAULT,H5P_DEFAULT);
    if (hflg < 0) dERROR(PETSC_ERR_LIB,"H5Lcreate_external");
  }
  herr = H5Gclose(meshgrp); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");
  /* Create group to hold this FS */
  err = dFSGetDHMLink(fs,viewer,&fsgrp);dCHK(err);

  /* write timestamp */
  dspace = H5Screate(H5S_SCALAR); if (dspace < 0) dERROR(PETSC_ERR_LIB,"H5Screate");
  att = H5Acreate(fsgrp,"time",dH5T_REAL,dspace,H5P_DEFAULT,H5P_DEFAULT); if (att < 0) dERROR(PETSC_ERR_LIB,"H5Acreate");
  herr = H5Awrite(att,dH5T_REAL,&dhm->time);if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Awrite");
  herr = H5Aclose(att);                     if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Aclose");
  herr = H5Sclose(dspace);                  if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Sclose");

  /* \todo We need a way to identify the active set in MOAB's file if the FS was only defined on a subset. */

  /* Write names of some mesh tags */
  err = WriteTagNameAsAttribute(viewer,fsgrp,"global_offset",fs->mesh,fs->gcoffsetTag);dCHK(err); /* The "closure" */
  err = WriteTagNameAsAttribute(viewer,fsgrp,"degree",fs->mesh,fs->degreetag);dCHK(err);
  /* Units, \todo not yet stored in the FS so we write defaults here */
  err = dViewerDHMWriteDimensions(viewer,fsgrp,"velocity","m s-1",exp(1));dCHK(err);

  herr = H5Gclose(fsgrp); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");
  dFunctionReturn(0);
}

static dErr VecView_Dohp_FSCont_DHM(Vec X,PetscViewer viewer)
{
  dFS      fs;
  hid_t    fslink,xlink,fspace,mspace,dset,plist;
  hsize_t  gdim[2],offset[2],count[2];
  herr_t   herr;
  dInt     m,low,high,bs;
  dErr     err;

  dFunctionBegin;
  err = dViewerDHMSetUp(viewer);dCHK(err);
  err = PetscObjectQuery((PetscObject)X,"dFS",(PetscObject*)&fs);dCHK(err);
  if (!fs) dERROR(PETSC_ERR_ARG_WRONG,"Vector not generated from a FS");
  err = dFSGetDHMLink(fs,viewer,&fslink);dCHK(err);
  err = dVecGetDHMLink(X,viewer,&xlink);dCHK(err);
  err = dFSView_Cont_DHM(fs,viewer);dCHK(err);
  err = VecGetSize(X,&m);dCHK(err);
  err = VecGetOwnershipRange(X,&low,&high);dCHK(err);
  err = VecGetBlockSize(X,&bs);dCHK(err);
  gdim[0]   = m/bs;
  gdim[1]   = bs;
  offset[0] = low/bs;
  offset[1] = 0;
  count[0]  = (high-low)/bs;
  count[1]  = bs;
  fspace = H5Screate_simple(2,gdim,NULL);dH5CHK(fspace,H5Screate_simple);
  dset = H5Dcreate(xlink,"values",dH5T_SCALAR,fspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dset,H5Dcreate);
  herr = H5Sselect_hyperslab(fspace,H5S_SELECT_SET,offset,NULL,count,NULL);dH5CHK(herr,H5Sselect_hyperslab);
  mspace = H5Screate_simple(2,count,NULL);dH5CHK(mspace,H5Screate_simple);

  plist = H5Pcreate(H5P_DATASET_XFER);dH5CHK(plist,H5Pcreate);
  herr = H5Pset_dxpl_mpio(plist,H5FD_MPIO_COLLECTIVE);dH5CHK(herr,H5Pset_dxpl_mpio);

  {
    Vec Xclosure;
    dScalar *x;
    err = VecDohpGetClosure(X,&Xclosure);dCHK(err);
    err = VecGetArray(Xclosure,&x);dCHK(err);
    herr = H5Dwrite(dset,dH5T_SCALAR,mspace,fspace,plist,x);dH5CHK(herr,H5Dwrite);
    err = VecRestoreArray(Xclosure,&x);dCHK(err);
    err = VecDohpRestoreClosure(X,&Xclosure);dCHK(err);
  }

  {
    char fsname[256];
    ssize_t namelen;
    namelen = H5Iget_name(fslink,fsname,sizeof fsname);dH5CHK(namelen,H5Iget_name);
    if (!namelen) dERROR(PETSC_ERR_LIB,"Could not get FS path");
    herr = H5Lcreate_soft(fsname,xlink,"fs",H5P_DEFAULT,H5P_DEFAULT);dH5CHK(herr,H5Lcreate_soft);
  }

  herr = H5Dclose(dset);dH5CHK(herr,H5Dclose);
  herr = H5Pclose(plist);dH5CHK(herr,H5Pclose);
  herr = H5Sclose(fspace);dH5CHK(herr,H5Sclose);
  herr = H5Sclose(mspace);dH5CHK(herr,H5Sclose);
  herr = H5Gclose(xlink);dH5CHK(herr,H5Gclose);
  herr = H5Gclose(fslink);dH5CHK(herr,H5Gclose);
  dFunctionReturn(0);
}


dErr VecView_Dohp_FSCont(Vec x,PetscViewer viewer)
{
  dFS fs;
  dTruth isdhm,isdraw;
  dErr err;

  dFunctionBegin;
  err = PetscObjectQuery((PetscObject)x,"dFS",(PetscObject*)&fs);dCHK(err);
  if (!fs) dERROR(PETSC_ERR_ARG_WRONG,"Vector not generated from a FS");
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_DHM,&isdhm);dCHK(err);
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_DRAW,&isdraw);dCHK(err);
  if (isdhm) {
    err = VecView_Dohp_FSCont_DHM(x,viewer);dCHK(err);
  } else if (isdraw) {
    dERROR(1,"not implemented");
  } else {
    dERROR(1,"not implemented");
  }
  dFunctionReturn(0);
}
