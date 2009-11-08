#include <dohpmesh.h>
#include "cont.h"
#include "../../../viewer/dhm.h"

#define H5CHK(hret,func) if (hret < 0) dERROR(PETSC_ERR_LIB, #func)

static dErr AttributeStringWrite(hid_t grp,const char *attname,const char *str)
{
  hid_t memtype,filetype,space,att;
  herr_t herr;

  dFunctionBegin;
  filetype = H5Tcopy(H5T_FORTRAN_S1);H5CHK(filetype,H5Tcopy);
  memtype = H5Tcopy(H5T_C_S1);H5CHK(memtype,H5Tcopy);
  herr = H5Tset_size(filetype,H5T_VARIABLE);H5CHK(herr,H5Tset_size);
  herr = H5Tset_size(memtype,H5T_VARIABLE);H5CHK(herr,H5Tset_size);
  space = H5Screate(H5S_SCALAR);H5CHK(space,H5Screate);
  att = H5Acreate(grp,attname,filetype,space,H5P_DEFAULT,H5P_DEFAULT);H5CHK(att,H5Acreate);
  herr = H5Awrite(att,memtype,&str);H5CHK(herr,H5Awrite);
  herr = H5Aclose(att);H5CHK(herr,H5Aclose);
  herr = H5Sclose(space);H5CHK(herr,H5Sclose);
  herr = H5Tclose(filetype);H5CHK(herr,H5Tclose);
  herr = H5Tclose(memtype);H5CHK(herr,H5Tclose);
  dFunctionReturn(0);
}

static dErr WriteTagNameAsAttribute(hid_t grp,const char *attname,dMesh mesh,dMeshTag tag)
{
  dErr err;
  char *tagname;

  dFunctionBegin;
  err = dMeshGetTagName(mesh,tag,&tagname);dCHK(err);
  err = AttributeStringWrite(grp,attname,tagname);dCHK(err);
  err = dFree(tagname);dCHK(err);
  dFunctionReturn(0);
}

static dErr WriteDimensions(hid_t grp,const char *name,dReal scale,const char *units)
{
  herr_t herr;
  hid_t dimgrp,satt,sspace;
  dErr err;

  dFunctionBegin;
  dimgrp = H5Gcreate(grp,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);                   if (dimgrp < 0) dERROR(PETSC_ERR_LIB,"H5Gcreate");
  sspace = H5Screate(H5S_SCALAR); if (sspace < 0) dERROR(PETSC_ERR_LIB,"H5Screate");  if (sspace < 0) dERROR(PETSC_ERR_LIB,"H5Screate"); 
  satt = H5Acreate(dimgrp,"scale",H5T_NATIVE_DOUBLE,sspace,H5P_DEFAULT,H5P_DEFAULT);if (satt < 0) dERROR(PETSC_ERR_LIB,"H5Acreate");
  herr = H5Awrite(satt,dH5T_REAL,&scale);  if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Awrite");
  err = AttributeStringWrite(dimgrp,"units",units);dCHK(err);
  herr = H5Aclose(satt);   if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Aclose");
  herr = H5Sclose(sspace); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Sclose");
  herr = H5Gclose(dimgrp); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");
  dFunctionReturn(0);
}


dErr dFSView_Cont_DHM(dFS fs,dViewer viewer)
{
  /* dFS_Cont *cont = fs->data; */
  dViewer_DHM *dhm = viewer->data;
  const char *meshname,*fsname;
  char mstatestr[16],fstatestr[16];
  dInt meshstate,fsstate;
  herr_t herr;
  htri_t hflg;
  hid_t fsgrp,fsgrp1,meshgrp,att,dspace;
  dErr err;
  dIInt ierr;

  dFunctionBegin;
  err = dViewerDHMSetUp(viewer);dCHK(err);
  /* Check if current mesh has been written */
  err = PetscObjectGetName((dObject)fs->mesh,&meshname);dCHK(err);
  err = PetscObjectStateQuery((dObject)fs->mesh,&meshstate);dCHK(err);
  hflg = H5Lexists(dhm->meshroot,meshname,H5P_DEFAULT); if (hflg < 0) dERROR(PETSC_ERR_LIB,"H5Lexists");
  if (!hflg) {
    meshgrp = H5Gcreate(dhm->meshroot,meshname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); if (meshgrp < 0) dERROR(PETSC_ERR_LIB,"H5Gopen");
  } else {
    meshgrp = H5Gopen(dhm->meshroot,meshname,H5P_DEFAULT); if (meshgrp < 0) dERROR(PETSC_ERR_LIB,"H5Gopen");
  }
  err = PetscSNPrintf(mstatestr,sizeof mstatestr,"%03d",meshstate);dCHK(err);
  hflg = H5Lexists(meshgrp,mstatestr,H5P_DEFAULT); if (hflg < 0) dERROR(PETSC_ERR_LIB,"H5Lexists");
  if (!hflg) {                  /* Save mesh to external file and create link to it */
    char imeshpath[dNAME_LEN];
    iMesh_Instance mi;
    err = PetscSNPrintf(imeshpath,sizeof imeshpath,"imesh-%s-%03d.h5m",meshname,meshstate);dCHK(err);
    err = dMeshGetInstance(fs->mesh,&mi);dCHK(err);
    iMesh_save(mi,0,imeshpath,"",&ierr,sizeof imeshpath,0);dICHK(mi,ierr);
    herr = H5Lcreate_external(imeshpath,"tstt",meshgrp,mstatestr,H5P_DEFAULT,H5P_DEFAULT);
    if (hflg < 0) dERROR(PETSC_ERR_LIB,"H5Lcreate_external");
  }
  herr = H5Gclose(meshgrp); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");
  /* Create group to hold this FS */
  err = PetscObjectGetName((dObject)fs,&fsname);dCHK(err);
  err = PetscObjectStateQuery((dObject)fs,&fsstate);dCHK(err);
  hflg = H5Lexists(dhm->fsroot,fsname,H5P_DEFAULT); if (hflg < 0) dERROR(PETSC_ERR_LIB,"H5Lexists");
  if (!hflg) {
    fsgrp = H5Gcreate(dhm->fsroot,fsname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); if (fsgrp < 0) dERROR(PETSC_ERR_LIB,"H5Gcreate");
  } else {
    fsgrp = H5Gopen(dhm->fsroot,fsname,H5P_DEFAULT); if (fsgrp < 0) dERROR(PETSC_ERR_LIB,"H5Gopen");
  }
  err = PetscSNPrintf(fstatestr,sizeof fstatestr,"%03d",fsstate);dCHK(err);
  fsgrp1 = H5Gcreate(fsgrp,fstatestr,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); if (fsgrp1 < 0) dERROR(PETSC_ERR_LIB,"H5Gcreate");

  /* write timestamp */
  dspace = H5Screate(H5S_SCALAR); if (dspace < 0) dERROR(PETSC_ERR_LIB,"H5Screate");
  att = H5Acreate(fsgrp1,"time",dH5T_REAL,dspace,H5P_DEFAULT,H5P_DEFAULT); if (att < 0) dERROR(PETSC_ERR_LIB,"H5Acreate");
  herr = H5Awrite(att,dH5T_REAL,&dhm->time);if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Awrite");
  herr = H5Aclose(att);                     if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Aclose");
  herr = H5Sclose(dspace);                  if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Sclose");

  /* \todo We need a way to identify the active set in MOAB's file if the FS was only defined on a subset. */

  /* Write names of some mesh tags */
  err = WriteTagNameAsAttribute(fsgrp1,"global_offset",fs->mesh,fs->gcoffsetTag);dCHK(err); /* The "closure" */
  err = WriteTagNameAsAttribute(fsgrp1,"degree",fs->mesh,fs->degreetag);dCHK(err);
  /* Units, \todo not yet stored in the FS so we write defaults here */
  err = WriteDimensions(fsgrp1,"velocity",1,"m s-1");dCHK(err);

  herr = H5Gclose(fsgrp1); if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");
  herr = H5Gclose(fsgrp);  if (herr < 0) dERROR(PETSC_ERR_LIB,"H5Gclose");
  dFunctionReturn(0);
}
