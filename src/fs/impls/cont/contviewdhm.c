#include <dohpfs.h>
#include <dohpvec.h>
#include "cont.h"
#include "../../../viewer/dhm.h"

EXTERN dErr VecView_Dohp_FSCont(Vec,dViewer);

/** Get data set and space to write the next element for this FS
* \note caller is responsible for closing both
**/
static dErr dFSGetDHMLink(dFS fs,dViewer viewer,hid_t *indset,hid_t *inspace)
{
  dViewer_DHM *dhm = viewer->data;
  const char  *fsname;
  htri_t       hflg;
  hsize_t dims[1] = {1},maxdims[1] = {H5S_UNLIMITED};
  hid_t        dset,space,h5t_fs;
  herr_t       herr;
  dErr         err;

  dFunctionBegin;
  err = PetscObjectGetName((dObject)fs,&fsname);dCHK(err);
  hflg = H5Lexists(dhm->fsroot,fsname,H5P_DEFAULT);dH5CHK(hflg,H5Lexists);
  if (!hflg) {
    hsize_t chunk[1] = {1};
    hid_t dcpl;
    space = H5Screate_simple(1,dims,maxdims);dH5CHK(space,H5Screate_simple);
    err = dViewerDHMGetFSType(viewer,&h5t_fs);dCHK(err);
    dcpl = H5Pcreate(H5P_DATASET_CREATE);dH5CHK(dcpl,H5Pcreate);
    herr = H5Pset_chunk(dcpl,1,chunk);dH5CHK(herr,H5Pset_chunk);
    dset = H5Dcreate(dhm->fsroot,fsname,h5t_fs,space,H5P_DEFAULT,dcpl,H5P_DEFAULT);dH5CHK(dset,H5Dcreate);
    herr = H5Pclose(dcpl);dH5CHK(herr,H5Pclose);
  } else {
    dset = H5Dopen(dhm->fsroot,fsname,H5P_DEFAULT);dH5CHK(dset,H5Dopen);
    space = H5Dget_space(dset);dH5CHK(space,H5Dget_space);
    herr = H5Sget_simple_extent_dims(space,dims,NULL);dH5CHK(herr,H5Sget_simple_extent_dims);
#if 0                           /* Handle case where it has not already been written in this state */
    herr = H5Sclose(space);dH5CHK(herr,H5Sclose);
    /* Extend by one */
    dims[0]++;
    herr = H5Dset_extent(dset,dims);dH5CHK(herr,H5Dset_extent);
    /* Select the last entry */
    dims[0]--;
#endif
    space = H5Dget_space(dset);dH5CHK(space,H5Dget_space);
    herr = H5Sselect_elements(space,H5S_SELECT_SET,1,dims);dH5CHK(herr,H5Sselect_elements);
  }
  *indset = dset;
  *inspace = space;
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
  hid_t meshgrp;
  dErr err;
  dIInt ierr;

  dFunctionBegin;
  err = dViewerDHMSetUp(viewer);dCHK(err);
  /* Check if current mesh has been written */
  err = PetscObjectGetName((dObject)fs->mesh,&meshname);dCHK(err);
  err = PetscObjectStateQuery((dObject)fs->mesh,&meshstate);dCHK(err);
  hflg = H5Lexists(dhm->meshroot,meshname,H5P_DEFAULT);dH5CHK(hflg,H5Lexists);
  if (!hflg) {
    meshgrp = H5Gcreate(dhm->meshroot,meshname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(meshgrp,H5Gcreate);
  } else {
    meshgrp = H5Gopen(dhm->meshroot,meshname,H5P_DEFAULT);dH5CHK(meshgrp,H5Gopen);
  }
  err = PetscSNPrintf(mstatestr,sizeof mstatestr,"%03d",meshstate);dCHK(err);
#if defined dH5_USE_EXTERNAL_LINK
  hflg = H5Lexists(meshgrp,mstatestr,H5P_DEFAULT);dH5CHK(hflg,H5Lexists);
#else
  hflg = H5Aexists(meshgrp,mstatestr);dH5CHK(hflg,H5Lexists);
  hflg = 0;                     /* Why can't I just check if the dataset exists? */
#endif
  if (!hflg) {                  /* Save mesh to external file and create link to it */
    char imeshpath[dNAME_LEN],*imeshpath_ptr = imeshpath;
    iMesh_Instance mi;
    err = PetscSNPrintf(imeshpath,sizeof imeshpath,"imesh-%s-%03d.h5m",meshname,meshstate);dCHK(err);
    err = dMeshGetInstance(fs->mesh,&mi);dCHK(err);
    iMesh_save(mi,0,imeshpath,"",&ierr,(int)strlen(imeshpath),0);dICHK(mi,ierr);
#if defined dH5_USE_EXTERNAL_LINK
    herr = H5Lcreate_external(imeshpath,"tstt",meshgrp,mstatestr,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(hflg,H5Lcreate_external);
#else
    {
      hid_t fstring,mstring,sspace,strattr;
      err = dViewerDHMGetStringTypes(viewer,&fstring,&mstring,&sspace);dCHK(err);
      strattr = H5Dcreate(meshgrp,mstatestr,fstring,sspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(strattr,H5Dcreate);
      herr = H5Dwrite(strattr,mstring,sspace,sspace,H5P_DEFAULT,&imeshpath_ptr);dH5CHK(herr,H5Dwrite);
      herr = H5Dclose(strattr);dH5CHK(strattr,H5Dclose);
    }
#endif
  }

  {
    dht_FS     fs5;
    dht_Field *field5;
    hid_t      h5t_fs,fsdset,fsspace;
    dInt       i;
    PetscMPIInt size;

    err = dViewerDHMGetFSType(viewer,&h5t_fs);dCHK(err);
    err = dMeshGetTagName(fs->mesh,fs->degreetag,&fs5.degree);dCHK(err);
    err = dMeshGetTagName(fs->mesh,fs->gcoffsetTag,&fs5.global_offset);dCHK(err);
    err = dMeshGetTagName(fs->mesh,fs->partitionTag,&fs5.partition);dCHK(err);
    err = dMeshGetTagName(fs->mesh,fs->orderedsubTag,&fs5.ordered_subdomain);dCHK(err);
    herr = H5Rcreate(&fs5.mesh,meshgrp,mstatestr,H5R_OBJECT,-1);dH5CHK(herr,H5Rcreate);
    fs5.time = dhm->time;
    err = PetscObjectStateQuery((PetscObject)fs,&fs5.internal_state);dCHK(err);
    err = MPI_Comm_size(((dObject)fs)->comm,&size);dCHK(err);
    fs5.number_of_subdomains = size;
    fs5.fields.len = fs->bs;
    err = dMallocA(fs5.fields.len,&field5);dCHK(err);
    for (i=0; i<fs->bs; i++) {
      field5[i].name = fs->fieldname[i];
      field5[i].units.dimensions = (char*)"m s-1"; /* we only use it as \c const */
      field5[i].units.scale = exp(1);
    }
    fs5.fields.p = field5;
    err = dFSGetDHMLink(fs,viewer,&fsdset,&fsspace);dCHK(err); /* Get location to write this FS */
    herr = H5Dwrite(fsdset,h5t_fs,H5S_ALL,H5S_ALL,H5P_DEFAULT,&fs5);dH5CHK(herr,H5Dwrite);
    err = dFree(field5);dCHK(err);
    err = dFree(fs5.partition);dCHK(err);
    err = dFree(fs5.ordered_subdomain);dCHK(err);
    err = dFree(fs5.global_offset);dCHK(err);
    err = dFree(fs5.degree);dCHK(err);
    herr = H5Dclose(fsdset);dH5CHK(herr,H5Dclose);
    herr = H5Sclose(fsspace);dH5CHK(herr,H5Sclose);
  }

  /* \todo We need a way to identify the active set in MOAB's file if the FS was only defined on a subset. */

  herr = H5Gclose(meshgrp);dH5CHK(herr,H5Gclose);
  dFunctionReturn(0);
}

static dErr VecView_Dohp_FSCont_DHM(Vec X,PetscViewer viewer)
{
  dViewer_DHM *dhm = viewer->data;
  const char  *xname;
  dFS          fs;
  Vec          Xclosure;
  dScalar     *x;
  hid_t        fsdset,fsspace,curstep,dset,fspace,mspace,plist;
  hsize_t      gdim[2],offset[2],count[2];
  herr_t       herr;
  dInt         m,low,high,bs;
  dErr         err;

  dFunctionBegin;
  err = dViewerDHMSetUp(viewer);dCHK(err);
  err = PetscObjectGetName((PetscObject)X,&xname);dCHK(err);
  err = PetscObjectQuery((PetscObject)X,"dFS",(PetscObject*)&fs);dCHK(err);
  if (!fs) dERROR(PETSC_ERR_ARG_WRONG,"Vector not generated from a FS");
  err = dFSGetDHMLink(fs,viewer,&fsdset,&fsspace);dCHK(err); /* we are responsible for closing */
  err = dViewerDHMGetStep(viewer,&curstep);dCHK(err);        /* leave curstep open */
  err = dFSView_Cont_DHM(fs,viewer);dCHK(err);
  err = VecDohpGetClosure(X,&Xclosure);dCHK(err);
  err = VecGetSize(Xclosure,&m);dCHK(err);
  err = VecGetOwnershipRange(Xclosure,&low,&high);dCHK(err);
  err = VecGetBlockSize(Xclosure,&bs);dCHK(err);
  gdim[0]   = m/bs;
  gdim[1]   = bs;
  offset[0] = low/bs;
  offset[1] = 0;
  count[0]  = (high-low)/bs;
  count[1]  = bs;
  fspace = H5Screate_simple(2,gdim,NULL);dH5CHK(fspace,H5Screate_simple);
  dset = H5Dcreate(curstep,xname,dH5T_SCALAR,fspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(dset,H5Dcreate);
  herr = H5Sselect_hyperslab(fspace,H5S_SELECT_SET,offset,NULL,count,NULL);dH5CHK(herr,H5Sselect_hyperslab);
  mspace = H5Screate_simple(2,count,NULL);dH5CHK(mspace,H5Screate_simple);

  plist = H5Pcreate(H5P_DATASET_XFER);dH5CHK(plist,H5Pcreate);
#if defined dUSE_PARALLEL_HDF5
  herr = H5Pset_dxpl_mpio(plist,H5FD_MPIO_COLLECTIVE);dH5CHK(herr,H5Pset_dxpl_mpio);
#endif

  err = VecGetArray(Xclosure,&x);dCHK(err);
  herr = H5Dwrite(dset,dH5T_SCALAR,mspace,fspace,plist,x);dH5CHK(herr,H5Dwrite);
  err = VecRestoreArray(Xclosure,&x);dCHK(err);
  err = VecDohpRestoreClosure(X,&Xclosure);dCHK(err);

  /* Write attributes on this dataset */
  {
    char fsname[256];
    ssize_t namelen;
    namelen = H5Iget_name(fsdset,fsname,sizeof fsname);dH5CHK(namelen,H5Iget_name);
    if (!namelen) dERROR(PETSC_ERR_LIB,"Could not get FS path");
    {
      dht_Vec vecatt[1];
      hsize_t dims[1] = {1};
      hid_t aspace,attr,vectype;

      err = dViewerDHMGetVecType(viewer,&vectype);dCHK(err);
      /* Initialize data vecatt[0] */
      herr = H5Rcreate(&vecatt[0].fs,dhm->file,fsname,H5R_DATASET_REGION,fsspace);dH5CHK(herr,H5Rcreate);
      vecatt[0].time = dhm->time;
      err = PetscObjectStateQuery((PetscObject)X,&vecatt[0].internal_state);dCHK(err);
      /* Write attribute */
      aspace = H5Screate_simple(1,dims,NULL);dH5CHK(aspace,H5Screate_simple);
      attr = H5Acreate(dset,"meta",vectype,aspace,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(attr,H5Acreate);
      herr = H5Awrite(attr,vectype,vecatt);dH5CHK(herr,H5Awrite);
      herr = H5Aclose(attr);dH5CHK(herr,H5Aclose);
      herr = H5Sclose(aspace);dH5CHK(herr,H5Aclose);
    }
  }

  herr = H5Dclose(dset);dH5CHK(herr,H5Dclose);
  herr = H5Pclose(plist);dH5CHK(herr,H5Pclose);
  herr = H5Sclose(fspace);dH5CHK(herr,H5Sclose);
  herr = H5Sclose(mspace);dH5CHK(herr,H5Sclose);
  herr = H5Dclose(fsdset);dH5CHK(herr,H5Dclose);
  herr = H5Sclose(fsspace);dH5CHK(herr,H5Sclose);
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

dErr dFSLoadIntoFS_Cont_DHM(PetscViewer viewer,const char fieldname[],dFS fs)
{
  dViewer_DHM *dhm = viewer->data;
  dErr        err;
  hid_t       curstep,vectype,vdset,vattr,fsobj,fsspace;
  herr_t      herr;
  dht_Vec     vecmeta;
  //dht_FS      fsmeta;

  dFunctionBegin;
  err = dViewerDHMSetUp(viewer);dCHK(err);
  err = dViewerDHMGetStep(viewer,&curstep);dCHK(err);
  err = dViewerDHMGetVecType(viewer,&vectype);dCHK(err);

  vdset = H5Dopen(curstep,fieldname,H5P_DEFAULT);dH5CHK(vdset,H5Dopen);
  vattr = H5Aopen(vdset,"meta",H5P_DEFAULT);dH5CHK(vattr,H5Aopen);
  herr = H5Aread(vattr,vectype,&vecmeta);dH5CHK(herr,H5Aread);
  herr = H5Aclose(vattr);dH5CHK(herr,H5Aclose);

  err = dPrintf(PETSC_COMM_SELF,"Vec name '%s'  time %g  internal_state %d\n",fieldname,vecmeta.time,vecmeta.internal_state);dCHK(err);

  fsobj = H5Rdereference(vdset,H5R_DATASET_REGION,vecmeta.fs);dH5CHK(fsobj,H5Rdereference);
  fsspace = H5Rget_region(vdset,H5R_DATASET_REGION,vecmeta.fs);dH5CHK(fsobj,H5Rget_region);
  {
    char fsobjname[256];
    ssize_t len;
    hssize_t nrec;
    len = H5Iget_name(fsobj,fsobjname,sizeof fsobjname);dH5CHK(len,H5Iget_name);
    nrec = H5Sget_select_npoints(fsspace);dH5CHK(nrec,H5Sget_select_npoints);
    err = dPrintf(PETSC_COMM_SELF,"fsobj name '%s', npoints %zd\n",len?fsobjname:"(no name)",nrec);dCHK(err);
  }
  {
    dht_FS fs5;
    hid_t memspace,h5t_FS,meshobj;
    err = dViewerDHMGetFSType(viewer,&h5t_FS);dCHK(err);
    memspace = H5Screate(H5S_SCALAR);
    err = H5Dread(fsobj,h5t_FS,memspace,fsspace,H5P_DEFAULT,&fs5);dH5CHK(herr,H5Dread);
    herr = H5Sclose(memspace);

    printf("degree = %s\nglobal_offset = %s\npartition = %s\nordered_subdomain = %s\n",fs5.degree,fs5.global_offset,fs5.partition,fs5.ordered_subdomain);
    meshobj = H5Rdereference(dhm->meshroot,H5R_OBJECT,&fs5.mesh);dH5CHK(meshobj,H5Rdereference);
    {
      char meshname[dNAME_LEN] = {0};
      ssize_t len;
      len = H5Iget_name(meshobj,meshname,sizeof meshname);dH5CHK(len,H5Iget_name);
      printf("mesh name = %s\n",meshname);
    }

    {
      hid_t mstring,strspace;
      dMesh mesh;
      char *imeshstr;           /* We are reading to a vlen string so HDF5 will allocate memory */
      err = dViewerDHMGetStringTypes(viewer,NULL,&mstring,&strspace);dCHK(err);
      herr = H5Dread(meshobj,mstring,H5S_ALL,H5S_ALL,H5P_DEFAULT,&imeshstr);dH5CHK(herr,H5Dread);
      printf("imeshstr = %s\n",imeshstr);
      err = dFSGetMesh(fs,&mesh);dCHK(err);
      err = dMeshSetInFile(mesh,imeshstr,NULL);dCHK(err);
      err = dMeshSetType(mesh,dMESHSERIAL);dCHK(err);
      err = dMeshLoad(mesh);dCHK(err);
      {
        dMeshTag tag;
        dJacobi jac;
        dIInt readrank = 0;         /* Hard-code the rank for now */
        err = dMeshGetTag(mesh,fs5.partition,&tag);dCHK(err);
        err = dMeshGetTaggedSet(mesh,tag,&readrank,&fs->activeSet);dCHK(err);
        err = dMeshGetTag(mesh,fs5.ordered_subdomain,&tag);dCHK(err);
        err = dMeshGetTaggedSet(mesh,tag,&readrank,&fs->orderedSet);dCHK(err);
        err = dMeshGetTag(mesh,fs5.degree,&tag);dCHK(err);
        err = dFSGetJacobi(fs,&jac);dCHK(err);
        err = dFSSetDegree(fs,jac,tag);dCHK(err);
        err = dMeshGetTag(mesh,fs5.global_offset,&fs->gcoffsetTag);dCHK(err);
      }
      herr = H5Dvlen_reclaim(mstring,strspace,H5P_DEFAULT,&imeshstr);dH5CHK(herr,H5Dvlen_reclaim);
    }
    herr = H5Dclose(meshobj);dH5CHK(herr,H5Aclose);
  }
  /** @note The FS has the layout and ordering tags set (@todo boundary conditions) so we are ready to build the
  * function space (mostly creating the LocalToGlobalMapping/update for VecDohp and setting up the element assembly
  * matrix E.
  **/

  /* @todo Call private dFSBuildSpace pieces (once they exist) */

  herr = H5Sclose(fsspace);dH5CHK(herr,H5Sclose);
  herr = H5Oclose(fsobj);dH5CHK(herr,H5Oclose);
  herr = H5Dclose(vdset);dH5CHK(herr,H5Dclose);
  dFunctionReturn(0);
}