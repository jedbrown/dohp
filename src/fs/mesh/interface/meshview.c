#include <dohpmeshimpl.h>
#include <dohpviewerdhm.h>
#include <dohp.h>
#include <dohpstring.h>
#include <ctype.h>              /* needed for isprint() */

static dErr dMeshView_DHM(dMesh mesh,PetscViewer viewer)
{
  dErr err;
  const char *meshname;
  char mstatestr[16];
  dInt meshstate;
  hid_t meshgrp;
  htri_t hflg;
  dIInt ierr;
  herr_t herr;

  dFunctionBegin;
  /* Check if current mesh has been written */
  err = PetscObjectGetName((dObject)mesh,&meshname);dCHK(err);
  err = PetscObjectStateQuery((dObject)mesh,&meshstate);dCHK(err);
  err = dViewerDHMGetMeshGroup(viewer,mesh,&meshgrp);dCHK(err);
  err = PetscSNPrintf(mstatestr,sizeof mstatestr,"%03d",meshstate);dCHK(err);
#if defined dH5_USE_EXTERNAL_LINK
  hflg = H5Lexists(meshgrp,mstatestr,H5P_DEFAULT);dH5CHK(hflg,H5Lexists);
#else
  hflg = H5Aexists(meshgrp,mstatestr);dH5CHK(hflg,H5Lexists);
  hflg = 0;                     /* Why can't I just check if the dataset exists? */
#endif
  if (!hflg) {                  /* Save mesh to external file and create link to it */
    const char *dhmpath;
    char dhmpathbuf[dMAX_PATH_LEN],imeshpath[dMAX_PATH_LEN],*imeshpath_ptr = imeshpath;
    iMesh_Instance mi;
    dInt slash,dot;
    err = PetscViewerFileGetName(viewer,&dhmpath);dCHK(err);
    err = dStrcpyS(dhmpathbuf,sizeof dhmpathbuf,dhmpath);dCHK(err);
    err = dFilePathSplit(dhmpathbuf,&slash,&dot);dCHK(err);
    dhmpathbuf[dot] = 0; // gets everything but the final '.'
    err = PetscSNPrintf(imeshpath,sizeof imeshpath,"%s-imesh-%s-%03d.h5m",dhmpathbuf,meshname,meshstate);dCHK(err);
    err = dMeshGetInstance(mesh,&mi);dCHK(err);
    iMesh_save(mi,0,imeshpath,"",&ierr,(int)strlen(imeshpath),0);dICHK(mi,ierr);
#if defined dH5_USE_EXTERNAL_LINK
    herr = H5Lcreate_external(imeshpath,"tstt",meshgrp,mstatestr,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(hflg,H5Lcreate_external);
#else
    {
      hid_t fstring,mstring,sspace,strattr;
      err = dViewerDHMGetStringTypes(viewer,&fstring,&mstring,&sspace);dCHK(err);
      hflg = H5Lexists(meshgrp,mstatestr,H5P_DEFAULT);dH5CHK(hflg,H5Lexists);
      if (!hflg) {
        strattr = H5Dcreate(meshgrp,mstatestr,fstring,sspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);dH5CHK(strattr,H5Dcreate);
        herr = H5Dwrite(strattr,mstring,sspace,sspace,H5P_DEFAULT,&imeshpath_ptr);dH5CHK(herr,H5Dwrite);
        herr = H5Dclose(strattr);dH5CHK(strattr,H5Dclose);
      }
    }
#endif
  }
  herr = H5Gclose(meshgrp);dH5CHK(herr,H5Gclose);
  dFunctionReturn(0);
}

dErr dMeshView(dMesh m,PetscViewer viewer)
{
  const char *type;
  dBool iascii,idhm;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(m,dMESH_CLASSID,1);
  if (!viewer) {
    err = PetscViewerASCIIGetStdout(((PetscObject)m)->comm,&viewer);dCHK(err);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(m,1,viewer,2);
  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);dCHK(err);
  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERDHM,&idhm);dCHK(err);
  if (iascii) {
    err = PetscObjectGetType((PetscObject)m,&type);dCHK(err);
    if (((PetscObject)m)->prefix) {
      err = PetscViewerASCIIPrintf(viewer,"dMesh object:(%s)\n",((PetscObject)m)->prefix);dCHK(err);
    } else {
      err = PetscViewerASCIIPrintf(viewer,"dMesh object:\n");dCHK(err);
    }
    err = PetscViewerASCIIPrintf(viewer,"Mesh type: %s\n",(type ? type : "not yet set"));dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Internal count by type: V=%d E=%d F=%d R=%d\n",m->v.s,m->e.s,m->f.s,m->r.s);dCHK(err);
    err = dMeshSetView(m,m->root,viewer);dCHK(err);
    if (m->ops->view) {
      err = PetscViewerASCIIPushTab(viewer);dCHK(err);
      err = (*m->ops->view)(m,viewer);dCHK(err);
      err = PetscViewerASCIIPopTab(viewer);dCHK(err);
    }
  } else if (idhm) {
    err = dMeshView_DHM(m,viewer);dCHK(err);
    if (m->ops->view) {err = (*m->ops->view)(m,viewer);dCHK(err);}
  } else {
    if (m->ops->view) {
      err = (*m->ops->view)(m,viewer);dCHK(err);
    } else dERROR(((PetscObject)m)->comm,PETSC_ERR_SUP,"Viewer type '%s'",((PetscObject)viewer)->type_name);
  }
  dFunctionReturn(0);
}

dErr dMeshSetView(dMesh m,dMeshESH root,PetscViewer viewer)
{
  size_t valuesLen = 256;
  char values[256];
  iMesh_Instance mi = m->mi;
  char *tagname,*z;
  int tagtype,tagsize,intdata;
  double dbldata;
  dMeshEH ehdata;
  MeshListTag tag=MLZ;
  MeshListData data=MLZ;
  MeshListESH esh=MLZ;
  dInt i,j,ntopo;
  dBool canprint,flg;
  dErr err;

  dFunctionBegin;
  if (!viewer) {err = PetscViewerASCIIGetStdout(((dObject)m)->comm,&viewer);dCHK(err);}
  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&flg);dCHK(err);
  if (!flg) dFunctionReturn(0);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  {
    for (i=iMesh_POINT; i<iMesh_ALL_TOPOLOGIES; i++) {
    iMesh_getNumOfTopo(mi,root,i,&ntopo,&err);dICHK(mi,err);
      if (ntopo) {
        err = PetscViewerASCIIPrintf(viewer,"%20s : %d\n",dMeshEntTopologyName(i),ntopo);dCHK(err);
      }
    }
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);

  iMesh_getAllEntSetTags(mi,root,&tag.v,&tag.a,&tag.s,&err);dICHK(mi,err);
  err = PetscViewerASCIIPrintf(viewer,"Number of tags %d\n",tag.s);dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  {
    for (i=0; i<tag.s; i++) {
      err = dMeshGetTagName(m,tag.v[i],&tagname);dCHK(err);
      iMesh_getTagType(mi,tag.v[i],&tagtype,&err);dICHK(mi,err);
      iMesh_getTagSizeValues(mi,tag.v[i],&tagsize,&err);dICHK(mi,err);
      switch (tagtype) {        /* this needs a refactor */
        case iBase_INTEGER:
          iMesh_getEntSetIntData(mi,root,tag.v[i],&intdata,&err);dICHK(mi,err);
          err = PetscSNPrintf(values,valuesLen,"%d",intdata);dCHK(err);
          break;
        case iBase_DOUBLE:
          iMesh_getEntSetDblData(mi,root,tag.v[i],&dbldata,&err);dICHK(mi,err);
          err = PetscSNPrintf(values,valuesLen,"%f",dbldata);dCHK(err);
          break;
        case iBase_ENTITY_HANDLE:
          iMesh_getEntSetEHData(mi,root,tag.v[i],&ehdata,&err);dICHK(mi,err);
          err = PetscSNPrintf(values,valuesLen,"%p",ehdata);dCHK(err);
          break;
        case iBase_BYTES:
          iMesh_getEntSetData(mi,root,tag.v[i],&data.v,&data.a,&data.s,&err);dICHK(mi,err);
          canprint = PETSC_TRUE;
          for (j=0; j<data.s && ((char*)data.v)[j]; j++) {
            if (!isprint(((char*)data.v)[i])) canprint = PETSC_FALSE;
          }
          if (canprint && false) {
            err = PetscSNPrintf(values,(size_t)data.s,"%s",data.v);dCHK(err); /* Just a copy, but ensures a NULL byte */
          } else {
            z = values;
            for (j=0; j<data.s && ((char*)data.v)[j] && (size_t)(z-values) < valuesLen-5; j++) {
              err = PetscSNPrintf(z,3,"%02uhhx ",((char*)data.v)[j]);dCHK(err);
              z += 3;
              if (j%4 == 0) {
                *(z++) = ' ';
              }
              *(z++) = '\0';       /* Terminate the string */
            }
          }
          err = MeshListFree(data);dCHK(err);
          break;
        default: dERROR(PETSC_COMM_SELF,1,"Invalid tag type, iMesh probably corrupt");
      }
      err = PetscViewerASCIIPrintf(viewer,"Tag: %30s : %20s [%3d] = %s\n",tagname,iBase_TagValueTypeName[tagtype],tagsize,values);dCHK(err);
      err = PetscFree(tagname);dCHK(err);
    }
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  err = MeshListFree(tag);dCHK(err);

  iMesh_getEntSets(mi,root,1,&esh.v,&esh.a,&esh.s,&err);dICHK(mi,err);
  err = PetscViewerASCIIPrintf(viewer,"Number of contained Entity Sets: %d\n",esh.s);dCHK(err);

  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  for (i=0; i<esh.s; i++) {
    err = PetscViewerASCIIPrintf(viewer,"Contained set %d/%d:\n",i+1,esh.s);dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = dMeshSetView(m,esh.v[i],viewer);dCHK(err);
    err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  err = MeshListFree(esh);dCHK(err);

  iMesh_getChldn(mi,root,1,&esh.v,&esh.a,&esh.s,&err);dICHK(mi,err);
  err = PetscViewerASCIIPrintf(viewer,"Number of child Entity Sets: %d\n",esh.s);dCHK(err);

  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  for (i=0; i<esh.s; i++) {
    err = PetscViewerASCIIPrintf(viewer,"Child %d/%d:\n",i+1,esh.s);dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = dMeshSetView(m,esh.v[i],viewer);dCHK(err);
    err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  err = MeshListFree(esh);dCHK(err);
  dFunctionReturn(0);
}

