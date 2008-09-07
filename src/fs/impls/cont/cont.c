#include "cont.h"

static dErr dFSView_Cont(dFS fs,dViewer viewer)
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (ascii) {
    err = PetscViewerASCIIPrintf(viewer,"dFSView_Cont()  nothing here yet\n");dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dFSCreate_Cont(dFS fs)
{
  dErr err;

  dFunctionBegin;
  err = dMalloc(sizeof(dFS_Cont),&fs->data);dCHK(err);
  err = dMemzero(fs->data,sizeof(dFS_Cont));dCHK(err);
  fs->ops->view = &dFSView_Cont;
  dFunctionReturn(0);
}
