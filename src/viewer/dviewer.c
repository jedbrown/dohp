#include <dohpviewer.h>

EXTERN PetscErrorCode PetscViewerCreate_DHM(PetscViewer);

dErr dViewerRegisterAll(const char *path)
{
  dErr err;

  dFunctionBegin;
  err = PetscViewerRegister(PETSC_VIEWER_DHM,path,"PetscViewerCreate_DHM",PetscViewerCreate_DHM);dCHK(err);
  dFunctionReturn(0);
}

dErr dViewerDHMSetFS(dViewer viewer,dFS fs)
{
  dErr err,(*r)(dViewer,dFS);

  dFunctionBegin;
  err = PetscObjectQueryFunction((PetscObject)viewer,"dViewerDHMSetFS_C",(void(**)(void))&r);dCHK(err);
  if (r) {
    err = (*r)(viewer,fs);dCHK(err);
  }
  dFunctionReturn(0);
}


dErr dViewerDHMSetTime(dViewer viewer,dReal time)
{
  dErr err,(*r)(dViewer,dReal);

  dFunctionBegin;
  err = PetscObjectQueryFunction((PetscObject)viewer,"dViewerDHMSetTime_C",(void(**)(void))&r);dCHK(err);
  if (r) {
    err = (*r)(viewer,time);dCHK(err);
  }
  dFunctionReturn(0);
}
