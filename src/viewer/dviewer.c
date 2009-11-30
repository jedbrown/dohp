#include <dohpviewer.h>
#include <dohp.h>

extern PetscErrorCode PetscViewerCreate_DHM(PetscViewer);

dErr dViewerRegisterAll(const char *path)
{
  dErr err;

  dFunctionBegin;
  err = PetscViewerRegister(PETSC_VIEWER_DHM,path,"PetscViewerCreate_DHM",PetscViewerCreate_DHM);dCHK(err);
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

dErr dViewerDHMSetTimeUnits(dViewer viewer,const char *units,dReal scale)
{
  dErr err,(*r)(dViewer,const char*,dReal);

  dFunctionBegin;
  err = PetscObjectQueryFunction((PetscObject)viewer,"dViewerDHMSetTimeUnits_C",(void(**)(void))&r);dCHK(err);
  if (r) {
    err = (*r)(viewer,units,scale);dCHK(err);
  }
  dFunctionReturn(0);
}
