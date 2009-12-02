#include <dohpviewer.h>
#include <dohpsys.h>
#include <dohp.h>

static bool DohpInitialized;

dErr dInitialize(int *argc,char ***argv,const char *file,const char *help)
{
  dErr err;

  dFunctionBegin;
  if (DohpInitialized) dFunctionReturn(0);
  err = PetscInitialize(argc,argv,file,help);dCHK(err);
  err = dViewerRegisterAll(NULL);dCHK(err);
  DohpInitialized = true;
  dFunctionReturn(0);
}

dErr dFinalize(void)
{
  dErr err;

  dFunctionBegin;
  if (!DohpInitialized) {
    (*PetscErrorPrintf)("dInitialize() must be called before dFinalize()\n");
    PetscFunctionReturn(0);
  }
  /* Don't unregister our viewers, PETSc will take care of it */
  err = PetscFinalize();dCHK(err);
  DohpInitialized = false;
  dFunctionReturn(0);
}
