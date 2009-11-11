#include <dohpviewer.h>
#include <dohp.h>

dErr dInitialize(int *argc,char ***argv,const char *file,const char *help)
{
  dErr err;

  dFunctionBegin;
  err = PetscInitialize(argc,argv,file,help);dCHK(err);
  err = dViewerRegisterAll(NULL);dCHK(err);
  dFunctionReturn(0);
}

dErr dFinalize(void)
{
  dErr err;

  dFunctionBegin;
  err = PetscFinalize();dCHK(err);
  dFunctionReturn(0);
}
