#include "private/dmeshimpl.h"

dErr dMeshCreate_Serial(dMesh mesh); /* The only exported function */

static dErr dMeshLoad_Serial(dMesh mesh)
{
  char options[dSTR_LEN];
  size_t fnamelen;
  dErr err;

  dFunctionBegin;
  err = PetscStrlen(mesh->infile,&fnamelen);dCHK(err);
  if (mesh->inoptions) {
    err = PetscSNPrintf(options,sizeof(options),"%s",mesh->inoptions);dCHK(err);
  } else options[0] = 0;
  iMesh_load(mesh->mi,0,mesh->infile,options,&err,(int)fnamelen,(int)strlen(options));dICHK(mesh->mi,err);
  dFunctionReturn(0);
}

dErr dMeshCreate_Serial(dMesh mesh)
{
  dErr err;

  dFunctionBegin;

  iMesh_newMesh("",&mesh->mi,&err,0);dICHK(mesh->mi,err);

  mesh->data = 0;
  mesh->ops->view = 0;
  mesh->ops->destroy = 0;
  mesh->ops->setfromoptions = 0;
  mesh->ops->load = dMeshLoad_Serial;
  mesh->ops->tagbcast = 0;
  dFunctionReturn(0);
}
