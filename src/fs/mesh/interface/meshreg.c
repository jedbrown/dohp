#include <dohpmeshimpl.h>
#include <dohp.h>

dClassId dMESH_CLASSID;
PetscBool dMeshRegisterAllCalled;
PetscLogEvent dLOG_MeshLoad;
PetscFList dMeshList;
static PetscBool dMeshPackageInitialized;

const char *const dMeshGenTypes[] = {"NONE","BLOCK","dMeshGenType","dMESHGEN_",0};

dErr dMeshCreate(MPI_Comm comm,dMesh *inm)
{
  dMesh m;
  dErr err;

  dFunctionBegin;
  dValidPointer(inm,2);
#if !defined(PETSC_USE_DYNAMIC_LIBRARIES)
  err = dMeshInitializePackage(PETSC_NULL);dCHK(err);
#endif
  *inm = 0;
  err = PetscHeaderCreate(m,_p_dMesh,struct _dMeshOps,dMESH_CLASSID,0,"dMesh",comm,dMeshDestroy,dMeshView);dCHK(err);
  *inm = m;
  dFunctionReturn(0);
}

dErr dMeshSetFromOptions(dMesh mesh)
{
  char type[256],str[dSTR_LEN];
  PetscMPIInt size;
  dBool flg;
  dMeshGenType gtype;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  err = MPI_Comm_size(((dObject)mesh)->comm,&size);dCHK(err);
  err = PetscStrncpy(type,size>1 ? dMESHPACK : dMESHSERIAL,sizeof(type));dCHK(err);
  err = PetscOptionsBegin(((PetscObject)mesh)->comm,((PetscObject)mesh)->prefix,"Mesh (dMesh) options","dMesh");dCHK(err);
  err = PetscOptionsList("-dmesh_type","Mesh type","dMeshSetType",dMeshList,type,type,sizeof(type),&flg);dCHK(err);
  if (flg || !((dObject)mesh)->type_name) {
    err = dMeshSetType(mesh,type);dCHK(err);
  }
  err = PetscOptionsString("-dmesh_in","File to read mesh from","dMeshSetInFile",mesh->infile,str,sizeof(str),&flg);dCHK(err);
  if (flg) {err = dMeshSetInFile(mesh,str,NULL);dCHK(err);}
  err = PetscOptionsString("-dmesh_inopts","Options for reading file","dMeshSetInFile",mesh->inoptions,str,sizeof(str),&flg);dCHK(err);
  if (flg) {err = dMeshSetInFile(mesh,NULL,str);dCHK(err);}
  err = PetscOptionsBool("-dmesh_intermediate_adjacencies","Infer face and edge adjacencies","dMeshInferIntermediateAdjacencies",flg=dFALSE,&flg,NULL);dCHK(err);
  if (flg) {err = dMeshInferIntermediateAdjacencies(mesh);dCHK(err);}
  err = PetscOptionsEnum("-dmesh_gen","Generate a mesh of a block at a resolution and size chosen by run-time parameters","dMeshSetGenerate",dMeshGenTypes,(PetscEnum)(gtype=dMESHGEN_NONE),(PetscEnum*)&gtype,&flg);dCHK(err);
  if (flg) {err = dMeshSetGenerate(mesh,gtype);dCHK(err);}
  if (mesh->ops->setfromoptions) {
    err = (*mesh->ops->setfromoptions)(mesh);dCHK(err);
  }
  err = PetscOptionsEnd();dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshSetType(dMesh mesh,const dMeshType type)
{
  dErr err,(*r)(dMesh);
  dBool     match;

  dFunctionBegin;
  PetscValidHeaderSpecific(mesh,dMESH_CLASSID,1);
  PetscValidCharPointer(type,2);
  err = PetscTypeCompare((PetscObject)mesh,type,&match);dCHK(err);
  if (match) dFunctionReturn(0);
  if (!dMeshRegisterAllCalled) {err = dMeshRegisterAll(NULL);dCHK(err);}
  err = PetscFListFind(dMeshList,((PetscObject)mesh)->comm,type,dTRUE,(void(**)(void))&r);dCHK(err);
  if (!r) dERROR(PETSC_COMM_SELF,1,"Unable to find requested dMesh type %s",type);
  if (mesh->ops->destroy) { err = (*mesh->ops->destroy)(mesh);dCHK(err); }
  err = (*r)(mesh);dCHK(err);
  err = PetscObjectChangeTypeName((PetscObject)mesh,type);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshRegister(const char name[],const char path[],const char cname[],dErr(*create)(dMesh))
{
  char fullname[dMAX_PATH_LEN];
  dErr err;

  dFunctionBegin;
  err = PetscFListConcat(path,cname,fullname);dCHK(err);
  err = PetscFListAdd(&dMeshList,name,fullname,(void (*)(void))create);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshInitializePackage(const char path[])
{
  dErr err;

  dFunctionBegin;
  if (dMeshPackageInitialized) dFunctionReturn(0);
  dMeshPackageInitialized = PETSC_TRUE;
  err = PetscClassIdRegister("Mesh",&dMESH_CLASSID);dCHK(err);
  err = dMeshRegisterAll(path);dCHK(err);
  err = PetscLogEventRegister("dMeshLoad", dMESH_CLASSID,&dLOG_MeshLoad);dCHK(err);
  err = PetscRegisterFinalize(dMeshFinalizePackage);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshFinalizePackage(void)
{
  dFunctionBegin;
  dMeshPackageInitialized = dFALSE;
  dMeshList = NULL;
  dMeshRegisterAllCalled = dFALSE;
  dFunctionReturn(0);
}

dEXTERN_C_BEGIN

extern dErr dMeshCreate_Pack(dMesh mesh);
extern dErr dMeshCreate_Serial(dMesh mesh);

dEXTERN_C_END

dErr dMeshRegisterAll(const char path[])
{
  dErr err;

  dFunctionBegin;
  if (dMeshRegisterAllCalled) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Function has already been called");
  err = dMeshRegisterDynamic(dMESHPACK    ,path,"dMeshCreate_Pack"     ,dMeshCreate_Pack);dCHK(err);
  err = dMeshRegisterDynamic(dMESHSERIAL  ,path,"dMeshCreate_Serial"   ,dMeshCreate_Serial);dCHK(err);
  dMeshRegisterAllCalled = PETSC_TRUE;
  dFunctionReturn(0);
}
