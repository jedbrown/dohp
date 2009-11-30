#include <MBTagConventions.hpp>
#include <dohpfsimpl.h>

PetscLogEvent dLOG_Q1HexComputeQuadrature,dLOG_FSMatSetValuesExpanded;
PetscCookie dFSROT_COOKIE;
static PetscFList FSList = 0;

/**
* These default operations are shared with the DM.  We are making a two-level inheritance since there may be different
* dFS implementations.  Certainly both continuous and discontinuous Galerkin are desirable.
*
* Since we use VecGhost, the local and global vectors are really the same.  Since there is no way to get a global vector
* from a local vector, the user cannot avoid seeing VecGhost.  We do have expanded vectors which are genuinely distinct,
* but I'm very skeptical of the usefulness of identifying them with the "local vector" of a DM.
**/
static const struct _dFSOps defaultFSOps = { .view = dFSView,
                                             .createglobalvector = dFSCreateGlobalVector,
                                             .createlocalvector  = dFSCreateExpandedVector,
                                             /* I think that these don't make sense with the current design.  In
                                             * particular, there may be points in the expanded space which are not
                                             * represented in the global vector.  In this configuration, an INSERT_MODE
                                             * scatter is not sufficient because it is not surjective.
                                             *
                                             .globaltolocalbegin = dFSGlobalToExpandedBegin,
                                             .globaltolocalend   = dFSGlobalToExpandedEnd,
                                             */
                                             .getmatrix          = dFSGetMatrix,
                                             .destroy            = dFSDestroy,
                                             .impldestroy        = 0};

dErr dFSCreate(MPI_Comm comm,dFS *infs)
{
  dFS fs;
  dErr err;

  dFunctionBegin;
  dValidPointer(infs,2);
  *infs = 0;
#if !defined(PETSC_USE_DYNAMIC_LIBRARIES)
  err = dFSInitializePackage(PETSC_NULL);dCHK(err);
#endif
  err = PetscHeaderCreate(fs,_p_dFS,struct _dFSOps,DM_COOKIE,0,"dFS",comm,dFSDestroy,dFSView);dCHK(err);

  err = dMemcpy(fs->ops,&defaultFSOps,sizeof defaultFSOps);dCHK(err);
  err = dStrcpyS(fs->bdyTagName,sizeof fs->bdyTagName,NEUMANN_SET_TAG_NAME);dCHK(err);
  /* RCM is a good default ordering because it improves the performance of smoothers and incomplete factorization */
  err = dStrcpyS(fs->orderingtype,sizeof fs->orderingtype,MATORDERING_RCM);dCHK(err);

  /* Defaults */
  fs->bs = 1;
  err = dCallocA(1,&fs->fieldname);dCHK(err);

  *infs = fs;
  dFunctionReturn(0);
}

dErr dFSSetType(dFS fs,const dFSType type)
{
  dErr err,(*r)(dFS);
  dBool     match;

  dFunctionBegin;
  PetscValidHeaderSpecific(fs,DM_COOKIE,1);
  PetscValidCharPointer(type,2);
  err = PetscTypeCompare((PetscObject)fs,type,&match);dCHK(err);
  if (match) dFunctionReturn(0);
  err = PetscFListFind(FSList,((PetscObject)fs)->comm,type,(void(**)(void))&r);dCHK(err);
  if (!r) dERROR(1,"Unable to find requested dFS type %s",type);
  if (fs->ops->impldestroy) { err = (*fs->ops->impldestroy)(fs);dCHK(err); }
  err = PetscMemcpy(fs->ops,&defaultFSOps,sizeof(defaultFSOps));dCHK(err);
  fs->spacebuilt = 0;
  err = (*r)(fs);dCHK(err);
  err = PetscObjectChangeTypeName((PetscObject)fs,type);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSSetFromOptions(dFS fs)
{
  static const char deft[] = dFSCONT;
  char type[256];
  dBool flg;
  dErr err;

  dFunctionBegin;
  dValidHeader(fs,DM_COOKIE,1);
  err = PetscOptionsBegin(((PetscObject)fs)->comm,((PetscObject)fs)->prefix,"Function Space (dFS) options","dFS");dCHK(err);
  err = PetscOptionsList("-dfs_type","Function Space type","dFSSetType",FSList,deft,type,256,&flg);dCHK(err);
  if (flg) {
    err = dFSSetType(fs,type);dCHK(err);
  } else if (!((dObject)fs)->type_name) {
    err = dFSSetType(fs,deft);dCHK(err);
  }
  err = PetscOptionsInt("-dfs_rule_strength","Choose rules that are stronger than necessary","dFSRuleStrength",fs->ruleStrength,&fs->ruleStrength,NULL);dCHK(err);
  err = PetscOptionsList("-dfs_ordering_type","Function Space ordering, usually to reduce bandwidth","dFSBuildSpace",MatOrderingList,fs->orderingtype,fs->orderingtype,256,NULL);dCHK(err);
  err = PetscOptionsTruth("-dfs_assemble_reduced","Assemble only diagonal part of blocks","",fs->assemblereduced,&fs->assemblereduced,NULL);dCHK(err);
  if (fs->ops->setfromoptions) {
    err = (*fs->ops->setfromoptions)(fs);dCHK(err);
  }
  err = PetscOptionsEnd();dCHK(err);
  err = dFSBuildSpace(fs);dCHK(err);

  /* View if requested */
  err = PetscOptionsHasName(((dObject)fs)->prefix,"-dfs_view",&flg);dCHK(err);
  if (flg) {
    dViewer viewer;
    err = PetscViewerASCIIGetStdout(((dObject)fs)->comm,&viewer);dCHK(err);
    err = dFSView(fs,viewer);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dFSRegister(const char name[],const char path[],const char cname[],dErr(*create)(dFS))
{
  char fullname[dMAX_PATH_LEN];
  dErr err;

  dFunctionBegin;
  err = PetscFListConcat(path,cname,fullname);dCHK(err);
  err = PetscFListAdd(&FSList,name,fullname,(void (*)(void))create);dCHK(err);
  dFunctionReturn(0);
}

dErr dFSInitializePackage(const char path[])
{
  static dBool initialized = PETSC_FALSE;
  dErr err;

  dFunctionBegin;
  if (initialized) dFunctionReturn(0);
  initialized = PETSC_TRUE;
  err = MatInitializePackage(path);dCHK(err);
  err = DMInitializePackage(path);dCHK(err);
  err = PetscCookieRegister("dFS Rotation",&dFSROT_COOKIE);dCHK(err);

  err = PetscLogEventRegister("dQ1HexComputQuad",DM_COOKIE,&dLOG_Q1HexComputeQuadrature);dCHK(err); /* only vaguely related */
  err = PetscLogEventRegister("dFSMatSetVExpand",DM_COOKIE,&dLOG_FSMatSetValuesExpanded);dCHK(err);
  err = dFSRegisterAll(path);dCHK(err);
  dFunctionReturn(0);
}



dEXTERN_C_BEGIN

extern dErr dFSCreate_Cont(dFS);

dEXTERN_C_END

dErr dFSRegisterAll(const char path[])
{
  static dBool called = PETSC_FALSE;
  dErr err;

  dFunctionBegin;
  err = dFSRegisterDynamic(dFSCONT,path,"dFSCreate_Cont",dFSCreate_Cont);dCHK(err);
  called = PETSC_TRUE;
  dFunctionReturn(0);
}
