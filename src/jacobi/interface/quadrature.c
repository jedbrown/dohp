#include <dohpjacimpl.h>
#include <dohp.h>

static PetscFList dQuadratureList;

/** Fill an array of dRule starting at \a firstrule.
*
* @param quad the context
* @param n number of elements
* @param topo topology of the element
* @param rsize number of points in each Cartesian direction
* @param firstrule place to put the newly constructed dRule, normally this will be \c s_dRule[]
*/
dErr dQuadratureGetRule(dQuadrature quad,dInt n,const dEntTopology topo[],const dInt rsize[],dRule firstrule)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(quad,dQUADRATURE_COOKIE,1);
  dValidPointer(topo,3);
  dValidPointer(rsize,4);
  dValidPointer(firstrule,5);
  err = quad->ops->GetRule(quad,n,topo,rsize,firstrule);dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureSetFromOptions(dQuadrature quad)
{
  char type[dNAME_LEN] = dQUADRATURE_TENSOR;
  dBool typeSet;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(quad,dQUADRATURE_COOKIE,1);
  err = PetscOptionsBegin(((PetscObject)quad)->comm,((PetscObject)quad)->prefix,"Quadrature options","dQuadrature");dCHK(err);
  err = PetscOptionsList("-dquad_type","Quadrature type","dQuadratureSetType",dQuadratureList,
                          (((PetscObject)quad)->type_name?((PetscObject)quad)->type_name:type),type,dNAME_LEN,&typeSet);dCHK(err);
  if (typeSet) {
    err = dQuadratureSetType(quad,type);dCHK(err);
  }
  if (!((PetscObject)quad)->type_name) {
    err = dQuadratureSetType(quad,type);dCHK(err);
  }
  if (quad->ops->SetFromOptions) {
    err = quad->ops->SetFromOptions(quad);dCHK(err);
  }
  err = PetscOptionsEnd();dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureDestroy(dQuadrature quad)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(quad,dQUADRATURE_COOKIE,1);
  if (--((PetscObject)quad)->refct > 0) dFunctionReturn(0);
  if (quad->ops->Destroy) {
    err = quad->ops->Destroy(quad);dCHK(err);
  }
  err = PetscHeaderDestroy(quad);dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureView(dQuadrature quad,PetscViewer viewer)
{
  dTruth iascii;
  dErr   err;

  dFunctionBegin;
  PetscValidHeaderSpecific(quad,dQUADRATURE_COOKIE,1);
  if (!viewer) {
    err = PetscViewerASCIIGetStdout(((PetscObject)quad)->comm,&viewer);dCHK(err);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_COOKIE,2);
  PetscCheckSameComm(quad,1,viewer,2);

  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);dCHK(err);
  if (iascii) {
    err = PetscViewerASCIIPrintf(viewer,"dQuadobi object:(%s)\n",
                                  ((PetscObject)quad)->prefix ? ((PetscObject)quad)->prefix : "no prefix");dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"type: %s\n",
                                  ((PetscObject)quad)->type_name ? ((PetscObject)quad)->type_name : "type not set");dCHK(err);
    if (quad->ops->View) {
      err = (*quad->ops->View)(quad,viewer);dCHK(err);
    } else {
      err = PetscViewerASCIIPrintf(viewer,"Internal info not available.\n");dCHK(err);
    }
    err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  } else if (quad->ops->View) {
    err = (*quad->ops->View)(quad,viewer);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dQuadratureSetType(dQuadrature quad,dQuadratureType type)
{
  dErr   err,(*r)(dQuadrature);
  dTruth match;

  dFunctionBegin;
  PetscValidHeaderSpecific(quad,dQUADRATURE_COOKIE,1);
  PetscValidCharPointer(type,2);
  err = PetscTypeCompare((PetscObject)quad,type,&match);dCHK(err);
  if (match) dFunctionReturn(0);
  err = PetscFListFind(dQuadratureList,((PetscObject)quad)->comm,type,(void(**)(void))&r);dCHK(err);
  if (!r) dERROR(1,"Unable to find requested dQuadrature type %s",type);
  if (quad->ops->Destroy) { err = (*quad->ops->Destroy)(quad);dCHK(err); }
  err = PetscMemzero(quad->ops,sizeof(quad->ops));dCHK(err);
  err = (*r)(quad);dCHK(err);
  err = PetscObjectChangeTypeName((PetscObject)quad,type);dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureRegister(const char name[],const char path[],const char cname[],dErr(*create)(dQuadrature))
{
  char fullname[dMAX_PATH_LEN];
  dErr err;

  dFunctionBegin;
  err = PetscFListConcat(path,cname,fullname);dCHK(err);
  err = PetscFListAdd(&dQuadratureList,name,fullname,(void (*)(void))create);dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureRegisterAll(const char path[])
{
  static dBool called = PETSC_FALSE;
  dErr err;

  dFunctionBegin;
  if (called) dFunctionReturn(0);
  err = dQuadratureRegisterDynamic(dQUADRATURE_TENSOR,path,"dQuadratureCreate_Tensor",dQuadratureCreate_Tensor);dCHK(err);
  called = PETSC_TRUE;
  dFunctionReturn(0);
}

dErr dQuadratureCreate(MPI_Comm comm,dQuadrature *inquad)
{
  dErr err;
  dQuadrature quad;

  dFunctionBegin;
  dValidPointer(inquad,2);
  *inquad = 0;
#if !defined PETSC_USE_DYNAMIC_LIBRARIES
  err = dJacobiInitializePackage(NULL);dCHK(err);
#endif
  err = PetscHeaderCreate(quad,p_dQuadrature,struct _dQuadratureOps,dQUADRATURE_COOKIE,0,"dQuadrature",comm,dQuadratureDestroy,dQuadratureView);dCHK(err);

  *inquad = quad;
  dFunctionReturn(0);
}