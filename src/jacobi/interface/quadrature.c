#include <dohpjacimpl.h>
#include <dohp.h>

PetscFList dQuadratureList;
PetscBool dQuadratureRegisterAllCalled;

const char *const dQuadratureMethods[] = {"FAST","SPARSE","SELF","dQuadratureMethod","dQUADRATURE_METHOD_",0};

/** Create an array of rules for integrating functions of given order on the reference element
*
* @param quad the context
* @param n number of elements
* @param topo topology of the element
* @param order order of polynomial to be integrated exactly
* @param rules pointer to new array of rules
*/
dErr dQuadratureGetRules(dQuadrature quad,dInt n,const dEntTopology topo[],const dPolynomialOrder order[],dRule **rules)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(quad,dQUADRATURE_CLASSID,1);
  dValidPointer(topo,3);
  dValidPointer(order,4);
  dValidPointer(rules,5);
  err = dMallocA(n,rules);dCHK(err);
  err = quad->ops->GetRule(quad,n,topo,order,*rules);dCHK(err);
  dFunctionReturn(0);
}

/** Create an array of rules for integrating functions of given order on a facet of the reference element
*
* @param quad the context
* @param n number of elements
* @param topo topology of the element relative to which we will integrate
* @param facet index of lower-dimensional entity upon which to integrate
* @param order order of polynomial to be integrated exactly (in the reference frame of the element)
* @param rules pointer to new array of rules
*/
dErr dQuadratureGetFacetRules(dQuadrature quad,dInt n,const dEntTopology topo[],const dInt facet[],const dPolynomialOrder order[],dRule **rules)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(quad,dQUADRATURE_CLASSID,1);
  dValidPointer(topo,3);
  dValidPointer(facet,4);
  dValidPointer(order,5);
  dValidPointer(rules,6);
  err = dMallocA(n,rules);dCHK(err);
  err = quad->ops->GetFacetRule(quad,n,topo,facet,order,*rules);dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureSetFromOptions(dQuadrature quad)
{
  char type[dNAME_LEN] = dQUADRATURE_TENSOR;
  dBool typeSet;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(quad,dQUADRATURE_CLASSID,1);
  err = PetscOptionsBegin(((PetscObject)quad)->comm,((PetscObject)quad)->prefix,"Quadrature options","dQuadrature");dCHK(err);
  err = PetscOptionsList("-dquad_type","Quadrature type","dQuadratureSetType",dQuadratureList,
                          (((PetscObject)quad)->type_name?((PetscObject)quad)->type_name:type),type,dNAME_LEN,&typeSet);dCHK(err);
  if (typeSet || !((PetscObject)quad)->type_name) {
    err = dQuadratureSetType(quad,type);dCHK(err);
  }
  if (quad->ops->SetFromOptions) {
    err = quad->ops->SetFromOptions(quad);dCHK(err);
  }
  err = PetscOptionsEnd();dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureSetMethod(dQuadrature quad,dQuadratureMethod method)
{
  dErr err;

  dFunctionBegin;
  err = (*quad->ops->SetMethod)(quad,method);dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureDestroy(dQuadrature *quad)
{
  dErr err;

  dFunctionBegin;
  if (!*quad) dFunctionReturn(0);
  dValidHeader(*quad,dQUADRATURE_CLASSID,1);
  if (--((PetscObject)*quad)->refct > 0) dFunctionReturn(0);
  if ((*quad)->ops->Destroy) {
    err = (*quad)->ops->Destroy(*quad);dCHK(err);
  }
  err = PetscHeaderDestroy(quad);dCHK(err);
  dFunctionReturn(0);
}

dErr dQuadratureView(dQuadrature quad,PetscViewer viewer)
{
  dBool  iascii;
  dErr   err;

  dFunctionBegin;
  PetscValidHeaderSpecific(quad,dQUADRATURE_CLASSID,1);
  if (!viewer) {
    err = PetscViewerASCIIGetStdout(((PetscObject)quad)->comm,&viewer);dCHK(err);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(quad,1,viewer,2);

  err = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);dCHK(err);
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
  dBool  match;

  dFunctionBegin;
  PetscValidHeaderSpecific(quad,dQUADRATURE_CLASSID,1);
  PetscValidCharPointer(type,2);
  err = PetscObjectTypeCompare((PetscObject)quad,type,&match);dCHK(err);
  if (match) dFunctionReturn(0);
  if (!dQuadratureRegisterAllCalled) {err = dQuadratureRegisterAll(NULL);dCHK(err);}
  err = PetscFListFind(dQuadratureList,((PetscObject)quad)->comm,type,dTRUE,(void(**)(void))&r);dCHK(err);
  if (!r) dERROR(PETSC_COMM_SELF,1,"Unable to find requested dQuadrature type %s",type);
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
  dErr err;

  dFunctionBegin;
  if (dQuadratureRegisterAllCalled) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Function has already been called");
  err = dQuadratureRegisterDynamic(dQUADRATURE_TENSOR,path,"dQuadratureCreate_Tensor",dQuadratureCreate_Tensor);dCHK(err);
  dQuadratureRegisterAllCalled = PETSC_TRUE;
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
  err = PetscHeaderCreate(quad,p_dQuadrature,struct _dQuadratureOps,dQUADRATURE_CLASSID,0,"dQuadrature","Quadrature rule service","Jacobi",comm,dQuadratureDestroy,dQuadratureView);dCHK(err);

  *inquad = quad;
  dFunctionReturn(0);
}
