/**
* @file   jacobi.c
* @author Jed Brown <jed@59A2.org>
* @date   Fri Aug 22 14:06:41 2008
* 
* @brief  Compute quadrature rules, interpolation, and differentiation matrices
* 
* 
*/

#include "petsc.h"
#include "dohpjacobi.h"

PetscCookie DOHP_JACOBI_COOKIE;
static PetscFList DohpJacobiList = 0;

static const struct _DohpJacobiOps _defaultOps = {
  .view = 0,
  .setup = 0,
  .setfromoptions = 0,
  .destroy = 0
};

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiCreate"
/** 
* Create a new Jacobi object and initialize with defaults.
* 
* @param comm 
* @param injacobi 
* 
* @return 
*/
PetscErrorCode DohpJacobiCreate(MPI_Comm comm,DohpJacobi *injacobi)
{
  DohpJacobi jac;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidPointer(injacobi,2);
  *injacobi = 0;
#if !defined(PETSC_USE_DYNAMIC_LIBRARIES)
  ierr = DohpJacobiInitializePackage(PETSC_NULL);CHKERRQ(ierr);
#endif
  ierr = PetscHeaderCreate(jac,_p_DohpJacobi,struct _DohpJacobiOps,DOHP_JACOBI_COOKIE,0,"DohpJacobi",comm,DohpJacobiDestroy,DohpJacobiView);CHKERRQ(ierr);

  jac->basisdegree = 10;
  jac->ruleexcess = 5;
  jac->setupcalled = 0;
  jac->impl = 0;
  ierr = PetscMemcpy(jac->ops,&_defaultOps,sizeof(struct _DohpJacobiOps));CHKERRQ(ierr);

  *injacobi = jac;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiSetType"
/** 
* Set the type for a DohpJacobi object.
* 
* @param jac 
* @param type 
* 
* @return 
*/
PetscErrorCode DohpJacobiSetType(DohpJacobi jac,DohpJacobiType type)
{
  PetscErrorCode ierr,(*r)(DohpJacobi);
  PetscTruth     match;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(jac,DOHP_JACOBI_COOKIE,1);
  PetscValidCharPointer(type,2);
  ierr = PetscTypeCompare((PetscObject)jac,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);
  ierr = PetscFListFind(DohpJacobiList,((PetscObject)jac)->comm,type,(void(**)(void))&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(1,"Unable to find requested DohpJacobi type %s",type);
  if (jac->ops->destroy) { ierr = (*jac->ops->destroy)(jac);CHKERRQ(ierr); }
  ierr = PetscMemcpy(jac->ops,&_defaultOps,sizeof(struct _DohpJacobiOps));CHKERRQ(ierr);
  jac->setupcalled = 0;
  ierr = (*r)(jac);CHKERRQ(ierr);
  ierr = PetscObjectChangeTypeName((PetscObject)jac,type);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiSetFromOptions"
/** 
* Set options from the options database.
* 
* @param jac
* 
* @return 
*/
PetscErrorCode DohpJacobiSetFromOptions(DohpJacobi jac)
{
  char type[256] = DOHP_JACOBI_LGL;
  PetscTruth typeSet;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(jac,DOHP_JACOBI_COOKIE,1);
  ierr = PetscOptionsBegin(((PetscObject)jac)->comm,((PetscObject)jac)->prefix,"Jacobi options (type and size of basis/quadrature rules)","DohpJacobi");CHKERRQ(ierr);
  ierr = PetscOptionsList("-djac_type","Basis/Quadrature type","DohpJacobiSetType",DohpJacobiList,
                          (((PetscObject)jac)->type_name?((PetscObject)jac)->type_name:type),type,256,&typeSet);CHKERRQ(ierr);
  if (typeSet) {
    ierr = DohpJacobiSetType(jac,type);CHKERRQ(ierr);
  }
  if (!((PetscObject)jac)->type_name) {
    ierr = DohpJacobiSetType(jac,type);CHKERRQ(ierr);
  }
  ierr = PetscOptionsInt("-djac_basis_degree","Max basis degree","DohpJacobiSetDegrees",jac->basisdegree,&jac->basisdegree,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-djac_rule_excess","Excess quadrature points","DohpJacobiSetDegrees",jac->ruleexcess,&jac->ruleexcess,PETSC_NULL);CHKERRQ(ierr);
  if (jac->ops->setfromoptions) {
    ierr = jac->ops->setfromoptions(jac);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiSetUp"
/** 
* Initialize the Jacobi object.
* 
* @param jac 
* 
* @return 
*/
PetscErrorCode DohpJacobiSetUp(DohpJacobi jac)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(jac,DOHP_JACOBI_COOKIE,1);
  if (!jac->setupcalled && jac->ops->setup) {
    ierr = jac->ops->setup(jac);CHKERRQ(ierr);
  }
  jac->setupcalled = 1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiDestroy"
/** 
* Destroy a Jacobi object.
* 
* @param jac 
* 
* @return 
*/
PetscErrorCode DohpJacobiDestroy(DohpJacobi jac)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(jac,DOHP_JACOBI_COOKIE,1);
  if (jac->ops->destroy) {
    ierr = jac->ops->destroy(jac);CHKERRQ(ierr);
  }
  ierr = PetscHeaderDestroy(jac);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiView"
/** 
* View the state of a DohpJacobi.
* 
* @param jac 
* @param viewer 
* 
* @return 
*/
PetscErrorCode DohpJacobiView(DohpJacobi jac,PetscViewer viewer)
{
  PetscTruth iascii;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(jac,DOHP_JACOBI_COOKIE,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(((PetscObject)jac)->comm,&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_COOKIE,2);
  PetscCheckSameComm(jac,1,viewer,2);

  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"DohpJacobi object:(%s)\n",
                                  ((PetscObject)jac)->prefix ? ((PetscObject)jac)->prefix : "no prefix");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"type: %s\n",
                                  ((PetscObject)jac)->type_name ? ((PetscObject)jac)->type_name : "type not set");CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"max basis degree: %d\n",jac->basisdegree);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"rule excess: %d\n",jac->basisdegree);CHKERRQ(ierr);
    if (!jac->setupcalled) {
      ierr = PetscViewerASCIIPrintf(viewer,"Object has not been set up.\n",jac->basisdegree);CHKERRQ(ierr);
    }
    if (jac->ops->view) {
      ierr = (*jac->ops->view)(jac,viewer);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"Internal info not available.\n");CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  } else if (jac->ops->view) {
    ierr = (*jac->ops->view)(jac,viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiRegister"
PetscErrorCode DohpJacobiRegister(const char name[],const char path[],const char cname[],PetscErrorCode(*create)(DohpJacobi))
{
  char fullname[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFListConcat(path,cname,fullname);CHKERRQ(ierr);
  ierr = PetscFListAdd(&DohpJacobiList,name,fullname,(void (*)(void))create);CHKERRQ(ierr); 
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DohpJacobiRegisterAll"
PetscErrorCode DohpJacobiRegisterAll(const char path[])
{
  static PetscTruth called = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DohpJacobiRegisterDynamic(DOHP_JACOBI_LGL,path,"DohpJacobiCreate_LGL",DohpJacobiCreate_LGL);CHKERRQ(ierr);
  called = PETSC_TRUE;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiInitializePackage"
PetscErrorCode DohpJacobiInitializePackage(const char path[])
{
  static PetscTruth initialized = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (initialized) PetscFunctionReturn(0);
  ierr = PetscCookieRegister("Jacobi context",&DOHP_JACOBI_COOKIE);CHKERRQ(ierr);
  ierr = DohpJacobiRegisterAll(path);CHKERRQ(ierr);
  initialized = PETSC_TRUE;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DohpJacobiSetDegrees"
/** 
* Set the maximum size of the approximation space generated by Jacobi.
*
* Jacobi will always generate quadrature rules up to the maximum order.  This is to save us from a degenerate case where
* one field has low order on an element but another has very high order.  In this case, a quadrature order close to \p
* basisdegree + \p ruleexcess will be required due to the second field.  We will not normally generate quadrature rules
* with fewer points than the number of functions in the basis because this makes the element mass matrix singular.
*
* @param jac The Jacobi context
* @param basisdegree The maximum number of functions in a 1D basis.
* @param ruleexcess The number of extra quadrature points to generate rules for.
* 
* @return 
*/
PetscErrorCode DohpJacobiSetDegrees(DohpJacobi jac,PetscInt basisdegree,PetscInt ruleexcess)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(jac,DOHP_JACOBI_COOKIE,1);
  if (jac->setupcalled) {
    if (jac->ops->destroy) { ierr = (*jac->ops->destroy)(jac);CHKERRQ(ierr); }
    ierr = PetscMemcpy(jac->ops,&_defaultOps,sizeof(struct _DohpJacobiOps));CHKERRQ(ierr);
  }
  jac->basisdegree = basisdegree;
  jac->ruleexcess = ruleexcess;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpJacobiGetRule"
/** 
* Get a pointer to the rule.
* 
* @param[in] jac context
* @param[in] n number of quadrature nodes in requested rule
* @param[out] rule rule handle
* 
* @return 
*/
PetscErrorCode DohpJacobiGetRule(DohpJacobi jac,PetscInt n,DohpRule *rule)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(jac,DOHP_JACOBI_COOKIE,1);
  ierr = jac->ops->getrule(jac,n,rule);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DohpJacobiGetBasis"
/** 
* Get a basis context.
* 
* @param jac 
* @param bsize number of basis functions
* @param qsize number of quadrature points
* @param[out] basis basis handle
* 
* @return 
*/
PetscErrorCode DohpJacobiGetBasis(DohpJacobi jac,PetscInt bsize,PetscInt qsize,DohpBasis *basis)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(jac,DOHP_JACOBI_COOKIE,1);
  ierr = jac->ops->getbasis(jac,bsize,qsize,basis);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
