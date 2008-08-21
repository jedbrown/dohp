#include "private/dohpimpl.h"

static const struct _DohpQuotientOps DohpQuotientDefaultOps = {
  .update = 0,
  .setup = 0,
  .setfromoptions = 0,
  .destroy = 0
};

PetscCookie DOHP_QUOTIENT_COOKIE,DOHP_MESH_COOKIE;
static PetscFList DohpQuotientList = 0;

#undef __FUNCT__
#define __FUNCT__ "DohpQuotientUpdate"
/*@
   DohpQuotientUpdate - update the quadrature order and element maps

   Collective on DohpQuotient

   Input parameter:
.  q - the quotient
@*/
PetscErrorCode DohpQuotientUpdate(DohpQuotient q)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(q,DOHP_QUOTIENT_COOKIE,1);
  ierr = (*q->ops->update)(q);CHKERRQ(ierr);
  PetscObjectStateIncrease((PetscObject)q);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpQuotientCreate"
/*@
   DohpQuotientCreate - set up a quotient map over the mesh with size given by the values in the tag

   Collective on DohpMesh

   Input Parameters:
.  m - the mesh over which to set up
.  loc - an entity set handle corresponding to the local portion of the domain
.  qsizetag - a tag handle on the mesh, valid at least for all regions in loc

   Output Parameters:
.  q - the new Quotient
@*/
PetscErrorCode DohpQuotientCreate(DohpMesh m,DohpESH loc,DohpTag qsizetag,DohpQuotient *inq)
{
  DohpQuotient   q;
  MPI_Comm       comm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(m,DOHP_MESH_COOKIE,1);
  PetscValidPointer(inq,4);
  *inq = 0;
  ierr = PetscObjectGetComm((PetscObject)m,&comm);CHKERRQ(ierr);
#if !defined(PETSC_USE_DYNAMIC_LIBRARIES)
  ierr = DohpQuotientInitializePackage(PETSC_NULL);CHKERRQ(ierr);
#endif
  ierr = PetscHeaderCreate(q,_p_DohpQuotient,struct _DohpQuotientOps,DOHP_QUOTIENT_COOKIE,0,"DohpQuotient",comm,DohpQuotientDestroy,DohpQuotientView);CHKERRQ(ierr);

  q->mesh        = m;
  q->loc         = loc;
  q->qsizetag    = qsizetag;
  q->ops->update = 0;

  *inq = q;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpQuotientSetType"
/*@
   DohpQuotientSetType - 

@*/
PetscErrorCode DohpQuotientSetType(DohpQuotient q,const DohpQuotientType type)
{
  PetscErrorCode (*r)(DohpQuotient);
  PetscTruth     match;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(q,DOHP_QUOTIENT_COOKIE,1);
  PetscValidCharPointer(type,2);
  ierr = PetscTypeCompare((PetscObject)q,type,&match);CHKERRQ(ierr);
  if (match) PetscFunctionReturn(0);
  ierr =  PetscFListFind(DohpQuotientList,((PetscObject)q)->comm,type,(void (**)(void)) &r);CHKERRQ(ierr);
  if (!r) SETERRQ1(PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested DohpQuotient type %s",type);
  if (q->ops->destroy) { ierr = (*q->ops->destroy)(q);CHKERRQ(ierr); }
  ierr = PetscMemcpy(q->ops,&DohpQuotientDefaultOps,sizeof(struct _DohpQuotientOps));CHKERRQ(ierr);
  q->setupcalled = 0;
  ierr = (*r)(q);CHKERRQ(ierr);
  ierr = PetscObjectChangeTypeName((PetscObject)q,type);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpQuotientSetFromOptions"
/*@
   DohpQuotientSetFromOptions - 

@*/
PetscErrorCode DohpQuotientSetFromOptions(DohpQuotient q)
{
  PetscTruth flg;
  const DohpQuotientType deft = DohpQuotientGauss;
  char type[256];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(q,DOHP_QUOTIENT_COOKIE,1);
  ierr = PetscOptionsBegin(((PetscObject)q)->comm,((PetscObject)q)->prefix,"Quotient (quadrature rule and element map) options","DohpQuotient");CHKERRQ(ierr);
  ierr = DohpQuotientRegisterAll(PETSC_NULL);CHKERRQ(ierr);
  if (((PetscObject)q)->type_name) { deft = ((PetscObject)q)->type_name; }
  ierr = PetscOptionsList("-dquot_type","Quotient type","DohpQuotientSetType",DohpQuotientList,deft,type,256,&flg);CHKERRQ(ierr);
  ierr = DohpQuotientSetType(q,type);CHKERRQ(ierr);
  if (q->ops->setfromoptions) {
    ierr = (*q->ops->setfromoptions)(q);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DohpQuotientSetUp"
/*@
   DohpQuotientSetUp - 

@*/
PetscErrorCode DohpQuotientSetUp(DohpQuotient q)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(q,DOHP_QUOTIENT_COOKIE,1);
  if (!q->setupcalled) {
    ierr = (*q->ops->setup)(q);CHKERRQ(ierr);
  }
  q->setupcalled = 1;
  PetscFunctionReturn(0);
}


#define DohpQuotientRegisterDynamic(a,b,c,d) DohpQuotientRegister(a,b,c,d)

#undef __FUNCT__
#define __FUNCT__ "DohpQuotientRegister"
/*@
   DohpQuotientRegister - 

@*/
PetscErrorCode DohpQuotientRegister(const char sname[],const char path[],const char name[],PetscErrorCode (*function)(DohpQuotient))
{
  char           fullname[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFListConcat(path,name,fullname);CHKERRQ(ierr);
  ierr = PetscFListAdd(&DohpQuotientList,sname,fullname,(void (*)(void))function);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpQuotientRegisterAll"
/*@
   DohpQuotientRegisterAll - 

@*/
PetscErrorCode DohpQuotientRegisterAll(const char path[])
{
  static PetscTruth registered = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (registered) PetscFunctionReturn(0);
  ierr = DohpQuotientRegisterDynamic(DohpQuotientGauss,path,"DohpQuotientCreate_Gauss",DohpQuotientCreate_Gauss);CHKERRQ(ierr);
  registered = PETSC_TRUE;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpQuotientDestroy"
/*@
   DohpQuotientDestroy - 

@*/
PetscErrorCode DohpQuotientDestroy(DohpQuotient q)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(q,DOHP_QUOTIENT_COOKIE,1);
  if (q->ops->destroy) {
    ierr = (*q->ops->destroy)(q);CHKERRQ(ierr);
  }
  ierr = PetscFree(q->quad);CHKERRQ(ierr);
  ierr = PetscFree(q->map);CHKERRQ(ierr);
  ierr = PetscHeaderDestroy(q);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpQuotientGetType"
/*@
   DohpQuotientGetType - 

@*/
PetscErrorCode DohpQuotientGetType(DohpQuotient q,const DohpQuotientType *type)
{

  PetscFunctionBegin;
  PetscValidHeaderSpecific(q,DOHP_QUOTIENT_COOKIE,1);
  PetscValidPointer(type,2);
  *type = ((PetscObject)q)->type_name;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DohpQuotientView"
/*@
   DohpQuotientView - 

@*/
PetscErrorCode DohpQuotientView(DohpQuotient q,PetscViewer viewer)
{
  const DohpQuotientType  type;
  char                   *tagname;
  PetscTruth              iascii;
  PetscErrorCode          ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(q,DOHP_QUOTIENT_COOKIE,1);
  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(((PetscObject)q)->comm,&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_COOKIE,2);
  PetscCheckSameComm(q,1,viewer,2);
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = DohpQuotientGetType(q,&type);CHKERRQ(ierr);
    if (((PetscObject)q)->prefix) {
      ierr = PetscViewerASCIIPrintf(viewer,"DohpQuotient object:(%s)\n",((PetscObject)q)->prefix);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"DohpQuotient object:\n");CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"type: %s\n",(type ? type : "not yet set"));CHKERRQ(ierr);
    if (q->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*q->ops->view)(q,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"numbor of elements: %d\n",q->nelems);CHKERRQ(ierr);
    if (q->qsizetag) {
      ierr = DohpMeshGetTagName(q->mesh,q->qsizetag,&tagname);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"quadrature size tag: %s\n",tagname);CHKERRQ(ierr);
    }
  } else {
    if (q->ops->view) {
      ierr = (*q->ops->view)(q,viewer);CHKERRQ(ierr);
    }
  }
  ierr = DohpMeshView(q->mesh,viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpQuotientInitializePackage"
/*@
   DohpQuotientInitializePackage - 

@*/
PetscErrorCode DohpQuotientInitializePackage(const char path[])
{
  static PetscTruth initialized = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (initialized) PetscFunctionReturn(0);
  ierr = PetscCookieRegister("Quotient map",&DOHP_QUOTIENT_COOKIE);CHKERRQ(ierr);
  ierr = PetscCookieRegister("Mesh",&DOHP_MESH_COOKIE);CHKERRQ(ierr);
  ierr = DohpMeshRegisterAll(path);CHKERRQ(ierr);
  ierr = DohpQuotientRegisterAll(path);CHKERRQ(ierr);
  initialized = PETSC_TRUE;
  PetscFunctionReturn(0);
}
