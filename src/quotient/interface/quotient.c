#include <dohpquotientimpl.h>
#include <dohp.h>

static const struct _dQuotientOps dQuotientDefaultOps = {
  .update = 0,
  .setup = 0,
  .setfromoptions = 0,
  .destroy = 0
};

dCookie dQUOTIENT_COOKIE;
static PetscFList dQuotientList = 0;

static dErr dQuotientSetUp_Private(dQuotient q);

#undef __FUNCT__
#define __FUNCT__ "dQuotientUpdate"
/*@
   dQuotientUpdate - update the quadrature order and element maps

   Collective on dQuotient

   Input parameter:
.  q - the quotient
@*/
dErr dQuotientUpdate(dQuotient q)
{
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(q,dQUOTIENT_COOKIE,1);
  if (q->ops->update) {
    err = (*q->ops->update)(q);dCHK(err);
  } else {
    dERROR(1,"No update function given");
  }
  PetscObjectStateIncrease((PetscObject)q);
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientCreate"
/*@
   dQuotientCreate - set up a quotient map over the mesh with size given by the values in the tag

   Collective on dMesh

   Input Parameters:
.  m - the mesh over which to set up
.  loc - an entity set handle corresponding to the local portion of the domain
.  qsizetag - a tag handle on the mesh, valid at least for all regions in loc

   Output Parameters:
.  q - the new Quotient
@*/
dErr dQuotientCreate(dMesh m,dMeshESH loc,dMeshTag qsizetag,dQuotient *inq)
{
  dQuotient   q;
  MPI_Comm       comm;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(m,dMESH_COOKIE,1);
  dValidPointer(inq,4);
  *inq = 0;
  err = PetscObjectGetComm((PetscObject)m,&comm);dCHK(err);
#if !defined(PETSC_USE_DYNAMIC_LIBRARIES)
  err = dQuotientInitializePackage(PETSC_NULL);dCHK(err);
#endif
  err = PetscHeaderCreate(q,p_dQuotient,struct _dQuotientOps,dQUOTIENT_COOKIE,0,"dQuotient",comm,dQuotientDestroy,dQuotientView);dCHK(err);

  q->mesh          = m;
  q->loc           = loc;
  q->qsizetag      = qsizetag;
  q->setdegreefunc = dQuotientSetDegreeConst;
  q->setdegreectx  = 0;
  q->setdegreeset  = PETSC_FALSE;

  *inq = q;
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientSetType"
/*@
   dQuotientSetType -

@*/
dErr dQuotientSetType(dQuotient q,const dQuotientType type)
{
  dErr (*r)(dQuotient);
  dBool     match;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(q,dQUOTIENT_COOKIE,1);
  PetscValidCharPointer(type,2);
  err = PetscTypeCompare((PetscObject)q,type,&match);dCHK(err);
  if (match) dFunctionReturn(0);
  err =  PetscFListFind(dQuotientList,((PetscObject)q)->comm,type,(void (**)(void)) &r);dCHK(err);
  if (!r) dERROR(PETSC_ERR_ARG_UNKNOWN_TYPE,"Unable to find requested dQuotient type %s",type);
  if (q->ops->destroy) { err = (*q->ops->destroy)(q);dCHK(err); }
  err = PetscMemcpy(q->ops,&dQuotientDefaultOps,sizeof(struct _dQuotientOps));dCHK(err);
  q->setupcalled = 0;
  err = (*r)(q);dCHK(err);
  err = PetscObjectChangeTypeName((PetscObject)q,type);dCHK(err);
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientSetFromOptions"
/**
* Use the options database to determine type and other settings for the quotient
*
* @param q the Quotient
*
* @return err
*/
dErr dQuotientSetFromOptions(dQuotient q)
{
  dBool typeSet,cdegSet;
  dInt cdeg;
  const dQuotientType deft = dQuotientGauss;
  char type[dNAME_LEN];
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(q,dQUOTIENT_COOKIE,1);
  err = PetscOptionsBegin(((PetscObject)q)->comm,((PetscObject)q)->prefix,"Quotient (quadrature rule and element map) options","dQuotient");dCHK(err);
  err = dQuotientRegisterAll(PETSC_NULL);dCHK(err);
  if (((PetscObject)q)->type_name) { deft = ((PetscObject)q)->type_name; }
  err = PetscOptionsList("-dquot_type","Quotient type","dQuotientSetType",dQuotientList,deft,type,dNAME_LEN,&typeSet);dCHK(err);
  err = PetscOptionsInt("-dquot_cdeg","Constant degree","dQuotientSetDegreeConst",5,&cdeg,&cdegSet);dCHK(err);
  if (typeSet) {
    err = dQuotientSetType(q,type);dCHK(err);
  } else if (!((PetscObject)q)->type_name) {
    err = dQuotientSetType(q,deft);dCHK(err);
  }
  if (cdegSet || !q->setdegreeset) {
    q->setdegreefunc = dQuotientSetDegreeConst;
    err = PetscMalloc(1*sizeof(void*),&q->setdegreectx);dCHK(err);
    ((dInt*)q->setdegreectx)[0] = cdeg ? cdeg : 5;
  }
  if (q->ops->setfromoptions) {
    err = (*q->ops->setfromoptions)(q);dCHK(err);
  }
  err = PetscOptionsEnd();
  dFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "dQuotientSetUp"
/*@
   dQuotientSetUp -

@*/
dErr dQuotientSetUp(dQuotient q)
{
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(q,dQUOTIENT_COOKIE,1);
  err = dQuotientSetUp_Private(q);dCHK(err);
  if (!q->setupcalled) {
    err = (*q->ops->setup)(q);dCHK(err);
  }
  if (q->ops->update) {
    err = (*q->ops->update)(q);
  }
  q->setupcalled = 1;
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientSetUp_Private"
/**
* This function should determine which subdomain of the mesh is in use.  For now, we assume that it is always the whole
* mesh.
*
* @param q The quotient
*
* @return err
*/
dErr dQuotientSetUp_Private(dQuotient q)
{
  dErr err;

  dFunctionBegin;
  err = PetscPrintf(((dObject)q)->comm,"%s()\n",__func__);dCHK(err);
  dFunctionReturn(0);
}


#define dQuotientRegisterDynamic(a,b,c,d) dQuotientRegister(a,b,c,d)

#undef __FUNCT__
#define __FUNCT__ "dQuotientRegister"
/**
* Register a new Quotient type in the quotient list.  It can then be selected with \fn dQuotientSetType or on the
* command line.
*
* @param name type name
* @param path source file containing implementation (\p cname)
* @param cname name of function to create an object of this type
* @param function pointer to function \p cname
*
* @return err
*/
dErr dQuotientRegister(const char name[],const char path[],const char cname[],dErr (*function)(dQuotient))
{
  char           fullname[PETSC_MAX_PATH_LEN];
  dErr err;

  dFunctionBegin;
  err = PetscFListConcat(path,cname,fullname);dCHK(err);
  err = PetscFListAdd(&dQuotientList,name,fullname,(void (*)(void))function);dCHK(err);
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientRegisterAll"
/*@
   dQuotientRegisterAll -

@*/
dErr dQuotientRegisterAll(const char path[])
{
  static dBool registered = PETSC_FALSE;
  dErr err;

  dFunctionBegin;
  if (registered) dFunctionReturn(0);
  err = dQuotientRegisterDynamic(dQuotientGauss,path,"dQuotientCreate_Gauss",dQuotientCreate_Gauss);dCHK(err);
  registered = PETSC_TRUE;
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientDestroy"
/*@
   dQuotientDestroy -

@*/
dErr dQuotientDestroy(dQuotient q)
{
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(q,dQUOTIENT_COOKIE,1);
  if (q->ops->destroy) {
    err = (*q->ops->destroy)(q);dCHK(err);
  }
  err = PetscFree(q->quad);dCHK(err);
  err = PetscFree(q->emap);dCHK(err);
  err = PetscHeaderDestroy(q);dCHK(err);
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientGetType"
/*@
   dQuotientGetType -

@*/
dErr dQuotientGetType(dQuotient q,const dQuotientType *type)
{

  dFunctionBegin;
  PetscValidHeaderSpecific(q,dQUOTIENT_COOKIE,1);
  dValidPointer(type,2);
  *type = ((PetscObject)q)->type_name;
  dFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "dQuotientView"
/*@
   dQuotientView -

@*/
dErr dQuotientView(dQuotient q,PetscViewer viewer)
{
  const dQuotientType  type;
  char                   *tagname;
  dBool              iascii;
  dErr          err;

  dFunctionBegin;
  PetscValidHeaderSpecific(q,dQUOTIENT_COOKIE,1);
  if (!viewer) {
    err = PetscViewerASCIIGetStdout(((PetscObject)q)->comm,&viewer);dCHK(err);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_COOKIE,2);
  PetscCheckSameComm(q,1,viewer,2);
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);dCHK(err);
  if (iascii) {
    err = dQuotientGetType(q,&type);dCHK(err);
    if (((PetscObject)q)->prefix) {
      err = PetscViewerASCIIPrintf(viewer,"dQuotient object:(%s)\n",((PetscObject)q)->prefix);dCHK(err);
    } else {
      err = PetscViewerASCIIPrintf(viewer,"dQuotient object:\n");dCHK(err);
    }
    err = PetscViewerASCIIPrintf(viewer,"type: %s\n",(type ? type : "not yet set"));dCHK(err);
    if (q->ops->view) {
      err = PetscViewerASCIIPushTab(viewer);dCHK(err);
      err = (*q->ops->view)(q,viewer);dCHK(err);
      err = PetscViewerASCIIPopTab(viewer);dCHK(err);
    }
    err = PetscViewerASCIIPrintf(viewer,"numbor of elements: %d\n",q->nelems);dCHK(err);
    if (q->qsizetag) {
      err = dMeshGetTagName(q->mesh,q->qsizetag,&tagname);dCHK(err);
      err = PetscViewerASCIIPrintf(viewer,"quadrature size tag: %s\n",tagname);dCHK(err);
    }
  } else {
    if (q->ops->view) {
      err = (*q->ops->view)(q,viewer);dCHK(err);
    }
  }
  err = dMeshView(q->mesh,viewer);dCHK(err);
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientInitializePackage"
/*@
   dQuotientInitializePackage -

@*/
dErr dQuotientInitializePackage(const char path[])
{
  static dBool initialized = PETSC_FALSE;
  dErr err;

  dFunctionBegin;
  if (initialized) dFunctionReturn(0);
  err = PetscCookieRegister("Quotient map",&dQUOTIENT_COOKIE);dCHK(err);
  err = dQuotientRegisterAll(path);dCHK(err);
  initialized = PETSC_TRUE;
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientSetSetDegree"
/*@
   dQuotientSetSetDegree -

@*/
dErr dQuotientSetSetDegree(dQuotient q,dQuotientSetDegreeFunc func,void *ctx)
{
  //dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(q,dQUOTIENT_COOKIE,1);
  dValidPointer(func,2);
  dValidPointer(ctx,3);
  q->setdegreefunc = func;
  q->setdegreectx  = ctx;
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dQuotientSetDegreeConst"
/*@
   dQuotientSetDegreeConst -

@*/
dErr dQuotientSetDegreeConst(dUNUSED dQuotient quot,void *vval,dInt n,dInt *degree)
{
  dInt i,*val=(dInt *)vval;

  dFunctionBegin;
  for (i=0; i<n; i++) {
    degree[3*i] = degree[3*i+1] = degree[3*i+2] = val[0];
  }
  dFunctionReturn(0);
}

dErr dQuotientGetMesh(dQuotient quot,dMesh *mesh)
{

  dFunctionBegin;
  dValidHeader(quot,dQUOTIENT_COOKIE,1);
  dValidPointer(mesh,2);
  *mesh = quot->mesh;
  dFunctionReturn(0);
}
