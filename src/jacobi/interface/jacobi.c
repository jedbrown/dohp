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
#include "private/jacimpl.h"

PetscCookie dJACOBI_COOKIE;
static PetscFList dJacobiList = 0;

static const struct _dJacobiOps _defaultOps = {
  .View = 0,
  .SetUp = 0,
  .SetFromOptions = 0,
  .Destroy = 0
};

/**
* Create a new Jacobi object and initialize with defaults.
*
* @param comm
* @param injacobi
*
* @return
*/
dErr dJacobiCreate(MPI_Comm comm,dJacobi *injacobi)
{
  dJacobi jac;
  dErr err;

  dFunctionBegin;
  dValidPointer(injacobi,2);
  *injacobi = 0;
#if !defined(PETSC_USE_DYNAMIC_LIBRARIES)
  err = dJacobiInitializePackage(PETSC_NULL);dCHK(err);
#endif
  err = PetscHeaderCreate(jac,p_dJacobi,struct _dJacobiOps,dJACOBI_COOKIE,0,"dJacobi",comm,dJacobiDestroy,dJacobiView);dCHK(err);

  jac->basisdegree = 10;
  jac->ruleexcess = 5;
  jac->setupcalled = 0;
  jac->impl = 0;
  err = PetscMemcpy(jac->ops,&_defaultOps,sizeof(struct _dJacobiOps));dCHK(err);

  *injacobi = jac;
  dFunctionReturn(0);
}

/**
* Set the type for a dJacobi object.
*
* @param jac
* @param type
*
* @return
*/
dErr dJacobiSetType(dJacobi jac,dJacobiType type)
{
  dErr err,(*r)(dJacobi);
  dBool     match;

  dFunctionBegin;
  PetscValidHeaderSpecific(jac,dJACOBI_COOKIE,1);
  PetscValidCharPointer(type,2);
  err = PetscTypeCompare((PetscObject)jac,type,&match);dCHK(err);
  if (match) dFunctionReturn(0);
  err = PetscFListFind(dJacobiList,((PetscObject)jac)->comm,type,(void(**)(void))&r);dCHK(err);
  if (!r) dERROR(1,"Unable to find requested dJacobi type %s",type);
  if (jac->ops->Destroy) { err = (*jac->ops->Destroy)(jac);dCHK(err); }
  err = PetscMemcpy(jac->ops,&_defaultOps,sizeof(struct _dJacobiOps));dCHK(err);
  jac->setupcalled = 0;
  err = (*r)(jac);dCHK(err);
  err = PetscObjectChangeTypeName((PetscObject)jac,type);dCHK(err);
  dFunctionReturn(0);
}

/**
* Set options from the options database.
*
* @param jac
*
* @return
*/
dErr dJacobiSetFromOptions(dJacobi jac)
{
  char type[dNAME_LEN] = dJACOBI_TENSOR;
  dBool typeSet;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(jac,dJACOBI_COOKIE,1);
  err = PetscOptionsBegin(((PetscObject)jac)->comm,((PetscObject)jac)->prefix,"Jacobi options (type and size of basis/quadrature rules)","dJacobi");dCHK(err);
  err = PetscOptionsList("-djac_type","Basis/Quadrature type","dJacobiSetType",dJacobiList,
                          (((PetscObject)jac)->type_name?((PetscObject)jac)->type_name:type),type,dNAME_LEN,&typeSet);dCHK(err);
  if (typeSet) {
    err = dJacobiSetType(jac,type);dCHK(err);
  }
  if (!((PetscObject)jac)->type_name) {
    err = dJacobiSetType(jac,type);dCHK(err);
  }
  err = PetscOptionsInt("-djac_basis_degree","Max basis degree","dJacobiSetDegrees",jac->basisdegree,&jac->basisdegree,PETSC_NULL);dCHK(err);
  err = PetscOptionsInt("-djac_rule_excess","Excess quadrature points","dJacobiSetDegrees",jac->ruleexcess,&jac->ruleexcess,PETSC_NULL);dCHK(err);
  if (jac->ops->SetFromOptions) {
    err = jac->ops->SetFromOptions(jac);dCHK(err);
  }
  err = PetscOptionsEnd();dCHK(err);
  dFunctionReturn(0);
}

/**
* Initialize the Jacobi object.
*
* @param jac
*
* @return
*/
dErr dJacobiSetUp(dJacobi jac)
{
  dBool flg;
  dErr  err;

  dFunctionBegin;
  PetscValidHeaderSpecific(jac,dJACOBI_COOKIE,1);
  if ((!jac->setupcalled) && jac->ops->SetUp) {
    err = jac->ops->SetUp(jac);dCHK(err);
  }
  jac->setupcalled = 1;

  /* View if requested */
  err = PetscOptionsHasName(((dObject)jac)->prefix,"-djac_view",&flg);dCHK(err);
  if (flg) {
    dViewer viewer;
    err = PetscViewerASCIIGetStdout(((dObject)jac)->comm,&viewer);dCHK(err);
    err = dJacobiView(jac,viewer);dCHK(err);
  }
  dFunctionReturn(0);
}

/**
* Destroy a Jacobi object.
*
* @param jac
*
* @return
*/
dErr dJacobiDestroy(dJacobi jac)
{
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(jac,dJACOBI_COOKIE,1);
  if (jac->ops->Destroy) {
    err = jac->ops->Destroy(jac);dCHK(err);
  }
  err = PetscHeaderDestroy(jac);dCHK(err);
  dFunctionReturn(0);
}

/**
* View the state of a dJacobi.
*
* @param jac
* @param viewer
*
* @return
*/
dErr dJacobiView(dJacobi jac,PetscViewer viewer)
{
  dBool iascii;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(jac,dJACOBI_COOKIE,1);
  if (!viewer) {
    err = PetscViewerASCIIGetStdout(((PetscObject)jac)->comm,&viewer);dCHK(err);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_COOKIE,2);
  PetscCheckSameComm(jac,1,viewer,2);

  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);dCHK(err);
  if (iascii) {
    err = PetscViewerASCIIPrintf(viewer,"dJacobi object:(%s)\n",
                                  ((PetscObject)jac)->prefix ? ((PetscObject)jac)->prefix : "no prefix");dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"type: %s\n",
                                  ((PetscObject)jac)->type_name ? ((PetscObject)jac)->type_name : "type not set");dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"max basis degree: %d\n",jac->basisdegree);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"rule excess: %d\n",jac->ruleexcess);dCHK(err);
    if (!jac->setupcalled) {
      err = PetscViewerASCIIPrintf(viewer,"Object has not been set up.\n",jac->basisdegree);dCHK(err);
    }
    if (jac->ops->View) {
      err = (*jac->ops->View)(jac,viewer);dCHK(err);
    } else {
      err = PetscViewerASCIIPrintf(viewer,"Internal info not available.\n");dCHK(err);
    }
    err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  } else if (jac->ops->View) {
    err = (*jac->ops->View)(jac,viewer);dCHK(err);
  }
  dFunctionReturn(0);
}

/**
* Propogate an the anisotropic values \a v from this entity to the facets \a adj which have current values \a av.
*
* This assumes a minimum rule, should that be optional?
*
* @param jac element function space
* @param topo topology of source entity
* @param econn connectivity (vertex handles) of source entity
* @param edeg values on source entity (array of length 3)
* @param conn connectivity of adjacent entities (packed, in order)
* @param ind index of entity in \a deg
* @param deg flat array of values on adjacent entities, indexed by \a ind
*
* @return error (if connectivity or adjacency is wrong)
*/
dErr dJacobiPropogateDown(dJacobi jac,const struct dMeshAdjacency *a,dInt deg[])
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_COOKIE,1);
  dValidPointer(a,2);
  dValidPointer(deg,3);
  err = jac->ops->PropogateDown(jac,a,deg);dCHK(err);
  dFunctionReturn(0);
}

dErr dJacobiRegister(const char name[],const char path[],const char cname[],dErr(*create)(dJacobi))
{
  char fullname[dMAX_PATH_LEN];
  dErr err;

  dFunctionBegin;
  err = PetscFListConcat(path,cname,fullname);dCHK(err);
  err = PetscFListAdd(&dJacobiList,name,fullname,(void (*)(void))create);dCHK(err);
  dFunctionReturn(0);
}


dErr dJacobiRegisterAll(const char path[])
{
  static dBool called = PETSC_FALSE;
  dErr err;

  dFunctionBegin;
  err = dJacobiRegisterDynamic(dJACOBI_TENSOR,path,"dJacobiCreate_Tensor",dJacobiCreate_Tensor);dCHK(err);
  called = PETSC_TRUE;
  dFunctionReturn(0);
}

dErr dJacobiInitializePackage(const char path[])
{
  static dBool initialized = PETSC_FALSE;
  dErr err;

  dFunctionBegin;
  if (initialized) dFunctionReturn(0);
  err = PetscCookieRegister("Jacobi context",&dJACOBI_COOKIE);dCHK(err);
  err = dJacobiRegisterAll(path);dCHK(err);
  initialized = PETSC_TRUE;
  dFunctionReturn(0);
}


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
dErr dJacobiSetDegrees(dJacobi jac,dInt basisdegree,dInt ruleexcess)
{
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(jac,dJACOBI_COOKIE,1);
  if (jac->setupcalled) {
    if (jac->ops->Destroy) { err = (*jac->ops->Destroy)(jac);dCHK(err); }
    err = PetscMemcpy(jac->ops,&_defaultOps,sizeof(struct _dJacobiOps));dCHK(err);
  }
  jac->basisdegree = basisdegree;
  jac->ruleexcess = ruleexcess;
  dFunctionReturn(0);
}

/** Fill an array of dRule starting at \a firstrule.
*
* @param jac the context
* @param n number of elements
* @param topo topology of the element
* @param rsize number of points in each Cartesian direction
* @param firstrule place to put the newly constructed dRule, normally this will be \c s_dRule[]
*/
dErr dJacobiGetRule(dJacobi jac,dInt n,const dEntTopology topo[],const dInt rsize[],dRule firstrule)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_COOKIE,1);
  dValidPointer(topo,3);
  dValidPointer(rsize,4);
  dValidPointer(firstrule,5);
  err = jac->ops->GetRule(jac,n,topo,rsize,firstrule);dCHK(err);
  dFunctionReturn(0);
}


/** Fill an array of EFS starting at \a firstefs compatible with the array of dRule starting at \a firstrule
*
* @param jac Jacobi
* @param n number of EFS to extract
* @param topo topology of elements
* @param bsize anisotropic basis degree
* @param firstrule handle of first Rule, usually an array
* @param firstefs handle of first EFS, usually an array
*/
dErr dJacobiGetEFS(dJacobi jac,dInt n,const dEntTopology topo[],const dInt bsize[],dRule firstrule,dEFS firstefs)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_COOKIE,1);
  dValidPointer(topo,3)
  dValidPointer(bsize,4);
  dValidPointer(firstrule,5);
  dValidPointer(firstefs,6);
  err = jac->ops->GetEFS(jac,n,topo,bsize,firstrule,firstefs);dCHK(err);
  dFunctionReturn(0);
}

/**
* Get the number of interior and expanded nodes for an array of entities with given topology and degree
*
* @param jac context, defines the function space we are constructing
* @param count number of entities
* @param top topology of entities
* @param deg degree, 3 values each, interlaced
* @param[out] inode number of interior nodes
* @param[out] xnode number of expanded nodes
*
* @return err
*/
dErr dJacobiGetNodeCount(dJacobi jac,dInt count,const dEntTopology top[],const dInt deg[],dInt inode[],dInt xnode[])
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_COOKIE,1);
  dValidPointer(top,3);
  dValidPointer(deg,4);
  err = jac->ops->GetNodeCount(jac,count,top,deg,inode,xnode);dCHK(err);
  dFunctionReturn(0);
}

dErr dRuleView(dRule rule,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  dValidHeader(viewer,PETSC_VIEWER_COOKIE,2);
  err = (*rule->ops->view)(rule,viewer);dCHK(err);
  dFunctionReturn(0);
}

dErr dRuleGetSize(dRule rule,dInt *dim,dInt *nnodes)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  err = (*rule->ops->getSize)(rule,dim,nnodes);dCHK(err);
  dFunctionReturn(0);
}

dErr dRuleGetNodeWeight(dRule rule,dReal *coord,dReal *weight)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  err = (*rule->ops->getNodeWeight)(rule,coord,weight);dCHK(err);
  dFunctionReturn(0);
}

dErr dRuleGetTensorNodeWeight(dRule rule,dInt *dim,dInt *nnodes,const dReal **coord,const dReal **weight)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  err = (*rule->ops->getTensorNodeWeight)(rule,dim,nnodes,coord,weight);dCHK(err);
  dFunctionReturn(0);
}


dErr dEFSView(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  dValidHeader(viewer,PETSC_VIEWER_COOKIE,2);
  err = (*efs->ops->view)(efs,viewer);dCHK(err);
  dFunctionReturn(0);
}

dErr dEFSGetSizes(dEFS efs,dInt *dim,dInt *inodes,dInt *total)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  err = (*efs->ops->getSizes)(efs,dim,inodes,total);dCHK(err);
  dFunctionReturn(0);
}

dErr dEFSGetTensorNodes(dEFS efs,dInt *dim,dInt *tsize,dReal *nodes[])
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  err = (*efs->ops->getTensorNodes)(efs,dim,tsize,nodes);
  dFunctionReturn(0);
}

dErr dEFSGetRule(dEFS efs,dRule *rule)
{

  dFunctionBegin;
  dValidPointer(efs,1);
  dValidPointer(rule,2);
  *rule = efs->rule;
  dFunctionReturn(0);
}

dErr dEFSApply(dEFS efs,const dReal mapdata[],dInt dofs,const dScalar in[],dScalar out[restrict],dApplyMode amode,InsertMode imode)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  err = (*efs->ops->apply)(efs,mapdata,dofs,in,out,amode,imode);dCHK(err);
  dFunctionReturn(0);
}

dErr dJacobiGetConstraintCount(dJacobi jac,dInt nx,const dInt xi[],const dInt xs[],const dInt deg[],const struct dMeshAdjacency *ma,dInt nnz[],dInt pnnz[])
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_COOKIE,1);
  dValidPointer(xi,3);
  dValidPointer(xs,4);
  dValidPointer(deg,5);
  dValidPointer(ma,6);
  dValidPointer(nnz,7);
  dValidPointer(pnnz,8);
  err = (*jac->ops->GetConstraintCount)(jac,nx,xi,xs,deg,ma,nnz,pnnz);dCHK(err);
  dFunctionReturn(0);
}

dErr dJacobiAddConstraints(dJacobi jac,dInt nx,const dInt xi[],const dInt xs[],const dInt is[],const dInt deg[],const struct dMeshAdjacency *ma,Mat C,Mat Cp)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_COOKIE,1);
  dValidPointer(xi,3);
  dValidPointer(xs,4);
  dValidPointer(deg,5);
  dValidPointer(ma,6);
  dValidHeader(C,MAT_COOKIE,7);
  dValidHeader(Cp,MAT_COOKIE,8);
  err = (*jac->ops->AddConstraints)(jac,nx,xi,xs,is,deg,ma,C,Cp);dCHK(err);
  dFunctionReturn(0);
}
