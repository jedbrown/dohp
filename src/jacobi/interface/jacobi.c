#include <dohpjacimpl.h>
#include <dohp.h>
#include <dohpgeom.h>

dClassId dJACOBI_CLASSID,dQUADRATURE_CLASSID;
PetscLogEvent dLOG_RuleComputeGeometry,dLOG_EFSApply;
PetscBool dJacobiRegisterAllCalled;

const char *dGaussFamilies[] = {"gauss","lobatto","radau","dGaussFamily","dGAUSS_",0};

static dBool dJacobiPackageInitialized;
PetscFList dJacobiList = NULL;

static const struct _dJacobiOps _defaultOps = {
  .View = 0,
  .SetUp = 0,
  .SetFromOptions = 0,
  .Destroy = 0
};

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
  err = PetscHeaderCreate(jac,p_dJacobi,struct _dJacobiOps,dJACOBI_CLASSID,0,"dJacobi",comm,dJacobiDestroy,dJacobiView);dCHK(err);

  jac->basisdegree = 10;
  jac->ruleexcess  = 5;
  jac->setupcalled = 0;
  jac->data        = 0;
  err = PetscMemcpy(jac->ops,&_defaultOps,sizeof(struct _dJacobiOps));dCHK(err);

  *injacobi = jac;
  dFunctionReturn(0);
}

dErr dJacobiSetType(dJacobi jac,const dJacobiType type)
{
  dErr err,(*r)(dJacobi);
  dBool      match;

  dFunctionBegin;
  PetscValidHeaderSpecific(jac,dJACOBI_CLASSID,1);
  PetscValidCharPointer(type,2);
  err = PetscTypeCompare((PetscObject)jac,type,&match);dCHK(err);
  if (match) dFunctionReturn(0);
  if (!dJacobiRegisterAllCalled) {err = dJacobiRegisterAll(NULL);dCHK(err);}
  err = PetscFListFind(dJacobiList,((PetscObject)jac)->comm,type,dTRUE,(void(**)(void))&r);dCHK(err);
  if (!r) dERROR(PETSC_COMM_SELF,1,"Unable to find requested dJacobi type %s",type);
  if (jac->ops->Destroy) { err = (*jac->ops->Destroy)(jac);dCHK(err); }
  err = PetscMemcpy(jac->ops,&_defaultOps,sizeof(_defaultOps));dCHK(err);
  jac->setupcalled = 0;
  err = (*r)(jac);dCHK(err);
  err = PetscObjectChangeTypeName((PetscObject)jac,type);dCHK(err);
  dFunctionReturn(0);
}

dErr dJacobiGetType(dJacobi jac,const dJacobiType *type)
{

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_CLASSID,1);
  dValidPointer(type,2);
  *type = ((PetscObject)jac)->type_name;
  dFunctionReturn(0);
}

dErr dJacobiSetFromOptions(dJacobi jac)
{
  char type[dNAME_LEN] = dJACOBI_TENSOR;
  dBool  typeSet;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(jac,dJACOBI_CLASSID,1);
  err = PetscOptionsBegin(((PetscObject)jac)->comm,((PetscObject)jac)->prefix,"Jacobi options (type and size of basis)","dJacobi");dCHK(err);
  err = PetscOptionsList("-djac_type","Basis type","dJacobiSetType",dJacobiList,
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

dErr dJacobiDestroy(dJacobi *jac)
{
  dErr err;

  dFunctionBegin;
  if (!*jac) dFunctionReturn(0);
  PetscValidHeaderSpecific(*jac,dJACOBI_CLASSID,1);
  if (--((PetscObject)*jac)->refct > 0) dFunctionReturn(0);
  if ((*jac)->ops->Destroy) {
    err = (*jac)->ops->Destroy(*jac);dCHK(err);
  }
  err = PetscHeaderDestroy(jac);dCHK(err);
  dFunctionReturn(0);
}

dErr dJacobiView(dJacobi jac,PetscViewer viewer)
{
  dBool  iascii;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(jac,dJACOBI_CLASSID,1);
  if (!viewer) {err = PetscViewerASCIIGetStdout(((PetscObject)jac)->comm,&viewer);dCHK(err);}
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(jac,1,viewer,2);

  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);dCHK(err);
  if (iascii) {
    err = PetscViewerASCIIPrintf(viewer,"dJacobi object:(%s)\n",
                                  ((PetscObject)jac)->prefix ? ((PetscObject)jac)->prefix : "no prefix");dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"type: %s\n",
                                  ((PetscObject)jac)->type_name ? ((PetscObject)jac)->type_name : "type not set");dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"max basis degree: %d\n",jac->basisdegree);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"rule excess: %d\n",jac->ruleexcess);dCHK(err);
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
* Propagate the anisotropic values \a deg from the high dimensional entities down to the lower dimensional ones.
*
* @note This assumes a minimum rule, should that be optional?
* @note This interface sucks because the Jacobi shouldn't know about adjacencies.
*
* @param[in] jac Jacobi object
* @param[in] a mesh adjacencies
* @param[in,out] deg orders of all entities, in the order implied by \a a
*/
dErr dJacobiPropagateDown(dJacobi jac,dMeshAdjacency a,dPolynomialOrder deg[])
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_CLASSID,1);
  dValidPointer(a,2);
  dValidPointer(deg,3);
  err = jac->ops->PropagateDown(jac,a,deg);dCHK(err);
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
  dErr err;

  dFunctionBegin;
  if (dJacobiRegisterAllCalled) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"Function has already been called");
  err = dJacobiRegisterDynamic(dJACOBI_TENSOR,path,"dJacobiCreate_Tensor",dJacobiCreate_Tensor);dCHK(err);
  err = dJacobiRegisterDynamic(dJACOBI_MODAL ,path,"dJacobiCreate_Modal" ,dJacobiCreate_Modal);dCHK(err);
  dJacobiRegisterAllCalled = PETSC_TRUE;
  dFunctionReturn(0);
}

dErr dJacobiInitializePackage(const char path[])
{
  dErr err;

  dFunctionBegin;
  if (dJacobiPackageInitialized) dFunctionReturn(0);
  dJacobiPackageInitialized = dTRUE;
  err = PetscClassIdRegister("Jacobi context",&dJACOBI_CLASSID);dCHK(err);
  err = PetscClassIdRegister("Quadrature context",&dQUADRATURE_CLASSID);dCHK(err);
  err = PetscLogEventRegister("dEFSApply",       dJACOBI_CLASSID,&dLOG_EFSApply);dCHK(err);
  err = PetscLogEventRegister("dRuleComputeGeom",dJACOBI_CLASSID,&dLOG_RuleComputeGeometry);dCHK(err);
  err = dJacobiRegisterAll(path);dCHK(err);
  err = dQuadratureRegisterAll(path);dCHK(err);
  err = PetscRegisterFinalize(dJacobiFinalizePackage);dCHK(err);
  dFunctionReturn(0);
}

dErr dJacobiFinalizePackage(void)
{
  dFunctionBegin;
  dJacobiPackageInitialized = dFALSE;
  dJacobiList = NULL;
  dQuadratureList = NULL;
  dJacobiRegisterAllCalled = dFALSE;
  dQuadratureRegisterAllCalled = dFALSE;
  dFunctionReturn(0);
}

/** Allocate and fill an array of dEFS with prescribed order, compatible with the given dRules.
 *
 * @bug This function needs another parameter to allow addressing faces instead of just regions.  That is, the rules
 * should correspond to arbitrary faces of the elements, instead of always the elements themselves.
 *
 * @param[in] jac Jacobi
 * @param[in] n number of EFS to extract
 * @param[in] topo topology of elements
 * @param[in] order polynomial order of the element
 * @param[in] rules array of quadrature rules
 * @param[out] efs new array of element function spaces, free with dFree()
 */
dErr dJacobiGetEFS(dJacobi jac,dInt n,const dEntTopology topo[],const dPolynomialOrder order[],const dRule rules[],dEFS **efs)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_CLASSID,1);
  dValidPointer(topo,3);
  dValidPointer(order,4);
  dValidPointer(rules,5);
  dValidPointer(efs,6);
  err = dMallocA(n,efs);dCHK(err);
  err = jac->ops->GetEFS(jac,n,topo,order,rules,*efs);dCHK(err);
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
dErr dJacobiGetNodeCount(dJacobi jac,dInt count,const dEntTopology top[],const dPolynomialOrder deg[],dInt inode[],dInt xnode[])
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_CLASSID,1);
  dValidPointer(top,3);
  dValidPointer(deg,4);
  err = jac->ops->GetNodeCount(jac,count,top,deg,inode,xnode);dCHK(err);
  dFunctionReturn(0);
}

/** Get preferred quadrature type for this Jacobi.
**/
dErr dJacobiGetQuadrature(dJacobi jac,dQuadratureMethod method,dQuadrature *quad)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_CLASSID,1);
  dValidPointer(quad,3);
  if (!jac->quad[method]) {
    if (jac->ops->GetQuadrature) {
      err = jac->ops->GetQuadrature(jac,method,&jac->quad[method]);dCHK(err);
    } else {
      err = dQuadratureCreate(((PetscObject)jac)->comm,&jac->quad[method]);dCHK(err);
      err = dQuadratureSetFromOptions(jac->quad[method]);dCHK(err);
      err = dQuadratureSetMethod(jac->quad[method],method);dCHK(err);
    }
  }
  *quad = jac->quad[method];
  dFunctionReturn(0);
}

dErr dRuleView(dRule rule,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  if (!viewer) {err = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer);dCHK(err);}
  dValidHeader(viewer,PETSC_VIEWER_CLASSID,2);
  err = (*rule->ops.view)(rule,viewer);dCHK(err);
  dFunctionReturn(0);
}

/** Get topological dimension of rule and total number of nodes
 */
dErr dRuleGetSize(dRule rule,dInt *dim,dInt *nnodes)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  err = (*rule->ops.getSize)(rule,dim,nnodes);dCHK(err);
  dFunctionReturn(0);
}

/** Decompose a composite rule into patches that preserve additional sparsity.
 *
 * @note Rationale: There is a one-many relationship between elements and patches.  High-order quadratures cannot be
 * decomposed and are used during traversals that require it.  When a high-order rule has embedded low-order rules,
 * additional sparsity can be preserved by using the low-order rules instead of the single high-order rule.  However,
 * such low-order quadrature involves smaller inner loops and may not be the fastest way to evaluate residuals (where
 * tensor products can be used and sparsity does not appear explicitly).
 */
dErr dRuleGetPatches(dRule rule,dInt *npatches,dInt *patchsize,const dInt **ind,const dReal **weights)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  if (rule->ops.getPatches) {
    err = (*rule->ops.getPatches)(rule,npatches,patchsize,ind,weights);dCHK(err);
  } else dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Rule does not have patches.  Maybe trivial patches should be implemented?");
  dFunctionReturn(0);
}

/** Get location of nodes and associated weights on reference element
 */
dErr dRuleGetNodeWeight(dRule rule,dReal *coord,dReal *weight)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  if (rule->ops.getNodeWeight) {
    err = (*rule->ops.getNodeWeight)(rule,coord,weight);dCHK(err);
  } else if (rule->ops.getTensorNodeWeight) {
    dInt dim,nnodes[3];
    const dReal *tcoord[3],*tweight[3];
    err = (*rule->ops.getTensorNodeWeight)(rule,&dim,nnodes,tcoord,tweight);dCHK(err);
    switch (dim) {
      case 3:
        for (dInt i=0; i<nnodes[0]; i++) {
          for (dInt j=0; j<nnodes[1]; j++) {
            for (dInt k=0; k<nnodes[2]; k++) {
              dInt idx = (i*nnodes[1]+j)*nnodes[2]+k;
              if (coord) {
                coord[idx*3+0] = tcoord[0][i];
                coord[idx*3+1] = tcoord[1][j];
                coord[idx*3+2] = tcoord[2][k];
              }
              if (weight) {
                weight[idx] = tweight[0][i]*tweight[1][j]*tweight[2][k];
              }
            }
          }
        }
        break;
      default:
        dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"dimension %d",dim);
    }
  } else dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"deficient rule");
  dFunctionReturn(0);
}

/** Get tensor-product representation of coordinate nodes and weights on reference element
 */
dErr dRuleGetTensorNodeWeight(dRule rule,dInt *dim,dInt nnodes[3],const dReal *coord[3],const dReal *weight[3])
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  err = (*rule->ops.getTensorNodeWeight)(rule,dim,nnodes,coord,weight);dCHK(err);
  dFunctionReturn(0);
}

/** Compute coordinate transforms required to differentiate and integrate in physical space.
* @note This should not be the responsibility of the Rule, it really doesn't depend on internal details, and instead
* depends on how the coordinates and derivatives are evaluated on the quadrature points.
*
* @param[in] rule to get geometries on
* @param[in] vtx coordinates of corners of this element
* @param[out] global coordinates of quadrature points
* @param[out] jinv derivative of reference coordinates with respect to physical coordinates
* @param[out] jdet determinant of Jacobian (derivative of physical w.r.t. reference coordinates)
**/
dErr dRuleComputeGeometry(dRule rule,const dReal vtx[restrict][3],dReal qg[restrict][3],dReal jinv[restrict][3][3],dReal jdet[restrict])
{
  dErr err;

  dFunctionBegin;
  dValidPointer(rule,1);
  err = PetscLogEventBegin(dLOG_RuleComputeGeometry,0,0,0,0);dCHK(err);
  err = rule->ops.computeGeometry(rule,vtx,qg,jinv,jdet);dCHK(err);
  err = PetscLogEventEnd(dLOG_RuleComputeGeometry,0,0,0,0);dCHK(err);
  dFunctionReturn(0);
}

dErr dRuleComputePhysical(dRule rule,const dScalar jac[],dReal jinv[],dReal jw[])
{
  dErr err;
  dInt n;

  dFunctionBegin;
  dValidPointer(rule,1);
  dValidPointer(jac,2);
  err = dRuleGetSize(rule,NULL,&n);dCHK(err);
  err = dRuleGetNodeWeight(rule,NULL,jw);dCHK(err);
  if (jac) {
    for (dInt i=0; i<n; i++) {
      dReal det;
      err = dGeomInvert3(jac+i*9,jinv+i*9,&det);dCHK(err);
      if (det < 0.) dERROR(PETSC_COMM_SELF,1,"Negative Jacobian at node %D");
      jw[i] *= det;
    }
  } else {
    for (dInt i=0; i<n; i++) {
      jinv[9*i+0] = 1.;
      jinv[9*i+1] = 0.;
      jinv[9*i+2] = 0.;
      jinv[9*i+3] = 0.;
      jinv[9*i+4] = 1.;
      jinv[9*i+5] = 0.;
      jinv[9*i+6] = 0.;
      jinv[9*i+7] = 0.;
      jinv[9*i+8] = 1.;
    }
  }
  dFunctionReturn(0);
}

dErr dEFSView(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  if (!viewer) {err = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer);dCHK(err);}
  dValidHeader(viewer,PETSC_VIEWER_CLASSID,2);
  err = (*efs->ops.view)(efs,viewer);dCHK(err);
  dFunctionReturn(0);
}

dErr dEFSGetSizes(dEFS efs,dInt *dim,dInt *inodes,dInt *total)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  err = (*efs->ops.getSizes)(efs,dim,inodes,total);dCHK(err);
  dFunctionReturn(0);
}

/** Get node locations (on reference element) as tensor product
*
* @param[in] dim[out] dimension of tensor product
* @param[in,out] tsize Array of length 3, on return, holds each size in each direction of tensor product
* @param[in,out] nodes Array of length 3, on return, holds pointers to node locations in each direction of tensor product
* @param[in,out] weight Array of length 3, on return, holds pointers to node weights in each direction of tensor product (quadrature on same nodes)
* @param[in,out] mscale Array of length 3 or NULL, on return if not NULL, holds scaling that can be applied within tensor product to improve conditioning of sparse mass matrix
* @param[in,out] lscale Array of length 3 or NULL, on return if not NULL, holds scaling that can be applied within tensor product to improve conditioning of sparse Laplacian matrix
**/
dErr dEFSGetTensorNodes(dEFS efs,dInt *dim,dInt *tsize,dReal *nodes[],dReal *weight[],const dReal *mscale[],const dReal *lscale[])
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  err = (*efs->ops.getTensorNodes)(efs,dim,tsize,nodes,weight,mscale,lscale);dCHK(err);
  dFunctionReturn(0);
}

/** Get the location of nodes for a nodal EFS (only works for tensor product)
*
* @param[in] x Geometry vector, usually obtained from dFSGetElements()
* @param[out] dim Spatial dimension on which the coordinates lie
* @param[out] P Tensor product size for nodes
* @param[out] qx Location of interpolation nodes
*
**/
dErr dEFSGetGlobalCoordinates(dEFS efs,const dReal x[restrict][3],dInt *dim,dInt P[3],dReal (*qx)[3])
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  dValidPointer(x,2);
  dValidPointer(dim,3);
  dValidPointer(P,4);
  dValidPointer(qx,5);
  err = (*efs->ops.getGlobalCoordinates)(efs,x,dim,P,qx);dCHK(err);
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
  if (dofs < 1) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"dofs %D, must be positive",dofs);
  dValidScalarPointer(in,4);
  dValidScalarPointer(out,5);
#if defined dUSE_LOG_FINEGRAIN
  err = PetscLogEventBegin(dLOG_EFSApply,0,0,0,0);dCHK(err);
#endif
  err = (*efs->ops.apply)(efs,mapdata,dofs,in,out,amode,imode);dCHK(err);
#if defined dUSE_LOG_FINEGRAIN
  err = PetscLogEventEnd(dLOG_EFSApply,0,0,0,0);dCHK(err);
#endif
  dFunctionReturn(0);
}

dErr dEFSGetExplicit(dEFS efs,const dReal jinv[],dInt *inQ,dInt *inP,const dReal **inbasis,const dReal **inderiv)
{
  dReal *newbasis,*newderiv;
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  if (jinv) dValidPointer(jinv,2);
  dValidPointer(inQ,3);
  dValidPointer(inP,4);
  dValidPointer(inbasis,5);
  dValidPointer(inderiv,6);
  if (efs->ops.getExplicit) {
    err = (*efs->ops.getExplicit)(efs,jinv,inQ,inP,inbasis,inderiv);dCHK(err);
  } else if (efs->ops.apply) {
    dRule rule;
    dScalar *tmp,*basis,*deriv;
    dInt P,Q;
    err = dEFSGetSizes(efs,NULL,NULL,&P);dCHK(err);
    err = dEFSGetRule(efs,&rule);dCHK(err);
    err = dRuleGetSize(rule,NULL,&Q);dCHK(err);
    err = dMallocA(P,&tmp);dCHK(err);
    err = dMallocA2(P*Q,&basis,P*Q*3,&deriv);dCHK(err);
    for (dInt i=0; i<P; i++) {
      err = dMemzero(tmp,P*sizeof(tmp[0]));dCHK(err);
      tmp[i] = 1.;
      err = dEFSApply(efs,jinv,1,tmp,basis+i*Q,dAPPLY_INTERP,INSERT_VALUES);dCHK(err);
      err = dEFSApply(efs,jinv,1,tmp,deriv+i*Q*3,dAPPLY_GRAD,INSERT_VALUES);dCHK(err);
    }
    err = dFree(tmp);dCHK(err);
    err = dMallocA2(P*Q,&newbasis,P*Q*3,&newderiv);dCHK(err);
    /* Transpose and place in user vectors */
    for (dInt i=0; i<Q; i++) {
      const PetscReal *ji = jinv ? &jinv[9*i] : ((const PetscReal[]){1,0,0,0,1,0,0,0,1});
      for (dInt j=0; j<P; j++) {
        PetscReal ux,uy,uz;
        newbasis[i*P+j] = PetscRealPart(basis[i+j*Q]);
        ux = PetscRealPart(deriv[(i+j*Q)*3+0]);
        uy = PetscRealPart(deriv[(i+j*Q)*3+1]);
        uz = PetscRealPart(deriv[(i+j*Q)*3+2]);
        /* chain rule: (\grad u)^T J^{-1} */
        newderiv[(i*P+j)*3+0] = ux*ji[0] + uy*ji[3] + uz*ji[6];
        newderiv[(i*P+j)*3+1] = ux*ji[1] + uy*ji[4] + uz*ji[7];
        newderiv[(i*P+j)*3+2] = ux*ji[2] + uy*ji[5] + uz*ji[8];
      }
    }
    err = dFree2(basis,deriv);dCHK(err);
    *inQ = Q;
    *inP = P;
    *inbasis = newbasis;
    *inderiv = newderiv;
  } else dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"this EFS is deficient");
  dFunctionReturn(0);
}

dErr dEFSRestoreExplicit(dEFS efs,const dReal jinv[],dInt *Q,dInt *P,const dReal **basis,const dReal **deriv)
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  if (jinv) dValidPointer(jinv,2);
  dValidPointer(Q,3);
  dValidPointer(P,4);
  dValidPointer(basis,5);
  dValidPointer(deriv,6);
  err = dFree2(*basis,*deriv);dCHK(err);
  *Q = 0;
  *P = 0;
  *basis = 0;
  *deriv = 0;
  dFunctionReturn(0);
}

/**
 * @param[in] efs Element function space
 * @param[in] npatches Number of patches on this element
 * @param[in] Q Number of quadrature points per patch
 * @param[in] qidx Element-natural indices of quadrature nodes for each patch (concatenated)
 * @param[in] cjinv Jacobian for the whole element in patch-natural ordering (concatenated)
 * @param[in] eoffset Offset of element in expanded ordering, \a eidx results will be offset by this amount
 * @param[out] P Number of basis functions with support on each patch
 * @param[in,out] eidx Element-natural indices of each basis function with support on each patch (concatenated)
 * @param[in,out] interp Explicit interpolation matrices for each patch (concatenated)
 * @param[in,out] deriv Explicit derivative evaluation matrices for each patch (concatenated)
 */
dErr dEFSGetExplicitSparse(dEFS efs,dInt npatches,dInt Q,const dInt qidx[],const dReal cjinv[],dInt eoffset,dInt *P,dInt eidx[],dReal interp[],dReal deriv[])
{
  dErr err;

  dFunctionBegin;
  dValidPointer(efs,1);
  dASSERT(npatches > 0);
  dASSERT(Q > 0);
  dValidIntPointer(qidx,4);
  if (cjinv) dValidRealPointer(cjinv,5);
  dValidPointer(P,7);
  dValidPointer(eidx,8);
  dValidPointer(interp,9);
  dValidPointer(deriv,10);
  if (!efs->ops.getExplicitSparse) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"Not available");
  err = (*efs->ops.getExplicitSparse)(efs,npatches,Q,qidx,cjinv,eoffset,P,eidx,interp,deriv);dCHK(err);
  dFunctionReturn(0);
}

dErr dEFSGetExplicitSparse_Basic(dEFS efs,dInt npatches,dInt Q,const dInt qidx[],const dReal cjinv[],dInt eoffset,dInt *P,dInt eidx[],dReal interp[],dReal deriv[])
{
  dErr err;
  dReal *ecjinv = NULL;
  const dReal *einterp,*ederiv;
  dInt eQ,eP;

  dFunctionBegin;
  if (cjinv) {
    /* Translate cjinv into element-natural ordering */
    err = dMallocA(npatches*Q*9,&ecjinv);dCHK(err);
    for (dInt i=0; i<npatches*Q; i++) {
      for (dInt j=0; j<9; j++) ecjinv[i*9+j] = cjinv[qidx[i]*9+j];
    }
  }
  err = dEFSGetExplicit(efs,ecjinv,&eQ,&eP,&einterp,&ederiv);dCHK(err);
  *P = eP;
  for (dInt i=0; i<npatches; i++) {
    for (dInt k=0; k<eP; k++) eidx[i*eP+k] = eoffset + k;
    for (dInt j=0; j<Q; j++) {
      const dInt
        pq = i*Q+j,             /* quadrature node in patch-natural ordering */
        eq = qidx[pq];          /* quadrature node in element-natural ordering */
      for (dInt k=0; k<eP; k++) {
        interp[pq*eP+k] = einterp[eq*eP+k];
        for (dInt l=0; l<3; l++) deriv[(pq*eP+k)*3+l] = ederiv[(eq*eP+k)*3+l];
      }
    }
  }
  err = dEFSRestoreExplicit(efs,ecjinv,&eQ,&eP,&einterp,&ederiv);dCHK(err);
  err = dFree(ecjinv);dCHK(err);
  dFunctionReturn(0);
}


/** Get constraint counts for element assembly matrices.
*
* @note This function (and dJacobiAddConstraints()) will be difficult to implement on mixed spaces.  Looking for a
* better interface.
*
* @param jac Jacobi object
* @param nx Number of entities in expanded vector
* @param xi Index of expanded entity in full adjacency (length \c nx)
* @param xs Offset in expanded vector of first node associated with each expanded entity (length \c nx+1)
* @param is Starting index of interior nodes for each entity in \a ma
* @param deg Basis degree for each entity in \a ma
* @param ma MeshAdjacency object
* @param nnz Number of nonzeros per row of element assembly matrix
* @param pnnz Number of nonzeros per row of preconditioning element assembly matrix
*/
dErr dJacobiGetConstraintCount(dJacobi jac,dInt nx,const dInt xi[],const dInt xs[],const dInt is[],const dPolynomialOrder deg[],dMeshAdjacency ma,dInt nnz[],dInt pnnz[])
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_CLASSID,1);
  dValidPointer(xi,3);
  dValidPointer(xs,4);
  dValidPointer(is,5);
  dValidPointer(deg,6);
  dValidPointer(ma,7);
  dValidPointer(nnz,8);
  dValidPointer(pnnz,9);
  err = (*jac->ops->GetConstraintCount)(jac,nx,xi,xs,is,deg,ma,nnz,pnnz);dCHK(err);
  dFunctionReturn(0);
}

/** Actually assemble the element assembly matrices, see dJacobiGetConstraintCount()
*
* @note This interface sucks for mixed spaces.  Looking for something better.
**/
dErr dJacobiAddConstraints(dJacobi jac,dInt nx,const dInt xi[],const dInt xs[],const dInt is[],const dPolynomialOrder deg[],dMeshAdjacency ma,Mat E,Mat Ep)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(jac,dJACOBI_CLASSID,1);
  dValidPointer(xi,3);
  dValidPointer(xs,4);
  dValidPointer(is,5);
  dValidPointer(deg,6);
  dValidPointer(ma,7);
  dValidHeader(E,MAT_CLASSID,8);
  dValidHeader(Ep,MAT_CLASSID,9);
  err = (*jac->ops->AddConstraints)(jac,nx,xi,xs,is,deg,ma,E,Ep);dCHK(err);
  dFunctionReturn(0);
}
