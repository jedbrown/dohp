/**
* @file   tensor.c
* @author Jed Brown <jed@59A2.org>
* @date   Fri Aug 22 20:39:12 2008
*
* @brief  A nodal Tensor product basis.
*
* See the Tensor section for details on the storage method used.
*
*/

#include "tensor.h"
#include "inlinepoly.h"
#include "dohpgeom.h"
#include "dohpmesh.h"           /* for iMesh_TopologyName */

static dErr TensorBuilderCreate(void*,TensorBuilder*);
static dErr TensorBuilderDestroy(TensorBuilder);

static dErr TensorRuleCreate(TensorBuilder,dInt,TensorRule*);
static dErr TensorRuleDestroy(TensorRule);

static dErr TensorBasisCreate(TensorBuilder,const TensorRule,dInt,TensorBasis*);
static dErr TensorBasisDestroy(TensorBasis);

static dErr dJacobiSetUp_Tensor(dJacobi);
static dErr dJacobiDestroy_Tensor(dJacobi);
static dErr dJacobiView_Tensor(dJacobi,PetscViewer);
static dErr dJacobiPropogateDown_Tensor(dJacobi jac,dEntTopology topo,const dMeshEH econn[],const dInt edeg[],const dMeshEH conn[],const dInt ind[],dInt deg[]);
static dErr dJacobiGetRule_Tensor(dJacobi jac,dEntTopology top,const dInt rsize[],dRule *rule,void **base,dInt *index);
static dErr dJacobiGetEFS_Tensor(dJacobi jac,dEntTopology top,const dInt bsize[],dRule rule,dEFS *efs,void **base,dInt *index);

static dErr TensorGetRule(Tensor this,dInt n,TensorRule *out);
static dErr TensorGetBasis(Tensor this,dInt m,dInt n,TensorBasis *out);

static dErr TensorJacobiHasBasis(dJacobi,dInt,dInt,dBool*);

static dErr dRealTableView(dInt m,dInt n,const dReal mat[],const char name[],dViewer viewer);


/**
* Initializes the ops table.
*
* @param jac
*
* @return
*/
dErr dJacobiCreate_Tensor(dJacobi jac)
{
  static const struct _dJacobiOps myops = {
    .setup = dJacobiSetUp_Tensor,
    .setfromoptions = 0,
    .destroy = dJacobiDestroy_Tensor,
    .view = dJacobiView_Tensor,
    .propogatedown = dJacobiPropogateDown_Tensor,
    .getrule = dJacobiGetRule_Tensor,
    .getefs = dJacobiGetEFS_Tensor
  };
  TensorRuleOptions ropt;
  TensorBasisOptions bopt;
  dErr err;

  dFunctionBegin;
  err = dPrintf(((PetscObject)jac)->comm,"dJacobiCreate_Tensor()\n");dCHK(err); /* diagnostic */
  err = dMemcpy(jac->ops,&myops,sizeof(struct _dJacobiOps));dCHK(err);
  err = dNew(struct s_Tensor,&jac->impl);dCHK(err);
  err = dNew(struct s_TensorRuleOptions,&ropt);dCHK(err);
  err = dNew(struct s_TensorBasisOptions,&bopt);dCHK(err);

  ropt->alpha  = 0.0;
  ropt->beta   = 0.0;
  ropt->family = GAUSS;
  ((Tensor)jac->impl)->ruleOpts = ropt;

  bopt->alpha  = 0.0;
  bopt->beta   = 0.0;
  bopt->family = GAUSS_LOBATTO;
  ((Tensor)jac->impl)->basisOpts = bopt;
  dFunctionReturn(0);
}

/**
* Prepare the dJacobi context to return dRule and dEFS objects (with dJacobiGetRule and dJacobiGetEFS).
*
* @param jac the context
*
* @return err
*/
static dErr dJacobiSetUp_Tensor(dJacobi jac)
{
  Tensor this = (Tensor)(jac->impl);
  dInt   M,N;
  dBool  has;
  dErr   err;

  dFunctionBegin;
  if (this->setupcalled) dFunctionReturn(0);
  err = TensorBuilderCreate((void*)this->ruleOpts,&this->ruleBuilder);dCHK(err);
  err = TensorBuilderCreate((void*)this->basisOpts,&this->basisBuilder);dCHK(err);
  this->N = N = jac->basisdegree;                     /* all valid basis degrees are < P */
  this->M = M = N + jac->ruleexcess;                  /* all valid rule degrees are < M */

  err = dMallocM(M,TensorRule,&this->rule);dCHK(err); /* Get space to store all the rule pointers */
  this->rule[0] = NULL;
  for (dInt i=1; i<M; i++) {
    err = TensorRuleCreate(this->ruleBuilder,i,&this->rule[i]);dCHK(err);
  }

  err = dMallocM(M*N,TensorBasis,&this->basis);dCHK(err);
  for (dInt i=0; i<M; i++) {
    TensorRule rule = this->rule[i];
    for (dInt j=0; j<N; j++) {
      err = TensorJacobiHasBasis(jac,i,j,&has);dCHK(err);
      if (!has || !rule) {
        this->basis[i*N+j] = NULL;
        continue;
      }
      err = TensorBasisCreate(this->basisBuilder,rule,j,&this->basis[i*N+j]);dCHK(err);
    }
  }
  err = dJacobiRuleOpsSetUp_Tensor(jac);dCHK(err);
  err = dJacobiEFSOpsSetUp_Tensor(jac);dCHK(err);
  this->setupcalled = true;
  dFunctionReturn(0);
}

static dErr dJacobiDestroy_Tensor(dJacobi jac)
{
  Tensor this = (Tensor)jac->impl;
  dErr err;

  dFunctionBegin;
  if (this->rule) {
    for (dInt i=0; i<this->M; i++) {
      if (this->rule[i]) { err = TensorRuleDestroy(this->rule[i]);dCHK(err); this->rule[i] = NULL; }
    }
    err = dFree(this->rule);dCHK(err);
  }
  if (this->basis) {
    for (dInt i=0; i<this->M; i++) {
      for (dInt j=0; j<this->N; j++) {
        dInt idx = i*this->N + j;
        if (this->basis[idx]) { err = TensorBasisDestroy(this->basis[idx]);dCHK(err); this->basis[idx] = NULL; }
      }
    }
    err = dFree(this->basis);dCHK(err);
  }
  if (this->ruleBuilder) { err = TensorBuilderDestroy(this->ruleBuilder);dCHK(err); }
  if (this->basisBuilder) { err = TensorBuilderDestroy(this->basisBuilder);dCHK(err); }
  if (this->ruleOpts) { err = dFree(this->ruleOpts);dCHK(err); }
  if (this->basisOpts) { err = dFree(this->basisOpts);dCHK(err); }
  err = dJacobiRuleOpsDestroy_Tensor(jac);dCHK(err);
  err = dJacobiEFSOpsDestroy_Tensor(jac);dCHK(err);
  err = dFree(this);dCHK(err);
  dFunctionReturn(0);
}

static dErr dJacobiView_Tensor(dJacobi jac,dViewer viewer)
{
  Tensor this = (Tensor)jac->impl;
  dBool ascii;
  TensorRule r;
  TensorBasis b;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  err = PetscViewerASCIIPrintf(viewer,"Tensor based Jacobi\n");dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"TensorRule database:\n");dCHK(err);
  for (dInt i=0; i<this->M; i++) {
    r = this->rule[i];
    if (r) {err = TensorRuleView(r,viewer);dCHK(err);}
  }
  err = PetscViewerASCIIPrintf(viewer,"TensorBasis database.\n");dCHK(err);
  for (dInt i=1; i<this->M; i++) {
    for (dInt j=1; j<this->N; j++) {
      b = 0;
      err = TensorGetBasis(this,i,j,&b);dCHK(err);
      if (b) {
        err = TensorBasisView(b,viewer);dCHK(err);
      }
    }
  }
  /* view the basis functions next */
  dFunctionReturn(0);
}

static dErr dJacobiPropogateDown_Tensor(dUNUSED dJacobi jac,dEntTopology topo,const dMeshEH econn[],const dInt edeg[],const dMeshEH conn[],const dInt ind[],dInt deg[])
{
  static const dInt quadperm[4] = {0,1,0,1};
  static const dInt hexperm[6][2] = {{0,2},{1,2},{0,2},{1,2},{1,0},{0,1}}; /* map natural axis of Quad to natural axis of Hex */
  static const dInt orient[8][2] = {{0,1},{1,0},{0,1},{1,0},            /* map actual Quad orientation to proper orientation */
                                    {1,0},{0,1},{1,0},{0,1}};
  dInt match;
  dErr err;

  dFunctionBegin;
  switch (topo) {
    case iMesh_HEXAHEDRON:
      for (dInt i=0; i<6; i++) {
        err = dGeomOrientFindPerm_HexQuad(econn,&conn[i*4],i,&match);dCHK(err);
        deg[ind[i]*3+0] = dMin(deg[ind[i]*3+0],edeg[hexperm[i][orient[match][0]]]);
        deg[ind[i]*3+1] = dMin(deg[ind[i]*3+1],edeg[hexperm[i][orient[match][1]]]);
        deg[ind[i]*3+2] = 1;
      }
      break;
    case iMesh_QUADRILATERAL:
      for (dInt i=0; i<4; i++) {
        match = dGeomMatchQuadLine(econn,&conn[i*2]);
        if (match < 0) dERROR(1,"Connectivity for edge %d does not match this Quad",i);
        if (match != i) dERROR(1,"Weird, edges aren't in order");
        deg[ind[i]*3+0] = dMin(deg[ind[i]*3+0],edeg[quadperm[match]]);
        deg[ind[i]*3+1] = deg[ind[i]*3+2] = 1;
      }
      break;
    case iMesh_LINE_SEGMENT:
      for (dInt i=0; i<2; i++) { /* Both endpoints are vertices, they always have degree 1 */
        if (!(conn[i] == econn[0] || conn[i] == econn[1])) dERROR(1,"Looks like incorrect connectivity");
        deg[ind[i]*3+0] = deg[ind[i]*3+1] = deg[ind[i]*3+2] = 1;
      }
      break;
    default: dERROR(1,"Topology %s not supported",iMesh_TopologyName[topo]);
  }
  dFunctionReturn(0);
}

static dErr dRealTableView(dInt m,dInt n,const dReal mat[],const char *name,dViewer viewer)
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  for (dInt i=0; i<m; i++) {
    if (name) {
      err = PetscViewerASCIIPrintf(viewer,"%10s[%2d][%2d:%2d] ",name,i,0,n);dCHK(err);
    }
    err = PetscViewerASCIIUseTabs(viewer,PETSC_NO);dCHK(err);
    for (dInt j=0; j<n; j++) {
      err = PetscViewerASCIIPrintf(viewer," % 9.5f",mat[i*n+j]);dCHK(err);
    }
    err = PetscViewerASCIIPrintf(viewer,"\n");dCHK(err);
    err = PetscViewerASCIIUseTabs(viewer,PETSC_YES);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr TensorJacobiHasBasis(dJacobi jac,dInt rule,dInt basis,dBool *has)
{
  Tensor this = jac->impl;

  dFunctionBegin;
  *has = (1 < basis && basis < this->N && basis <= rule && rule < this->M);
  dFunctionReturn(0);
}

/**
* Get the dRule corresponding to the given sizes.
*
* @param jac Jacobi context
* @param top topology
* @param rsize rule size in each dimension (should become polynomial degree to integrate exactly)
* @param rule space to write the rule
* @param base base of array to write the private data, if NULL, write nothing
* @param[in,out] index into \a data to write the private data, incremented to point to the next free space
*
* @return
*/
static dErr dJacobiGetRule_Tensor(dJacobi jac,dEntTopology top,const dInt rsize[],dRule *ruleout,void **base,dInt *index)
{
  Tensor this = (Tensor)jac->impl;
  void **start = &base[*index];
  struct s_dRule_Tensor_Line *line;
  struct s_dRule_Tensor_Quad *quad;
  struct s_dRule_Tensor_Hex *hex;
  dErr err;

  dFunctionBegin;
  switch (top) {
    case iMesh_LINE_SEGMENT:
      if (base) {
        line = (struct s_dRule_Tensor_Line*)start;
        line->ops = this->ruleOpsLine;
        err = TensorGetRule(this,rsize[0],&line->trule[0]);dCHK(err);
        *ruleout = (dRule)line;
      }
      if (index) *index += (dInt)sizeof(*line)/(dInt)sizeof(base[0]);
      break;
    case iMesh_QUADRILATERAL:
      if (base) {
        quad = (struct s_dRule_Tensor_Quad*)start;
        quad->ops = this->ruleOpsQuad;
        err = TensorGetRule(this,rsize[0],&quad->trule[0]);dCHK(err);
        err = TensorGetRule(this,rsize[1],&quad->trule[1]);dCHK(err);
        *ruleout = (dRule)quad;
      }
      if (index) *index += (dInt)sizeof(*quad)/(dInt)sizeof(base[0]);
      break;
    case iMesh_HEXAHEDRON:
      if (base) {
        hex = (struct s_dRule_Tensor_Hex*)start;
        hex->ops = this->ruleOpsHex;
        err = TensorGetRule(this,rsize[0],&hex->trule[0]);dCHK(err);
        err = TensorGetRule(this,rsize[1],&hex->trule[1]);dCHK(err);
        err = TensorGetRule(this,rsize[2],&hex->trule[2]);dCHK(err);
        *ruleout = (dRule)hex;
      }
      if (index) *index += (dInt)sizeof(*hex)/(dInt)sizeof(base[0]);
      break;
    default:
      dERROR(1,"no rule available for given topology");
  }
  dFunctionReturn(0);
}

/**
* Fill in an EFS of the specified order.
*
* @param jac context
* @param top topology (from the iMesh_Topology enum)
* @param bsize vector of basis sizes in each direction
* @param rule rule on which the dEFS resides
* @param efs dEFS context, must be allocated
* @param base start of array in which to put private data, or NULL to only calculate private space requirement
* @param index index into base to start putting private data, updated on exit
*
* @return
*/
dErr dJacobiGetEFS_Tensor(dJacobi jac,dEntTopology top,const dInt bsize[],dRule rule,dEFS *efsout,void **base,dInt *index)
{
  Tensor this = (Tensor)jac->impl;
  void **start = &base[*index];
  struct s_dEFS_Tensor_Line *line;
  struct s_dEFS_Tensor_Quad *quad;
  struct s_dEFS_Tensor_Hex *hex;
  dInt rdim,rsize[3];
  dErr err;

  dFunctionBegin;
  err = dRuleGetTensorNodeWeight(rule,&rdim,rsize,NULL,NULL);dCHK(err);
  switch (top) {
    case iMesh_LINE_SEGMENT:
      if (rdim != 1) dERROR(1,"Incompatible Rule size %d, expected 1",rdim);
      if (base) {
        line = (struct s_dEFS_Tensor_Line*)start;
        line->ops = this->efsOpsLine;
        line->rule = rule;
        err = TensorGetBasis(this,rsize[0],bsize[0],&line->basis[0]);dCHK(err);
        *efsout = (dEFS)line;
      }
      if (index) *index += (dInt)sizeof(*line)/(dInt)sizeof(base[0]);
      break;
    case iMesh_QUADRILATERAL:
      if (rdim != 2) dERROR(1,"Incompatible Rule size %d, expected 2",rdim);
      if (base) {
        quad = (struct s_dEFS_Tensor_Quad*)start;
        quad->ops = this->efsOpsQuad;
        quad->rule = rule;
        for (dInt i=0; i<2; i++) {
          err = TensorGetBasis(this,rsize[i],bsize[i],&quad->basis[i]);dCHK(err);
        }
        *efsout = (dEFS)quad;
      }
      if (index) *index += (dInt)sizeof(*quad)/(dInt)sizeof(base[0]);
      break;
    case iMesh_HEXAHEDRON:
      if (rdim != 3) dERROR(1,"Incompatible Rule size %d, expected 3",rdim);
      if (base) {
        hex = (struct s_dEFS_Tensor_Hex*)start;
        hex->ops = this->efsOpsHex;
        hex->rule = rule;
        for (dInt i=0; i<3; i++) {
          err = TensorGetBasis(this,rsize[i],bsize[i],&hex->basis[i]);dCHK(err);
        }
        *efsout = (dEFS)hex;
      }
      if (index) *index += (dInt)sizeof(*hex)/(dInt)sizeof(base[0]);
      break;
    default:
      dERROR(1,"no basis available for given topology");
  }
  dFunctionReturn(0);
}


static dErr TensorBuilderCreate(void *opts,TensorBuilder *out)
{
  dErr err;
  TensorBuilder new;

  dFunctionBegin;
  err = dNew(struct s_TensorBuilder,&new);dCHK(err);
  new->options = opts;
  new->workLength = 0;
  new->work = NULL;
  *out = new;
  dFunctionReturn(0);
}

static dErr TensorBuilderGetArray(TensorBuilder build,dInt size,dReal **a)
{
  dErr err;

  dFunctionBegin;
  if (build->workLength < size) { /* make sure we have enough work space, we're currently not using this space. */
    err = dFree(build->work);dCHK(err);
    build->workLength = 2*size;
    err = dMalloc(build->workLength*sizeof(dReal),&build->work);dCHK(err);
  }
  *a = build->work;
  dFunctionReturn(0);
}

static dErr TensorBuilderDestroy(TensorBuilder build)
{
  dErr err;

  dFunctionBegin;
  err = dFree(build->work);dCHK(err);
  err = dFree(build);dCHK(err);
  dFunctionReturn(0);
}

static dErr TensorRuleCreate(TensorBuilder build,dInt size,TensorRule *rule)
{
  TensorRuleOptions opt = (TensorRuleOptions)build->options;
  dReal *work;
  TensorRule r;
  dErr err;

  dFunctionBegin;
  if (!(0 < size && size < 50)) dERROR(1,"rule size out of bounds.");
  /* Check options */
  if (!opt) dERROR(1,"TensorRuleOptions not set.");
  if (opt->family != GAUSS) dERROR(1,"GaussFamily %d not supported",opt->family);
  err = dNew(struct s_TensorRule,rule);dCHK(err);

  r = *rule;
  err = PetscMalloc2(size,dReal,&r->weight,size,dReal,&r->coord);dCHK(err);
  err = TensorBuilderGetArray(build,size,&work);dCHK(err); /* We're not using this now */

  r->size = size;
  if (size == 1) {              /* Polylib function fails for this size. */
    r->weight[0] = 2.0;
    r->coord[0] = 0.0;
  } else {
    zwgj(r->coord,r->weight,size,opt->alpha,opt->beta); /* polylib function */
  }
  dFunctionReturn(0);
}

static dErr TensorRuleDestroy(TensorRule rule)
{
  dErr err;

  dFunctionBegin;
  if (!rule) dFunctionReturn(0);
  err = PetscFree2(rule->weight,rule->coord);dCHK(err);
  err = dFree(rule);dCHK(err);
  dFunctionReturn(0);
}

dErr TensorRuleView(const TensorRule rule,PetscViewer viewer) /* exported so that topology implementation can use it */
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  err = PetscViewerASCIIPrintf(viewer,"TensorRule with %d nodes.\n",rule->size);dCHK(err);
  err = dRealTableView(1,rule->size,rule->coord,"q",viewer);dCHK(err);
  err = dRealTableView(1,rule->size,rule->weight,"w",viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr TensorBasisCreate(TensorBuilder build,const TensorRule rule,dInt P,TensorBasis *basis)
{
  TensorBasisOptions opt = (TensorBasisOptions)build->options;
  TensorBasis b;
  const dInt Q=rule->size;
  dReal *work;
  dErr err;

  dFunctionBegin;
  if (!(0 < P && P <= Q)) dERROR(1,"Requested TensorBasis size %d out of bounds.",P);
  if (!opt) dERROR(1,"TensorRuleOptions not set.");
  if (opt->family != GAUSS_LOBATTO) dERROR(1,"GaussFamily %d not supported",opt->family);
  err = dNew(struct s_TensorBasis,&b);dCHK(err);
  err = PetscMalloc3(P*Q,dReal,&b->interp,P*Q,dReal,&b->deriv,P,dReal,&b->node);dCHK(err);
  err = TensorBuilderGetArray(build,2*P*Q,&work);dCHK(err); /* We're not using this now */
  b->Q = Q;
  b->P = P;

  if (P == 1) {         /* degenerate case */
    b->interp[0] = 1.0;
    b->deriv[0] = 0.0;
    b->node[0] = 0.0;
    dFunctionReturn(0);
  } else {
    const dReal alpha=opt->alpha,beta=opt->beta;
    dReal *cDeriv = work;        /* collocation derivative at Gauss-Lobatto points */
    dReal *cDerivT = work + P*Q; /* useless matrix spewed out of Dglj() */
    dReal *interp = b->interp;
    dReal *node = b->node;
    node[0] = -1.0; node[P-1] = 1.0;
    jacobz(P-2,node+1,alpha+1.0,beta+1.0);        /* Gauss-Lobatto nodes */
    Imglj(interp,node,rule->coord,P,Q,alpha,beta); /* interpolation matrix */
    Dglj(cDeriv,cDerivT,node,P,alpha,beta);       /* collocation derivative matrix */
    for (dInt i=0; i<Q; i++) {
      for (dInt j=0; j<P; j++) {
        dReal z = 0;
        for (dInt k=0; k<P; k++) {
          z += interp[i*P+k] * cDeriv[k*P+j];
        }
        b->deriv[i*P+j] = z;
      }
    }
  }
  *basis = b;
  dFunctionReturn(0);
}

static dErr TensorBasisDestroy(TensorBasis basis)
{
  dErr err;

  dFunctionBegin;
  if (!basis) dFunctionReturn(0);
  err = PetscFree3(basis->interp,basis->deriv,basis->node);dCHK(err);
  err = dFree(basis);dCHK(err);
  dFunctionReturn(0);
}

dErr TensorBasisView(const TensorBasis basis,PetscViewer viewer) /* exported so that topology implementations can see it */
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  err = PetscViewerASCIIPrintf(viewer,"TensorBasis with rule=%d basis=%d.\n",basis->Q,basis->P);dCHK(err);
  err = dRealTableView(basis->Q,basis->P,basis->interp,"interp",viewer);dCHK(err);
  err = dRealTableView(basis->Q,basis->P,basis->deriv,"deriv",viewer);dCHK(err);
  err = dRealTableView(1,basis->P,basis->node,"node",viewer);dCHK(err);
  dFunctionReturn(0);
}


/**
* Just an error checking indexing function.
*
* @param this
* @param m
* @param out
*
* @return
*/
static dErr TensorGetRule(Tensor this,dInt m,TensorRule *out)
{

  dFunctionBegin;
  if (!this->setupcalled) dERROR(1,"Attempt to get rule before Tensor setup.");
  if (!(0 < m && m < this->M)) dERROR(1,"Rule %d not less than limit %d",m,this->M);
  *out = this->rule[m];
  dFunctionReturn(0);
}

/**
* An error checking lookup function.
*
* @param this
* @param m
* @param n
* @param out
*
* @return
*/
static dErr TensorGetBasis(Tensor this,dInt m,dInt n,TensorBasis *out)
{

  dFunctionBegin;
  if (!this->setupcalled) dERROR(1,"Attempt to get basis before Tensor setup.");
  if (!(0 < m && m < this->M)) dERROR(1,"Rule size %d not less than limit %d",m,this->M);
  if (!(0 < n && n < this->N)) dERROR(1,"Basis size %d not less than limit %d",n,this->N);
  *out = this->basis[m*this->N+n];
  dFunctionReturn(0);
}
