#include <dohpfsimpl.h>
#include <stdarg.h>

/** This is a fairly high-level function to get the preferred quadrature for this FS.  Usually the FS with the highest
* order elements, or the physics with the most aliasing problems, is used to define the quadrature on which all
* components are integrated.
*
* @note This function helps to hide the low-level dQuadrature object from the user, since it is almost always the case
* that the user wants the "best" quadrature for a particular function space they are working with.
**/
dErr dFSGetPreferredQuadratureRuleSet(dFS fs,dMeshESH set,dEntType etype,dEntTopology etopo,dQuadratureMethod method,dRuleset *ruleset)
{
  dInt             ents_a,ents_s;
  dEntTopology     *topo;
  dPolynomialOrder *order;
  dQuadrature      quad;
  dRuleset         rset;
  dMeshEH          *ents;
  dErr             err;

  dFunctionBegin;
  *ruleset = NULL;
  err = dMeshGetNumEnts(fs->mesh,set,etype,etopo,&ents_a);dCHK(err);
  err = dMallocA3(ents_a,&ents,ents_a,&topo,ents_a,&order);dCHK(err);
  err = dMeshGetEnts(fs->mesh,set,etype,etopo,ents,ents_a,&ents_s);dCHK(err);
  err = dMeshGetTopo(fs->mesh,ents_s,ents,topo);dCHK(err);
  err = dMeshTagGetData(fs->mesh,fs->tag.degree,ents,ents_s,order,ents_s,dDATA_INT);dCHK(err);

  /* Request exact integration of a mass matrix.  More generally, the required order should be based on the order of the
   * adjacent elements, and perhaps also the physics.
   */
  for (dInt i=0; i<ents_s; i++) {
    order[i] = dPolynomialOrderCreate(2*(dPolynomialOrderMax(order[i])),
                                      2*(dPolynomialOrder1D(order[i],0)),
                                      2*(dPolynomialOrder1D(order[i],1)),
                                      2*(dPolynomialOrder1D(order[i],2)));
  }

  err = dNew(struct _n_dRuleset,&rset);dCHK(err);
  err = dFSGetMesh(fs,&rset->mesh);dCHK(err);
  rset->set = set;
  rset->type = etype;
  rset->topo = etopo;
  rset->n = ents_s;
  err = dJacobiGetQuadrature(fs->jacobi,method,&quad);dCHK(err);
  err = dQuadratureGetRules(quad,ents_s,topo,order,&rset->rules);dCHK(err);
  err = dFree3(ents,topo,order);dCHK(err);
  *ruleset = rset;
  dFunctionReturn(0);
}

static dErr dRulesetWorkspaceDestroy(struct dRulesetWorkspace *ws)
{
  dErr err;

  dFunctionBegin;
  if (!ws) dFunctionReturn(0);
  err = dFree4(ws->q,ws->cjac,ws->cjinv,ws->jw);dCHK(err);
  for (struct dRulesetWorkspaceLink *link = ws->link,*next; link; link = next) {
    err = dFree4(link->u,link->v,link->du,link->dv);dCHK(err);
    next = link->next;
    err = dFree(link);dCHK(err);
  }
  err = dFree(ws);dCHK(err);
  dFunctionReturn(0);
}


dErr dRulesetDestroy(dRuleset rset)
{
  dErr err;

  dFunctionBegin;
  err = dFree(rset->rules);dCHK(err);
  err = dRulesetWorkspaceDestroy(rset->workspace);dCHK(err);
  err = dFree(rset);dCHK(err);
  dFunctionReturn(0);
}

dErr dRulesetGetMaxQ(dRuleset rset,dInt *maxQ,dInt *maxnpatches,dInt *maxQelem)
{
  dErr err;

  dFunctionBegin;
  if (!rset->maxQ) {
    for (dInt i=0,Q,npatches; i<rset->n; i++) {
      err = dRuleGetSize(rset->rules[i],NULL,&Q);dCHK(err);
      err = dRuleGetPatches(rset->rules[i],&npatches,NULL,NULL,NULL);dCHK(err);
      if (Q % npatches) dERROR(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Number of quadrature points in element %D not divisible by number of patches %D",Q,npatches);
      rset->maxQ = dMaxInt(rset->maxQ,Q/npatches);dCHK(err);
      rset->maxnpatches = dMaxInt(rset->maxnpatches,npatches);
      rset->maxQelem = dMaxInt(rset->maxQelem,Q);
    }
  }
  if (maxQ) *maxQ = rset->maxQ;
  if (maxnpatches) *maxnpatches = rset->maxnpatches;
  if (maxQelem) *maxQelem = rset->maxQelem;
  dFunctionReturn(0);
}

dErr dRulesetGetSize(dRuleset rset,dInt *size)
{
  dFunctionBegin;
  *size = rset->n;
  dFunctionReturn(0);
}

/**
 * @note pass dof=0 to terminate
 */
dErr dRulesetGetWorkspace(dRuleset rset,dScalar **q,dScalar **cjac,dScalar **cjinv,dScalar **jw,dInt dof,...)
{
  dErr err;
  dInt Q;
  va_list ap;
  struct dRulesetWorkspace *ws;

  dFunctionBegin;
  err = dRulesetGetMaxQ(rset,&Q,NULL,NULL);dCHK(err);
  if (!rset->workspace) {
    err = dCallocA(1,&rset->workspace);dCHK(err);
  }
  ws = rset->workspace;
  if (!ws->q) {
    err = dMallocA4(Q*3,&ws->q,Q*9,&ws->cjac,Q*9,&ws->cjinv,Q*1,&ws->jw);dCHK(err);
  }
  if (q) *q         = ws->q;
  if (cjac) *cjac   = ws->cjac;
  if (cjinv) *cjinv = ws->cjinv;
  if (jw) *jw       = ws->jw;

  va_start(ap,dof);
  for (struct dRulesetWorkspaceLink *next,**linkp=&ws->link; dof; *linkp = next,linkp = &next->next) {
    dScalar **u,**v,**du,**dv;
    next = *linkp;
    if (!next) {
      err = dCallocA(1,&next);dCHK(err);
      err = dMallocA4(dof*Q,&next->u,dof*Q,&next->v,dof*dof*Q,&next->du,dof*dof*Q,&next->dv);dCHK(err);
      next->dof = dof;
      next->next = NULL;
    }
    if (next->dof != dof) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"changing size of requested link");
    if (next->checkedout) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"already checked out");
    next->checkedout = dTRUE;
    u  = va_arg(ap,dScalar**);
    v  = va_arg(ap,dScalar**);
    du = va_arg(ap,dScalar**);
    dv = va_arg(ap,dScalar**);
    if (u) *u   = next->u;
    if (v) *v   = next->v;
    if (du) *du = next->du;
    if (dv) *dv = next->dv;
    dof = va_arg(ap,dInt);
  }
  va_end(ap);
  dFunctionReturn(0);
}

/** dRulesetRestoreWorkspace - return a workspace managed by the ruleset
 *
 */
dErr dRulesetRestoreWorkspace(dRuleset rset,dScalar **q,dScalar **cjac,dScalar **cjinv,dScalar **jw,dInt dof,...)
{
  va_list ap;

  dFunctionBegin;
  for (struct dRulesetWorkspaceLink *next = rset->workspace->link; next; next = next->next) {
    next->checkedout = dFALSE;
  }
  if (q) *q         = NULL;
  if (cjac) *cjac   = NULL;
  if (cjinv) *cjinv = NULL;
  if (jw) *jw       = NULL;
  va_start(ap,dof);
  for (; dof; dof = va_arg(ap,dInt)) {
    dScalar **u,**v,**du,**dv;
    u  = va_arg(ap,dScalar**);
    v  = va_arg(ap,dScalar**);
    du = va_arg(ap,dScalar**);
    dv = va_arg(ap,dScalar**);
    if (u) *u   = NULL;
    if (v) *v   = NULL;
    if (du) *du = NULL;
    if (dv) *dv = NULL;
  }
  va_end(ap);
  dFunctionReturn(0);
}
