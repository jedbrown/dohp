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
static dErr dJacobiGetRule_Tensor(dJacobi jac,dEntTopology top,const dInt rsize[],dRule *rule,void **base,dInt *index);
static dErr dJacobiGetEFS_Tensor(dJacobi jac,dEntTopology top,const dInt bsize[],dRule rule,dEFS *efs,void **base,dInt *index);

static dErr TensorGetRule(Tensor this,dInt n,TensorRule *out);
static dErr TensorGetBasis(Tensor this,dInt m,dInt n,TensorBasis *out);

static dErr TensorJacobiHasBasis(dJacobi,dInt,dInt,dBool*);

static dErr dRealTableView(dInt m,dInt n,const dReal mat[],const char name[],dViewer viewer);

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

static dErr dJacobiPropogateDown_Tensor(dUNUSED dJacobi jac,const struct dMeshAdjacency *a,dInt deg[])
{
  static const dInt quadperm[4] = {0,1,0,1};
  static const dInt hexperm[6][2] = {{0,2},{1,2},{0,2},{1,2},{1,0},{0,1}}; /* map natural axis of Quad to natural axis of Hex */
  static const dInt orient[8][2] = {{0,1},{1,0},{0,1},{1,0},            /* map actual Quad orientation to proper orientation */
                                    {1,0},{0,1},{1,0},{0,1}};
  dInt e,i,ai,aind,match;
  dErr err;

  dFunctionBegin;
  for (e=a->toff[dTYPE_REGION]; e<a->toff[dTYPE_REGION+1]; e++) {
    switch (a->topo[e]) {
      case iMesh_HEXAHEDRON:
        for (i=0; i<6; i++) {
          ai = a->adjoff[e]+i; aind = a->adjind[ai]; match = a->adjperm[ai];
          deg[aind*3+0] = dMin(deg[aind*3+0],deg[e*3+hexperm[i][orient[match][0]]]);
          deg[aind*3+1] = dMin(deg[aind*3+1],deg[e*3+hexperm[i][orient[match][1]]]);
          deg[aind*3+2] = 1;
          if (a->topo[aind] != dTOPO_QUAD) dERROR(1,"corrupt adjacency");
        }
        break;
      default: dERROR(1,"Region topology %d not supported",a->topo[e]);dCHK(err);
    }
  }
  for (e=a->toff[dTYPE_FACE]; e<a->toff[dTYPE_FACE+1]; e++) {
    switch (a->topo[e]) {
      case iMesh_QUADRILATERAL:
        for (i=0; i<4; i++) {
          ai = a->adjoff[e]+i; aind = a->adjind[ai];
          deg[aind*3+0] = dMin(deg[aind*3+0],deg[e*3+quadperm[i]]);
          deg[aind*3+1] = deg[aind*3+2] = 1;
          if (a->topo[aind] != dTOPO_LINE) dERROR(1,"corrupt adjacency");
        }
        break;
      default: dERROR(1,"Face topology %d not supported",a->topo[e]);dCHK(err);
    }
  }
  for (e=a->toff[dTYPE_EDGE]; e<a->toff[dTYPE_EDGE+1]; e++) {
    for (i=0; i<2; i++) {
      ai = a->adjoff[e]+i; aind = a->adjind[ai];
      deg[aind*3+0] = deg[aind*3+1] = deg[aind*3+2] = 1; /* should always be vertices */
      if (a->topo[aind] != dTOPO_POINT) dERROR(1,"corrupt adjacency");
    }
  }
  dFunctionReturn(0);
}

static dErr dJacobiGetNodeCount_Tensor(dUNUSED dJacobi jac,dInt count,const dEntTopology top[],const dInt deg[],dInt inode[],dInt xnode[])
{
  dInt ideg[3];
  const dInt *xdeg;

  dFunctionBegin;
  for (dInt i=0; i<count; i++) {
    for (dInt j=0; j<3; j++) {
      ideg[j] = dMax(0,deg[3*i+j]-2);
    }
    xdeg = deg+3*i;
    switch (top[i]) {
      case dTOPO_HEX:
        if (inode) inode[i] = ideg[0]*ideg[1]*ideg[2];
        if (xnode) xnode[i] = xdeg[0]*xdeg[1]*xdeg[2];
        break;
      case dTOPO_QUAD:
        if (inode) inode[i] = ideg[0]*ideg[1];
        if (xnode) xnode[i] = xdeg[0]*xdeg[1];
        break;
      case dTOPO_LINE:
        if (inode) inode[i] = ideg[0];
        if (xnode) xnode[i] = xdeg[0];
        break;
      case dTOPO_POINT:
        if (inode) inode[i] = 1;
        if (xnode) xnode[i] = 1;
        break;
      default:
        dERROR(1,"Topology %d not supported",top[i]);
    }
  }
  dFunctionReturn(0);
}

static dErr dJacobiGetConstraintCount_Tensor(dUNUSED dJacobi jac,dInt nx,const dInt xi[],const dUNUSED dInt xs[],const dInt deg[],const struct dMeshAdjacency *ma,dInt nnz[],dInt pnnz[])
{
  dInt i,j;

  dFunctionBegin;
  for (i=0; i<nx; i++) {
    j = xi[i];
    if (0) {
      switch (ma->topo[i]) {
        case dTOPO_HEX:
          for (dInt i0=1; i0<deg[3*j]-1; i0++) {
          }
          if (1  ) {          /* Conforming */
            nnz[i] = pnnz[i] = 1;
          } else {                /* Nonconforming */

          }
          break;
        default: dERROR(1,"not implemented for expanded topology %d",ma->topo[i]);
      }
    } else {
      for (j=xs[i]; j<xs[i+1]; j++) {
        nnz[j] = pnnz[j] = 1;
      }
    }
  }
  dFunctionReturn(0);
}

#define ASSERT(cond) if (!(cond)) dERROR(1,"assert " #cond "failed")

static inline dInt same3(dInt a,dInt b,dInt c)
{
  dFunctionBegin;
  ASSERT(a == b);
  ASSERT(a == c);
  dFunctionReturn(a);
}

static inline dInt same2(dInt a,dInt b)
{
  dFunctionBegin;
  ASSERT(a == b);
  dFunctionReturn(a);
}

/** Get the index into the interior of a face with orientation \a perm, (full) dimensions \a dim, and index \a ij
*
* @param perm Permutation (0,...,7)
* @param dim  Full dimension of face in canonical orientation, the interior dimensions are two less
* @param ij   Interior index in canonical orientation
* @param[out] ind  Index into interior of face
*/
static inline dErr dGeomPermQuadIndex(dInt perm,const dInt dim[],const dInt ij[2],dInt *ind)
{
  const dInt M = dim[0]-2,N = dim[1]-2,i = ij[0],j = ij[1];

  dFunctionBegin;
  switch (perm) {
    case 0: *ind = i*N + j; break;
    case 1: *ind = (M-1-i) + j*M; break;
    case 2: *ind = (M-1-i)*N + (N-1-j); break;
    case 3: *ind = i + (N-1-j)*M; break;
    case 4: *ind = i + j*M; break;
    case 5: *ind = i*N + (N-1-j); break;
    case 6: *ind = (M-1-i) + (N-1-j)*M; break;
    case 7: *ind = (M-1-i)*N + j; break;
    default: dERROR(1,"Invalid permutation");
  }
  dFunctionReturn(0);
}

static dErr dJacobiAddConstraints_Tensor(dJacobi dUNUSED jac,dInt nx,const dInt xi[],const dInt xs[],const dInt is[],const dInt deg[],const struct dMeshAdjacency *ma,Mat C,Mat Cp)
{
  dInt elem,i,j,k,l,e[12],eP[12],v[8];
  dInt nrow,ncol,irow[10],icol[30];
  dScalar interp[30]; /* FIXME: check bounds */
  const dInt *f,*fP;
  const dInt *aI=ma->adjind,*aO=ma->adjoff,*aP=ma->adjperm;
  dErr err;

  dFunctionBegin;
  for (elem=0; elem<nx; elem++) {
    const dInt ei = xi[elem]; /* Element index, \a deg and everything in \a ma is addressed by \a ei. */
    const dInt *d = deg+3*ei,d0 = d[0],d1 = d[1],d2 = d[2];
    const dInt scan[3] = {d1*d2,d2,1};
    switch (ma->topo[ei]) {
      case dTOPO_HEX:                              /* ****************************** HEX ************************** */
        f = aI+aO[ei];         /* Array of indices for adjacent faces */
        fP = aP+aO[ei];        /* Permutations for faces */
        for (i=0; i<12; i++) { /* Extract and orient edges */
          const struct {dInt f,e,o;} fe[12][2] = { /* Each edge is adjacent to two faces, f=face, e=edge_of_face, o=orientation */
            {{0,0,0},{4,3,1}}, {{1,0,0},{4,2,1}}, {{2,0,0},{4,1,1}}, {{3,0,0},{4,0,1}},  /* bottom */
            {{0,2,1},{5,0,0}}, {{1,2,1},{5,1,0}}, {{2,2,1},{5,2,0}}, {{3,2,1},{5,3,0}},  /* top */
            {{0,3,1},{3,1,0}}, {{0,1,0},{1,3,1}}, {{1,1,0},{2,3,1}}, {{2,1,0},{3,3,1}}}; /* vertical */
          //const dInt fperm[8][4] = {{0,1,2,3},{1,2,3,0},{2,3,0,1},{3,0,1,2},{3,2,1,0},{0,3,2,1},{1,0,3,2},{0,3,2,1}};
          const dInt iperm[8][4] = {{0,1,2,3},{3,0,1,2},{2,3,0,1},{1,2,3,0},{3,2,1,0},{0,3,2,1},{1,0,3,2},{2,1,0,3}};
          const dInt ef0 = fe[i][0].f,ef1 = fe[i][1].f;
          e[i] = aI[aO[f[ef0]] + iperm[fP[ef0]][fe[i][0].e]];
          if (e[i] != aI[aO[f[ef1]] + iperm[fP[ef1]][fe[i][1].e]])
            dERROR(1,"faces don't agree about edge");
          eP[i] = aP[aO[f[ef0]] + iperm[fP[ef0]][fe[i][0].e]] ^ (((fP[ef0]>>2) & 1) ^ fe[i][0].o);
          if (eP[i] != (aP[aO[f[ef1]] + iperm[fP[ef1]][fe[i][1].e]] ^ (((fP[ef1]>>2) & 1) ^ fe[i][1].o)))
            dERROR(1,"orientations do not agree");
        }
        for (i=0; i<8; i++) { /* Extract vertices */
          const dInt edge_flip[2][2] = {{0,1},{1,0}};
          const struct {dInt e,v;} ev_common[8][3] = {
            {{0,0},{3,1},{8,0}}, {{0,1},{1,0},{9,0}}, {{1,1},{2,0},{10,0}}, {{2,1},{3,0},{11,0}},
            {{4,0},{7,1},{8,1}}, {{4,1},{5,0},{9,1}}, {{5,1},{6,0},{10,1}}, {{6,1},{7,0},{11,1}}};
#define E(j) ev_common[i][j].e
#define EV(j) aI[aO[e[E(j)]] + edge_flip[eP[E(j)]][ev_common[i][j].v]]
          v[i] = EV(0);
          if (v[i] != EV(1)) dERROR(1,"first two edges don't agree about vertex");
          if (v[i] != EV(2)) dERROR(1,"first and third edges don't agree about vertex");
#undef E
#undef EV
        }
        for (i=0; i<8; i++) { /* Set vertices, always conforming until we have h-nonconforming meshes */
          const dInt T[8][3] = {{0,0,0},{d0-1,0,0},{d0-1,d1-1,0},{0,d1-1,0},
                                {0,0,d2-1},{d0-1,0,d2-1},{d0-1,d1-1,d2-1},{0,d1-1,d2-1}};
          err = MatSetValue(C, xs[elem]+(T[i][0]*d1+T[i][1])*d2+T[i][2],is[v[i]],1,INSERT_VALUES);dCHK(err);
          err = MatSetValue(Cp,xs[elem]+(T[i][0]*d1+T[i][1])*d2+T[i][2],is[v[i]],1,INSERT_VALUES);dCHK(err);
        }
        for (i=0; i<12; i++) { /* Set edges */
          const struct {dInt start[3],incd,inci,end;} E[12] = { /* How to traverse the interior of the edge in forward order */
            {{1,0,0},0,1,d0-1}, {{d0-1,1,0},1,1,d1-1}, {{d0-2,d1-1,0},0,-1,0}, {{0,d1-2,0},1,-1,0},
            {{1,0,d2-1},0,1,d0-1}, {{d0-1,1,d2-1},1,1,d1-1}, {{d0-2,d1-1,d2-1},0,-1,0}, {{0,d1-2,d2-1},1,-1,0},
            {{0,0,1},2,1,d2-1}, {{d0-1,0,1},2,1,d2-1}, {{d0-1,d1-1,1},2,1,d2-1}, {{0,d1-1,1},2,1,d2-1}};
          const dInt *start = E[i].start,incd = E[i].incd,inci = E[i].inci,end = E[i].end;
          if (deg[3*e[i]] != deg[3*ei+E[i].incd]) dERROR(1,"degree does not agree, p-nonconforming");
          for (j=start[incd]; j!=end; j += inci) {
            nrow = 0; ncol = 0;
            irow[nrow++] = xs[elem] + (start[0]*d1+start[1])*d2+start[2] + (j-start[incd])*scan[incd];
            switch (eP[i]) {
              case 0: icol[ncol++] = is[e[i]] + (j-start[incd])/inci; break; /* traverse the edge forwards */
              case 1: icol[ncol++] = is[e[i]] - (j-(end-inci))/inci; break;  /* traverse the edge in reverse */
            }
            interp[0] = 1;
            err = MatSetValues(C,nrow,irow,ncol,icol,interp,INSERT_VALUES);dCHK(err);
            err = MatSetValues(Cp,nrow,irow,ncol,icol,interp,INSERT_VALUES);dCHK(err);
          }
        }
        for (i=0; i<6; i++) { /* Faces */
          const struct {dInt start[3],incd[2],inc[2],end[2];} F[6] = {
            {{1,0,1},{0,2},{1,1},{d0-1,d2-1}},     {{d0-1,1,1},{1,2},{1,1},{d1-1,d2-1}},
            {{d0-2,d1-1,1},{0,2},{-1,1},{0,d2-1}}, {{0,d1-2,1},{1,2},{-1,1},{0,d2-1}},
            {{1,1,0},{1,0},{1,1},{d1-1,d0-1}},     {{1,1,d2-1},{0,1},{1,1},{d0-1,d1-1}}};
          const dInt *start = F[i].start,*incd = F[i].incd,*inc = F[i].inc,*end = F[i].end;
          for (j=start[incd[0]]; j!=end[0]; j+=inc[0]) {
            for (k=start[incd[1]]; k!=end[1]; k+=inc[1]) {
              const dInt faceDim[2] = {d[incd[0]],d[incd[1]]};
              dInt facejk[2],faceIndex;
              nrow = ncol = 0;
              irow[nrow++] = (xs[elem] + (start[0]*d1+start[1])*d2+start[2]
                              + (j-start[incd[0]])*scan[incd[0]] + (k-start[incd[1]])*scan[incd[1]]);
              facejk[0] = (j-start[incd[0]])/inc[0];
              facejk[1] = (k-start[incd[1]])/inc[1];
              err = dGeomPermQuadIndex(fP[i],faceDim,facejk,&faceIndex);dCHK(err);
              icol[ncol++] = is[f[i]] + faceIndex;
              interp[0] = 1;
              err = MatSetValues(C, nrow,irow,ncol,icol,interp,INSERT_VALUES);dCHK(err);
              err = MatSetValues(Cp,nrow,irow,ncol,icol,interp,INSERT_VALUES);dCHK(err);
            }
          }
        }
        /* Set trivial map for region interior dofs */
        for (j=1; j<d0-1; j++) {
          for (k=1; k<d1-1; k++) {
            for (l=1; l<d2-1; l++) {
              irow[0] = xs[elem] + (j*d1+k)*d2+l;
              icol[0] = is[ei] + ((j-1)*(d1-2)+(k-1))*(d2-2)+(l-1);
              interp[0] = 1;
              err = MatSetValues(C, 1,irow,1,icol,interp,INSERT_VALUES);dCHK(err);
              err = MatSetValues(Cp,1,irow,1,icol,interp,INSERT_VALUES);dCHK(err);
            }
          }
        }
        break;
      default: dERROR(1,"not implemented for expanded topology %d",ma->topo[ei]);
    }
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
  Tensor this = jac->impl;
  void **start = &base[*index];
  dRule_Tensor *rule;
  dErr err;

  dFunctionBegin;
  switch (top) {
    case iMesh_LINE_SEGMENT:
      if (base) {
        rule = (dRule_Tensor*)start;
        rule->ops = this->ruleOpsLine;
        err = TensorGetRule(this,rsize[0],&rule->trule[0]);dCHK(err);
        *ruleout = (dRule)rule;
      }
      if (index) *index += (dInt)sizeof(*rule)/(dInt)sizeof(base[0]);
      break;
    case iMesh_QUADRILATERAL:
      if (base) {
        rule = (dRule_Tensor*)start;
        rule->ops = this->ruleOpsQuad;
        err = TensorGetRule(this,rsize[0],&rule->trule[0]);dCHK(err);
        err = TensorGetRule(this,rsize[1],&rule->trule[1]);dCHK(err);
        *ruleout = (dRule)rule;
      }
      if (index) *index += (dInt)sizeof(*rule)/(dInt)sizeof(base[0]);
      break;
    case iMesh_HEXAHEDRON:
      if (base) {
        rule = (dRule_Tensor*)start;
        rule->ops = this->ruleOpsHex;
        err = TensorGetRule(this,rsize[0],&rule->trule[0]);dCHK(err);
        err = TensorGetRule(this,rsize[1],&rule->trule[1]);dCHK(err);
        err = TensorGetRule(this,rsize[2],&rule->trule[2]);dCHK(err);
        *ruleout = (dRule)rule;
      }
      if (index) *index += (dInt)sizeof(*rule)/(dInt)sizeof(base[0]);
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
  Tensor this = jac->impl;
  void **start = &base[*index];
  dEFS_Tensor *efs;
  dInt rdim,rsize[3];
  dErr err;

  dFunctionBegin;
  err = dRuleGetTensorNodeWeight(rule,&rdim,rsize,NULL,NULL);dCHK(err);
  switch (top) {
    case iMesh_LINE_SEGMENT:
      if (rdim != 1) dERROR(1,"Incompatible Rule size %d, expected 1",rdim);
      if (base) {
        efs = (dEFS_Tensor*)start;
        efs->ops = this->efsOpsLine;
        efs->rule = rule;
        err = TensorGetBasis(this,rsize[0],bsize[0],&efs->basis[0]);dCHK(err);
        *efsout = (dEFS)efs;
      }
      if (index) *index += (dInt)sizeof(*efs)/(dInt)sizeof(base[0]);
      break;
    case iMesh_QUADRILATERAL:
      if (rdim != 2) dERROR(1,"Incompatible Rule size %d, expected 2",rdim);
      if (base) {
        efs = (dEFS_Tensor*)start;
        efs->ops = this->efsOpsQuad;
        efs->rule = rule;
        for (dInt i=0; i<2; i++) {
          err = TensorGetBasis(this,rsize[i],bsize[i],&efs->basis[i]);dCHK(err);
        }
        *efsout = (dEFS)efs;
      }
      if (index) *index += (dInt)sizeof(*efs)/(dInt)sizeof(base[0]);
      break;
    case iMesh_HEXAHEDRON:
      if (rdim != 3) dERROR(1,"Incompatible Rule size %d, expected 3",rdim);
      if (base) {
        efs = (dEFS_Tensor*)start;
        efs->ops = this->efsOpsHex;
        efs->rule = rule;
        for (dInt i=0; i<3; i++) {
          err = TensorGetBasis(this,rsize[i],bsize[i],&efs->basis[i]);dCHK(err);
        }
        *efsout = (dEFS)efs;
      }
      if (index) *index += (dInt)sizeof(*efs)/(dInt)sizeof(base[0]);
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



/**
* Initializes the ops table, this is the only non-static function in this file.
*
* @param jac
*
* @return
*/
dErr dJacobiCreate_Tensor(dJacobi jac)
{
  static const struct _dJacobiOps myops = {
    .SetUp              = dJacobiSetUp_Tensor,
    .SetFromOptions     = 0,
    .Destroy            = dJacobiDestroy_Tensor,
    .View               = dJacobiView_Tensor,
    .PropogateDown      = dJacobiPropogateDown_Tensor,
    .GetRule            = dJacobiGetRule_Tensor,
    .GetEFS             = dJacobiGetEFS_Tensor,
    .GetNodeCount       = dJacobiGetNodeCount_Tensor,
    .GetConstraintCount = dJacobiGetConstraintCount_Tensor,
    .AddConstraints     = dJacobiAddConstraints_Tensor
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
