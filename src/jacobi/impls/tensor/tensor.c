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
#include "polylib.h"
#include "optimalscale.h"
#include <dohpgeom.h>
#include "inlinetmulthex.h"     /* Unrolled variants of TensorMult_Hex_ */

static dErr dJacobiDestroy_Tensor(dJacobi);
static dErr dJacobiView_Tensor(dJacobi,PetscViewer);

static dErr TensorBasisCreate(Tensor tnsr,dInt rsize,const dReal rcoord[],dInt P,TensorBasis *basis)
{
  TensorBasis b;
  const dInt Q=rsize;
  dErr err;

  dFunctionBegin;
  if (!(0 < P && P <= Q)) dERROR(1,"Requested TensorBasis size %d out of bounds.",P);
  if (tnsr->family != GAUSS_LOBATTO) dERROR(1,"GaussFamily %s not supported",GaussFamilies[tnsr->family]);
  err = dNew(struct s_TensorBasis,&b);dCHK(err);
  err = dMallocA6(P*Q,&b->interp,P*Q,&b->deriv,P*Q,&b->interpTranspose,P*Q,&b->derivTranspose,P,&b->node,P,&b->weight);dCHK(err);
  b->Q = Q;
  b->P = P;

  if (P == 1) {         /* degenerate case */
    b->interp[0] = 1.0;
    b->deriv[0] = 0.0;
    b->node[0] = 0.0;
    goto matrices_computed;
  } else {
    const dReal alpha=tnsr->alpha,beta=tnsr->beta;
    dReal *cDeriv;        /* collocation derivative at Gauss-Lobatto points */
    dReal *cDerivT;       /* useless matrix spewed out of Dglj() */
    dReal *interp = b->interp;
    dReal *node = b->node;
    dReal *weight = b->weight;
    err = dMallocA2(P*Q,&cDeriv,P*Q,&cDerivT);dCHK(err);
    zwglj(node,weight,P,alpha,beta);               /* Gauss-Lobatto nodes */
    Imglj(interp,node,rcoord,P,Q,alpha,beta);      /* interpolation matrix */
    Dglj(cDeriv,cDerivT,node,P,alpha,beta);        /* collocation derivative matrix */
    for (dInt i=0; i<Q; i++) {
      for (dInt j=0; j<P; j++) {
        dReal z = 0;
        for (dInt k=0; k<P; k++) {
          z += interp[i*P+k] * cDeriv[k*P+j];
        }
        b->deriv[i*P+j] = z;
      }
    }
    err = dFree2(cDeriv,cDerivT);dCHK(err);
  }
  /* Storing the transposed version explicitly is sort of lame because it costs the same to multiply dense
  * transposed matrices with vectors, however it makes things simple for now. */
  for (dInt i=0; i<Q; i++) {
    for (dInt j=0; j<P; j++) {
      b->interpTranspose[j*Q+i] = b->interp[i*P+j];
      b->derivTranspose[j*Q+i] = b->deriv[i*P+j];
    }
  }

  matrices_computed:
  switch (P) {
#define _C(p) case p: b->mscale = optimal_mscale_ ## p; b->lscale = optimal_lscale_ ## p; break
    _C(2);
    _C(3);
    _C(4);
    _C(5);
    _C(6);
    _C(7);
    _C(8);
    _C(9);
    _C(10);
    _C(11);
    _C(12);
    _C(13);
    _C(14);
#undef _C
    default:
      b->mscale = optimal_ones;
      b->lscale = optimal_ones;
      dERROR(1,"optimal scaling not available for this order, this should just be a PetscInfo warning");
  }
  {
    b->multhex[0] = &TensorMult_Hex_nounroll;
    b->multhex[1] = &TensorMult_Hex_nounroll;
    b->multhex[2] = &TensorMult_Hex_nounroll;
#if defined __SSE3__
    if (P == 4 && Q == 4) {
      b->multhex[0] = &TensorMult_Hex_P4_Q4_D1;
    }
#endif
  }
  *basis = b;
  dFunctionReturn(0);
}

static dErr TensorBasisDestroy(TensorBasis basis)
{
  dErr err;

  dFunctionBegin;
  if (!basis) dFunctionReturn(0);
  err = dFree6(basis->interp,basis->deriv,basis->interpTranspose,basis->derivTranspose,basis->node,basis->weight);dCHK(err);
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

static dErr dJacobiSetFromOptions_Tensor(dJacobi jac)
{
  Tensor this = jac->data;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsHead("Tensor options");dCHK(err);
  {
    err = PetscOptionsTruth("-djac_tensor_mass_scale","Use optimal scaling for Q1 mass matrices","Jacobi",this->usemscale,&this->usemscale,NULL);dCHK(err);
    err = PetscOptionsTruth("-djac_tensor_laplace_scale","Use optimal scaling for Q1 Laplacian matrices","Jacobi",this->uselscale,&this->uselscale,NULL);dCHK(err);
    err = PetscOptionsTruth("-djac_tensor_no_unroll","Do not use unrolled versions of tensor operations","Jacobi",this->nounroll,&this->nounroll,NULL);dCHK(err);
  }
  err = PetscOptionsTail();dCHK(err);
  dFunctionReturn(0);
}

static dErr dJacobiDestroy_Tensor(dJacobi jac)
{
  Tensor tnsr = (Tensor)jac->data;
  khash_t(basis) *bases = tnsr->bases;
  dErr err;

  dFunctionBegin;
  for (khiter_t k=kh_begin(bases); k!=kh_end(bases); k++) {
    if (!kh_exist(bases,k)) continue;
    err = TensorBasisDestroy(kh_val(bases,k));dCHK(err);
  }
  if (jac->quad) {err = dQuadratureDestroy(jac->quad);dCHK(err);}
  kh_destroy_basis(bases);
  err = dJacobiEFSOpsDestroy_Tensor(jac);dCHK(err);
  err = dFree(tnsr);dCHK(err);
  dFunctionReturn(0);
}

static dErr dJacobiView_Tensor(dJacobi jac,dViewer viewer)
{
  Tensor tnsr = (Tensor)jac->data;
  dTruth ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (!ascii) SETERRQ(PETSC_ERR_SUP,"only ASCII");
  err = PetscViewerASCIIPrintf(viewer,"Tensor based Jacobi\n");dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"TensorBasis database.\n");dCHK(err);
  for (khiter_t k=kh_begin(tnsr->bases); k!=kh_end(tnsr->bases); k++) {
    err = TensorBasisView(kh_val(tnsr->bases,k),viewer);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr dJacobiPropogateDown_Tensor(dUNUSED dJacobi jac,dMeshAdjacency a,dInt deg[])
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
      case dTOPO_HEX:
        for (i=0; i<6; i++) {
          ai = a->adjoff[e]+i; aind = a->adjind[ai]; match = a->adjperm[ai];
          deg[aind*3+0] = dMinInt(deg[aind*3+0],deg[e*3+hexperm[i][orient[match][0]]]);
          deg[aind*3+1] = dMinInt(deg[aind*3+1],deg[e*3+hexperm[i][orient[match][1]]]);
          deg[aind*3+2] = 1;
          if (a->topo[aind] != dTOPO_QUAD) dERROR(1,"corrupt adjacency");
        }
        break;
      default: dERROR(1,"Region topology %d not supported",a->topo[e]);dCHK(err);
    }
  }
  for (e=a->toff[dTYPE_FACE]; e<a->toff[dTYPE_FACE+1]; e++) {
    switch (a->topo[e]) {
      case dTOPO_QUAD:
        for (i=0; i<4; i++) {
          ai = a->adjoff[e]+i; aind = a->adjind[ai];
          deg[aind*3+0] = dMinInt(deg[aind*3+0],deg[e*3+quadperm[i]]);
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
      ideg[j] = dMaxInt(0,deg[3*i+j]-2);
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

static dErr dJacobiGetConstraintCount_Tensor(dUNUSED dJacobi jac,dInt nx,const dInt xi[],const dUNUSED dInt xs[],const dInt dUNUSED is[],
                                             const dInt dUNUSED deg[],dMeshAdjacency ma,dInt nnz[],dInt pnnz[])
{

  dFunctionBegin;
  for (dInt i=0; i<nx; i++) {
    const dInt ei = xi[i];
    switch (ma->topo[ei]) {
      case dTOPO_HEX:
        /* \todo Check whether adjacent entities are explicit or Dirichlet when determining the number of entries */
        /* \todo Handle nonconforming elements */
        for (dInt j=xs[i]; j<xs[i+1]; j++) {
          nnz[j] = pnnz[j] = 1;
        }
        break;
      default: dERROR(1,"not implemented for expanded topology %d",ma->topo[i]);
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
    case 5: *ind = (M-1-i)*N + j; break;
    case 6: *ind = (M-1-i) + (N-1-j)*M; break;
    case 7: *ind = i*N + (N-1-j); break;
    default: dERROR(1,"Invalid permutation");
  }
  dFunctionReturn(0);
}

/** Insert the value in both the element assembly matrix and preconditioning assembly matrix */
static dErr PrivateMatSetValue(Mat E,Mat Ep,dInt row,dInt col,MatScalar v,InsertMode imode)
{
  dErr err;

  dFunctionBegin;
  err = MatSetValue(E,row,col,v,imode);dCHK(err);
  if (E != Ep) {
    err = MatSetValue(Ep,row,col,v,imode);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr dJacobiAddConstraints_Tensor(dJacobi dUNUSED jac,dInt nx,const dInt xi[],const dInt xs[],const dInt is[],const dInt deg[],dMeshAdjacency ma,Mat matE,Mat matEp)
{
  dInt elem,i,j,k,l,e[12],eP[12],v[8];
  dInt nrow,ncol,irow[10],icol[30];
  dScalar interp[30]; /* \todo implement nonconforming elements so that this array is actually used (instead of just the first entry)  */
  const dInt *f,*fP;
  const dInt *aI=ma->adjind,*aO=ma->adjoff,*aP=ma->adjperm;
  dErr err;

  dFunctionBegin;
  for (elem=0; elem<nx; elem++) {
    const dInt ei = xi[elem]; /* Element index, \a is, \a deg and everything in \a ma is addressed by \a ei. */
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
          err = PrivateMatSetValue(matE,matEp,xs[elem]+(T[i][0]*d1+T[i][1])*d2+T[i][2],is[v[i]],1,INSERT_VALUES);dCHK(err);
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
            err = PrivateMatSetValue(matE,matEp,irow[0],icol[0],interp[0],INSERT_VALUES);dCHK(err);
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
              dInt facejk[2],faceIndex = -1;
              nrow = ncol = 0;
              irow[nrow++] = (xs[elem] + (start[0]*d1+start[1])*d2+start[2]
                              + (j-start[incd[0]])*scan[incd[0]] + (k-start[incd[1]])*scan[incd[1]]);
              facejk[0] = (j-start[incd[0]])/inc[0];
              facejk[1] = (k-start[incd[1]])/inc[1];
              err = dGeomPermQuadIndex(fP[i],faceDim,facejk,&faceIndex);dCHK(err);
              icol[ncol++] = is[f[i]] + faceIndex;
              interp[0] = 1;
              err = PrivateMatSetValue(matE,matEp,irow[0],icol[0],interp[0],INSERT_VALUES);dCHK(err);
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
              err = PrivateMatSetValue(matE,matEp,irow[0],icol[0],interp[0],INSERT_VALUES);dCHK(err);
            }
          }
        }
        break;
      default: dERROR(1,"not implemented for expanded topology %d",ma->topo[ei]);
    }
  }
  dFunctionReturn(0);
}

static dErr TensorGetBasis(Tensor tnsr,dInt rsize,const dReal rcoord[],dInt bsize,TensorBasis *out)
{
  dErr err;
  int key,new;
  khiter_t k;

  dFunctionBegin;
  if (rsize <= 0) dERROR(1,"Rule size %d must be positive",rsize);
  if (bsize <= 0) dERROR(1,"Basis size %d must be positive",bsize);
  *out = 0;
  key = ((int)tnsr->family << 24) | (rsize << 16) | bsize;
  k = kh_put_basis(tnsr->bases,key,&new);
  if (new) {
    TensorBasis b;
    err = TensorBasisCreate(tnsr,rsize,rcoord,bsize,&b);dCHK(err);
    if (!tnsr->usemscale) b->mscale = optimal_ones;
    if (!tnsr->uselscale) b->lscale = optimal_ones;
    if (tnsr->nounroll) {
      b->multhex[0] = TensorMult_Hex_nounroll;
      b->multhex[1] = TensorMult_Hex_nounroll;
      b->multhex[2] = TensorMult_Hex_nounroll;
    }
    kh_val(tnsr->bases,k) = b;
  }
  *out = kh_val(tnsr->bases,k);
  dFunctionReturn(0);
}

static dErr dJacobiGetEFS_Tensor(dJacobi jac,dInt n,const dEntTopology topo[],const dInt bsize[],dRule rules,dEFS firstefs)
{
  Tensor       this = jac->data;
  dEFS_Tensor  *efs  = (dEFS_Tensor*)firstefs;
  dInt rdim,rsize[3];
  const dReal *rcoord[3];
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<n; i++) {
    err = dRuleGetTensorNodeWeight(&rules[i],&rdim,rsize,rcoord,NULL);dCHK(err);
    switch (topo[i]) {
      case dTOPO_LINE:
        if (rdim != 1) dERROR(1,"Incompatible Rule size %d, expected 1",rdim);
        efs[i].ops = this->efsOpsLine;
        efs[i].rule = &rules[i];
        err = TensorGetBasis(this,rsize[0],rcoord[0],bsize[3*i+0],&efs[i].basis[0]);dCHK(err);
        if (bsize[3*i+1] != 1 || bsize[3*i+2] != 1)
          dERROR(1,"Invalid basis size for Line");
        efs[i].basis[1] = efs[i].basis[2] = NULL;
        break;
      case dTOPO_QUAD:
        if (rdim != 2) dERROR(1,"Incompatible Rule size %d, expected 2",rdim);
        efs[i].ops = this->efsOpsQuad;
        efs[i].rule = &rules[i];
        err = TensorGetBasis(this,rsize[0],rcoord[0],bsize[3*i+0],&efs[i].basis[0]);dCHK(err);
        err = TensorGetBasis(this,rsize[1],rcoord[1],bsize[3*i+1],&efs[i].basis[1]);dCHK(err);
        if (bsize[3*i+2] != 1) dERROR(1,"Invalid basis size for Quad");
        efs[i].basis[2] = NULL;
        break;
      case dTOPO_HEX:
        if (rdim != 3) dERROR(1,"Incompatible Rule size %d, expected 3",rdim);
        efs[i].ops = this->efsOpsHex;
        efs[i].rule = &rules[i];
        err = TensorGetBasis(this,rsize[0],rcoord[0],bsize[3*i+0],&efs[i].basis[0]);dCHK(err);
        err = TensorGetBasis(this,rsize[1],rcoord[1],bsize[3*i+1],&efs[i].basis[1]);dCHK(err);
        err = TensorGetBasis(this,rsize[2],rcoord[2],bsize[3*i+2],&efs[i].basis[2]);dCHK(err);
        break;
      default:
        dERROR(1,"no basis available for given topology");
    }
  }
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
    .SetFromOptions     = dJacobiSetFromOptions_Tensor,
    .Destroy            = dJacobiDestroy_Tensor,
    .View               = dJacobiView_Tensor,
    .PropogateDown      = dJacobiPropogateDown_Tensor,
    .GetEFS             = dJacobiGetEFS_Tensor,
    .GetNodeCount       = dJacobiGetNodeCount_Tensor,
    .GetConstraintCount = dJacobiGetConstraintCount_Tensor,
    .AddConstraints     = dJacobiAddConstraints_Tensor
  };
  Tensor tensor;
  dErr err;

  dFunctionBegin;
  err = dMemcpy(jac->ops,&myops,sizeof(struct _dJacobiOps));dCHK(err);
  err = dNew(struct s_Tensor,&tensor);dCHK(err);

  tensor->usemscale = dFALSE;
  tensor->uselscale = dFALSE;

  tensor->alpha  = 0.0;
  tensor->beta   = 0.0;
  tensor->family = GAUSS_LOBATTO;
  tensor->bases  = kh_init_basis();

  jac->data = tensor;

  err = dJacobiEFSOpsSetUp_Tensor(jac);dCHK(err);
  dFunctionReturn(0);
}
