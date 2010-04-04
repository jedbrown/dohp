#include "modalimpl.h"
#include <dohp.h>
#include <dohpmesh.h>           /* iMesh_TypeFromTopology */

static dErr ModalBasisView(dInt Q,dInt P,const dReal interp[],const dReal deriv[],PetscViewer viewer)
{
  dTruth ascii;
  dErr   err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  err = PetscViewerASCIIPrintf(viewer,"Modal basis with %d point rule and %d point basis.\n",Q,P);dCHK(err);
  err = dRealTableView(Q,P,interp,"interp",viewer);dCHK(err);
  err = dRealTableView(Q,3*P,deriv,"deriv",viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dEFSView_Modal_Private(const dEFS_Modal *modal,const char *name,PetscViewer viewer)
{
  dTruth ascii;
  dErr   err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (ascii) {
    err = PetscViewerASCIIPrintf(viewer,"dEFS type %s\n",name);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"based on dRule:\n");dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = dRuleView(modal->rule,viewer);dCHK(err);
    err = PetscViewerASCIIPopTab(viewer);
    err = ModalBasisView(modal->Q,modal->P,modal->interp,modal->deriv,viewer);
  }
  dFunctionReturn(0);
}

static dErr dEFSView_Modal_Line(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Modal_Private((dEFS_Modal*)efs,"Modal_Line",viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr dEFSView_Modal_Quad(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Modal_Private((dEFS_Modal*)efs,"Modal_Quad",viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr dEFSView_Modal_Hex(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Modal_Private((dEFS_Modal*)efs,"Modal_Hex",viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dEFSGetSizes_Modal_All(dEFS efs,dInt *dim,dInt *inodes,dInt *total)
{
  const dEFS_Modal *e = (dEFS_Modal*)efs;

  dFunctionBegin;
  if (dim)    *dim    = dMeshEntTypeFromTopology(e->topo);
  if (inodes) *inodes = e->P;
  if (total)  *total  = e->P;
  dFunctionReturn(0);
}

static dErr dEFSApply_Modal_All(dEFS efs_generic,const dReal dUNUSED jinv[restrict],dInt D,const dScalar in_flat[],dScalar out_flat[],dApplyMode amode,InsertMode imode)
{
  dEFS_Modal *efs = (dEFS_Modal*)efs_generic;
  dInt       P = efs->P,Q = efs->Q,i,j,c,d;
  const dReal (*interp)[P] = (const dReal(*)[P])efs->interp;
  const dReal (*deriv)[P][3] = (const dReal(*)[P][3])efs->deriv;
  dErr err;

  dFunctionBegin;
  switch (amode) {
    case dAPPLY_INTERP: {
      const dScalar (*in)[D] = (const dScalar(*)[D])in_flat;
      dScalar (*restrict out)[D] = (dScalar(*)[D])out_flat;
      if (imode == INSERT_VALUES) {err = dMemzero(&out[0][0],Q*D*sizeof(dScalar));dCHK(err);}
      for (i=0; i<Q; i++) {
        for (j=0; j<P; j++) {
          for (c=0; c<D; c++) {
            out[i][c] += interp[i][j] * in[j][c];
          }
        }
      }
    } break;
    case dAPPLY_INTERP_TRANSPOSE: {
      const dScalar (*in)[D] = (const dScalar(*)[D])in_flat;
      dScalar (*restrict out)[D] = (dScalar(*)[D])out_flat;
      if (imode == INSERT_VALUES) {err = dMemzero(&out[0][0],P*D*sizeof(dScalar));dCHK(err);}
      for (i=0; i<Q; i++) {
        for (j=0; j<P; j++) {
          for (c=0; c<D; c++) {
            out[j][c] += interp[i][j] * in[i][c];
          }
        }
      }
    } break;
    case dAPPLY_GRAD: {
      const dScalar (*in)[D] = (const dScalar(*)[D])in_flat;
      dScalar (*restrict out)[D][3] = (dScalar(*)[D][3])out_flat;
      if (imode == INSERT_VALUES) {err = dMemzero(&out[0][0],Q*D*3*sizeof(dScalar));dCHK(err);}
      for (i=0; i<Q; i++) {
        for (j=0; j<P; j++) {
          for (c=0; c<D; c++) {
            for (d=0; d<3; d++) {
              out[i][c][d] += deriv[i][j][d] * in[j][c];
            }
          }
        }
      }
    } break;
    default: dERROR(PETSC_ERR_SUP,"Apply-mode not implemented");
  }
  dFunctionReturn(0);
}

dErr dJacobiEFSOpsSetUp_Modal(dJacobi jac)
{
    static const struct _dEFSOps efsOpsLine = {
    .view     = dEFSView_Modal_Line,
    .getSizes = dEFSGetSizes_Modal_All,
    .apply    = dEFSApply_Modal_All,
  };
  static const struct _dEFSOps efsOpsQuad = {
    .view     = dEFSView_Modal_Quad,
    .getSizes = dEFSGetSizes_Modal_All,
    .apply    = dEFSApply_Modal_All,
  };
  static const struct _dEFSOps efsOpsHex = {
    .view     = dEFSView_Modal_Hex,
    .getSizes = dEFSGetSizes_Modal_All,
    .apply    = dEFSApply_Modal_All,
  };
  dJacobi_Modal *modal = jac->data;
  dErr err;

  dFunctionBegin;
  err = dMallocA3(1,&modal->efsOpsLine,1,&modal->efsOpsQuad,1,&modal->efsOpsHex);dCHK(err);
  err = dMemcpy(modal->efsOpsLine,&efsOpsLine,sizeof(struct _dEFSOps));dCHK(err);
  err = dMemcpy(modal->efsOpsQuad,&efsOpsQuad,sizeof(struct _dEFSOps));dCHK(err);
  err = dMemcpy(modal->efsOpsHex,&efsOpsHex,sizeof(struct _dEFSOps));dCHK(err);
  dFunctionReturn(0);
}
