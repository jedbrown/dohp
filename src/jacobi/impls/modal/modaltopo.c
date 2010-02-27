#include "modalimpl.h"
#include <dohp.h>

static dErr dEFSView_Modal_Private(const char *name,dRule rule,ModalBasis basis,PetscViewer viewer)
{
  dTruth ascii;
  dErr   err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&ascii);dCHK(err);
  if (ascii) {
    err = PetscViewerASCIIPrintf(viewer,"dEFS type %s\n",name);dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"based on dRule:\n");dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = dRuleView(rule,viewer);dCHK(err);
    err = PetscViewerASCIIPopTab(viewer);
    err = ModalBasisView(basis,viewer);
  }
  dFunctionReturn(0);
}

static dErr dEFSView_Modal_Line(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Modal_Private("Modal_Line",efs->rule,((dEFS_Modal*)efs)->basis,viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr dEFSView_Modal_Quad(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Modal_Private("Modal_Quad",efs->rule,((dEFS_Modal*)efs)->basis,viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr dEFSView_Modal_Hex(dEFS efs,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = dEFSView_Modal_Private("Modal_Hex",efs->rule,((dEFS_Modal*)efs)->basis,viewer);dCHK(err);
  dFunctionReturn(0);
}

static dErr dEFSGetSizes_Modal_All(dEFS efs,dInt *dim,dInt *inodes,dInt *total)
{
  ModalBasis basis = ((dEFS_Modal*)efs)->basis;

  dFunctionBegin;
  if (dim)    *dim = basis->dim;
  if (inodes) *inodes = basis->P;
  if (total)  *total = basis->P;
  dFunctionReturn(0);
}

static dErr dEFSApply_Modal_All(dEFS efs,const dReal dUNUSED jinv[restrict],dInt D,const dScalar in_flat[],dScalar out_flat[],dApplyMode amode,InsertMode imode)
{
  ModalBasis basis = ((dEFS_Modal*)efs)->basis;
  dInt       P = basis->P,Q = basis->Q,dUNUSED dim=basis->dim,i,j,c;
  const dScalar (*interp)[P] = (const dScalar(*)[P])basis->interp;
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
