/* Implementation of Gauss quadrature with affine coordinate mapping. */

#include "private/dohpimpl.h"

typedef struct {
  PetscReal *nodes;
  PetscReal *weights;
  PetscInt size;
} EQuad_Base;

typedef struct {
  EQuad_Base base[1];
} EQuad_Line;

typedef struct {
  EQuad_Base base[2];
} EQuad_Quad;

typedef struct {
  EQuad_Base base[3];
} EQuad_Hex;

typedef struct {
  PetscReal jac[9];
  PetscReal jinv[9];
  PetscReal jdet;
} EMap_Affine3;

#undef __FUNCT__
#define __FUNCT__ "EQuotGetJacobian_Affine3"
PetscErrorCode EQuotGetJacobian_Affine3(EQuot q, PetscInt qsize, PetscReal *jac, PetscReal *jinv, PetscReal *jdet)
{
  const EMap_Affine3 emap = (EMap_Affine3 *)q->emap;
  PetscInt size, i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (i=0; i<qsize; i++) {
    ierr = PetscMemcpy(&jac[i*9], emap->jac, 9*sizeof(PetscReal));CHKERRQ(ierr);
    ierr = PetscMemcpy(&jinv[i*9], emap->jinv, 9*sizeof(PetscReal));CHKERRQ(ierr);
    jdet[i] = emap->jdet;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EQuadGetSize_Line"
PetscErrorCode EQuadGetSize_Line(void *q, PetscInt *size)
{
  EQuad_Line *quad = (EQuad_Line *)q;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *size = quad->base[0]->size;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EQuadGetSize_Quad"
PetscErrorCode EQuadGetSize_Quad(void *q, PetscInt *size)
{
  EQuad_Quad *quad = (EQuad_Quad *)q;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *size = quad->base[0]->size * quad->base[1]->size;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EQuadGetSize_Hex"
PetscErrorCode EQuadGetSize_Hex(void *q, PetscInt *size)
{
  EQuad_Hex *quad = (EQuad_Hex *)q;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *size = quad->base[0]->size * quad->base[1]->size * quad->base[2]->size;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MQuotSetSize"
