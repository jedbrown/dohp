/* Implementation of Gauss quadrature with affine coordinate mapping. */

#include "private/dohpimpl.h"

typedef struct {
  dReal *nodes;
  dReal *weights;
  dInt size;
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
  dReal jac[9];
  dReal jinv[9];
  dReal jdet;
} EMap_Affine3;

static dErr dQuotientSetUp_Gauss(dQuotient quot)
{
  dInt nelems;

  dFunctionBegin;
  PetscValidHeaderSpecific(quot,dQUOTIENT_COOKIE,1);
  nelems = quot->nelems;
  dERROR(1,"This is broken");
  //err = PetscMalloc2(nelems,EQuad_Hex,&quot->quad,nelems,EMap_Affine3,&quot->emap);dCHK(err);
  dFunctionReturn(0);
}

static dErr dQuotientUpdate_Gauss(dQuotient q)
{
  dInt *newdegree,i;
  dErr err;

  dFunctionBegin;
  err = PetscMalloc(3*q->nelems*sizeof(dInt),&newdegree);dCHK(err);
  if (q->setdegreefunc) {
    err = (*q->setdegreefunc)(q,q->setdegreectx,q->nelems,newdegree);dCHK(err);
  } else {
    dERROR(1,"SetDegreeFunc not set");
  }
  for (i=0; i<q->nelems; i++) {
    //err = DohpGeomCompute
  }
  err = PetscFree(newdegree);dCHK(err);
  dFunctionReturn(0);
}

static dErr dQuotientDestroy_Gauss(dQuotient q)
{
  dErr err;

  dFunctionBegin;
  err = PetscFree2(q->quad,q->emap);dCHK(err);
  dFunctionReturn(0);
}

#if 0

#undef __FUNCT__
#define __FUNCT__ "EQuotGetJacobian_Affine3"
dErr EQuotGetJacobian_Affine3(EQuot q, dInt qsize, dReal *jac, dReal *jinv, dReal *jdet)
{
  const EMap_Affine3 emap = (EMap_Affine3 *)q->emap;
  dInt size, i;
  dErr err;

  dFunctionBegin;
  for (i=0; i<qsize; i++) {
    err = PetscMemcpy(&jac[i*9], emap->jac, 9*sizeof(dReal));dCHK(err);
    err = PetscMemcpy(&jinv[i*9], emap->jinv, 9*sizeof(dReal));dCHK(err);
    jdet[i] = emap->jdet;
  }
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EQuadGetSize_Line"
dErr EQuadGetSize_Line(void *q, dInt *size)
{
  EQuad_Line *quad = (EQuad_Line *)q;
  dErr err;

  dFunctionBegin;
  *size = quad->base[0]->size;
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EQuadGetSize_Quad"
dErr EQuadGetSize_Quad(void *q, dInt *size)
{
  EQuad_Quad *quad = (EQuad_Quad *)q;
  dErr err;

  dFunctionBegin;
  *size = quad->base[0]->size * quad->base[1]->size;
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EQuadGetSize_Hex"
dErr EQuadGetSize_Hex(void *q, dInt *size)
{
  EQuad_Hex *quad = (EQuad_Hex *)q;
  dErr err;

  dFunctionBegin;
  *size = quad->base[0]->size * quad->base[1]->size * quad->base[2]->size;
  dFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MQuotSetSize"

#endif

dErr dQuotientCreate_Gauss(dQuotient quot)
{

  dFunctionBegin;
  quot->ops->setup = dQuotientSetUp_Gauss;
  quot->ops->update = dQuotientUpdate_Gauss;
  quot->ops->destroy = dQuotientDestroy_Gauss;
  dFunctionReturn(0);
}
