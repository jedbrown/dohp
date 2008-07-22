#ifndef __DOHPMESH_H
#define __DOHPMESH_H

#include "petsc.h"

PETSC_EXTERN_CXX_BEGIN

#define USE_ORIENT_ENUM 0
#if USE_ORIENT_ENUM
typedef enum {
  ORIENT_ID = 0x0,              /* Identity */
  ORIENT_X_REV = 0x1,           /* Reverse the X axis */
  ORIENT_Y_REV = 0x2,
  ORIENT_XY_SWAP = 0x4,
  ORIENT_X_LOW = 0x8,
  ORIENT_X_HIGH = 0x10,
  ORIENT_Y_LOW = 0x40,
  ORIENT_Y_HIGH = 0x80
} DohpOrient;
#else
typedef unsigned short DohpOrient;
#endif

typedef struct {
  double *v;
  int a, s;
} MeshListDouble;

typedef struct {
  int *v, a, s;
} MeshListInt;

typedef struct {
  iBase_EntityHandle *v;
  int a, s;
} MeshListEH;

#define MeshListFree(m) { if ((m).a) free((m).v); memset(&(m),0,sizeof(m)); }
#define MeshListMalloc(m,n) { (m).a = (m).s = n; (m).v = malloc((m).a*sizeof(m.v[0])); }

const char *DOHP_ORIENT_REGION_FACE = "dohp_orient_region_face";
const char *DOHP_ORIENT_FACE_EDGE = "dohp_orient_face_edge";
const char *DOHP_ADJ_REGION_FACE = "dohp_adj_region_face";
const char *DOHP_ADJ_FACE_EDGE = "dohp_adj_face_edge";

#undef __FUNCT__
#define __FUNCT__ "MeshListIntView"
PetscErrorCode MeshListIntView(MeshListInt *m, const char *name) {
  CHKERRQ(PetscPrintf(PETSC_COMM_SELF,"# %s [%d]\n", name, m->s));
  CHKERRQ(PetscIntView(m->s,m->v,PETSC_VIEWER_STDOUT_SELF));
  PetscFunctionReturn(0);
}

PETSC_EXTERN_CXX_END
#endif
