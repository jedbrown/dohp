#ifndef __DOHPMESH_H
#define __DOHPMESH_H

#include "petsc.h"
#include "iMesh.h"
#include <stdlib.h>
#include <string.h>

PETSC_EXTERN_CXX_BEGIN

typedef int MeshInt;
typedef double MeshReal;
typedef iBase_EntityHandle MeshEH;

EXTERN const char *iBase_ErrorString[];

#define ICHKERRQ(n) if (n) { SETERRQ1(1,"ITAPS error: %s", iBase_ErrorString[n]); }

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
  char *v;
  MeshInt a,s;
} MeshListData;

typedef struct {
  MeshReal *v;
  MeshInt a, s;
} MeshListReal;

typedef struct {
  MeshInt *v, a, s;
} MeshListInt;

typedef struct {
  MeshEH *v;
  MeshInt a, s;
} MeshListEH;

typedef struct {
  iBase_EntitySetHandle *v;
  MeshInt a,s;
} MeshListESH;

typedef struct {
  PetscInt start, stride, end;
} DohpLoopBounds;

#define MeshListFree(m) { if ((m).a) free((m).v); memset(&(m),0,sizeof(m)); }
#define MeshListMalloc(m,n) { (m).s = 0; (m).a = n; (m).v = malloc((m).a*sizeof(m.v[0])); }

/* These tags are used to give full adjacencies because the builtin adjacencies are broken for nonconforming meshes */
#define DOHP_TAG_ADJ_REGION_FACE    "dohp_adj_region_face"
#define DOHP_TAG_ADJ_FACE_EDGE      "dohp_adj_face_edge"
#define DOHP_TAG_ORIENT_REGION_FACE "dohp_orient_region_face"
#define DOHP_TAG_ORIENT_FACE_EDGE   "dohp_orient_face_edge"

#define DOHP_ENT_SET_NAME           "dohp_ent_set_name"
#define DOHP_BDY_ROOT               "dohp_bdy_root"
#define DOHP_TAG_BDY_NUM            "dohp_bdy_num"

/* Assign a canonical orientation of the faces on a hex.  The ordering of faces
* and vertices of the hex is defined by iMesh.  The ordering of edges and
* vertices of the faces are also defined.  We need to find out which rotation of
* the face (in 3D) is required so that the face is in standard orientation.
* Equivalently, we need to know how to compute loop bounds on the face dofs so
* that they are traversed in the forward order when the hex face dofs are
* traversed in forward order.  We will need the permutation P such that for face 'i'
*   region_vertex[DohpHexQuad[i][j]] = face_vertex[P[j]]
* */
static const int DohpHexQuad[6][4] = {{0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7}, {0,3,2,1}, {4,5,6,7}};

/* Assign a canonical orientation of the edges of a quad.  This is a much simpler analogue of DohpHexQuad.
* We will need the permutation P such that for edge 'i'
*   face_vertex[DohpQuadLine[i][j]] = edge_vertex[P[j]]
* */
static const int DohpQuadLine[4][2] = {{0,1},{1,2},{2,3},{3,4}};

EXTERN PetscErrorCode MeshListIntView(MeshListInt *,const char *);
EXTERN PetscErrorCode DohpOrientFindPerm_HexQuad(const iBase_EntityHandle *,const iBase_EntityHandle *,PetscInt,DohpOrient*);
EXTERN PetscErrorCode DohpOrientFindPerm_QuadLine(const iBase_EntityHandle *,const iBase_EntityHandle *,PetscInt,DohpOrient*);

PETSC_EXTERN_CXX_END
#endif
