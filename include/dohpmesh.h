#ifndef __DOHPMESH_H
#define __DOHPMESH_H

#include "petsc.h"
#include <stdlib.h>
#include <string.h>
#include "dohptype.h"
#include "dohpjacobi.h"

PETSC_EXTERN_CXX_BEGIN

extern PetscCookie dMESH_COOKIE;

EXTERN const char *const iBase_ErrorString[];
EXTERN const char *const iMesh_TopologyName[];
EXTERN const char *const iBase_TypeName[];
EXTERN const char *const iBase_TagValueTypeName[];
EXTERN const int iMesh_TypeFromTopology[];

// #define ICHKERRQ(n) if (n) { dERROR(1,"ITAPS error: %s", iBase_ErrorString[n]); }

/* Unfortunately the explicit `mesh' is necessary to get a useful error string */
#define dICHK(m,e) ICHKERRQ((m),(e))
#define ICHKERRQ(mesh,err)                                     \
  if (err) {                                                   \
    dErr _l_ret = err;                               \
    char           _l_desc[512];                                \
    iMesh_getDescription(mesh,_l_desc,&err,512);dCHK(err); \
    dERROR(1,"%s: %s",iBase_ErrorString[_l_ret],_l_desc);     \
  }

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
} dMeshOrient;
#else
typedef unsigned char dMeshOrient;
#endif

typedef struct {
  char *v;
  dMeshInt a,s;
} MeshListData;

typedef struct {
  dMeshReal *v;
  dMeshInt a, s;
} MeshListReal;

typedef struct {
  dMeshInt *v, a, s;
} MeshListInt;

typedef struct {
  dMeshEH *v;
  dMeshInt a, s;
} MeshListEH;

typedef struct {
  dMeshESH *v;
  dMeshInt a,s;
} MeshListESH;

typedef struct {
  dMeshTag *v;
  dMeshInt a,s;
} MeshListTag;

/* Use this macro to zero a MeshListXXX, i.e. MeshListInt a=MLZ; */
#define MLZ {0,0,0}
#define MeshListFree(m) ((m).a ? (free((m).v),(m).v=0,(m).a=0,(m).s=0,0) : 0)
#define MeshListMalloc(m,n) ( m ? ((n).s=0,(n).a=m,(n).v=malloc((n).a*sizeof(n.v[0])),!(n.v)) : MeshListFree(n))

typedef struct {
  dInt start, stride, end;
} DohpLoopBounds;

/* These tags are used to give full adjacencies because the builtin adjacencies are broken for nonconforming meshes */
#define dTAG_ADJ_REGION_FACE    "dohp_adj_region_face"
#define dTAG_ADJ_FACE_EDGE      "dohp_adj_face_edge"
#define dTAG_ORIENT_REGION_FACE "dohp_orient_region_face"
#define dTAG_ORIENT_FACE_EDGE   "dohp_orient_face_edge"

#define dENT_SET_NAME           "dohp_ent_set_name"
#define dBDY_ROOT               "dohp_bdy_root"
#define dTAG_BDY_NUM            "dohp_bdy_num"
#define dTAG_BDY_NORMAL         "dohp_bdy_normal"

/**
* Assign a canonical orientation of the faces on a hex.  The ordering of faces
* and vertices of the hex is defined by iMesh.  The ordering of edges and
* vertices of the faces are also defined.  We need to find out which rotation of
* the face (in 3D) is required so that the face is in standard orientation.
* Equivalently, we need to know how to compute loop bounds on the face dofs so
* that they are traversed in the forward order when the hex face dofs are
* traversed in forward order.  We will need the permutation P such that for face 'i'
*   region_vertex[dMeshConnectHexQuad[i][j]] = face_vertex[P[j]]
* */
static const int dMeshConnectHexQuad[6][4] = {{0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7}, {0,3,2,1}, {4,5,6,7}};

/**
* Assign a canonical orientation of the edges of a quad.  This is a much simpler analogue of DohpHexQuad.
* We will need the permutation P such that for edge \a i
*   face_vertex[dMeshConnectQuadLine[i][j]] = edge_vertex[P[j]]
* */
static const int dMeshConnectQuadLine[4][2] = {{0,1},{1,2},{2,3},{3,0}};

extern dErr dMeshListIntView(MeshListInt*,const char*);
extern dErr dMeshListEHView(MeshListEH*,const char*);
extern dErr dMeshOrientFindPerm_HexQuad(const iBase_EntityHandle[],const iBase_EntityHandle[],dInt,dMeshOrient*);
extern dErr dMeshOrientFindPerm_QuadLine(const iBase_EntityHandle[],const iBase_EntityHandle[],dInt,dMeshOrient*);

extern dErr dMeshOrientLoopBounds_Quad(dMeshOrient orient, const dInt *size, DohpLoopBounds *l);
extern dErr dMeshOrientLoopBounds_Line(dMeshOrient orient, const dInt *size, DohpLoopBounds *l);

typedef struct p_dMesh *dMesh;

EXTERN dErr dMeshGetLocalNodeNumbering(dMesh,dInt,dInt*,dInt*);
EXTERN dErr dMeshGetTagName(dMesh m,dMeshTag tag,char **name);
EXTERN dErr dMeshLoad(dMesh m,const char fname[],const char opt[]);
EXTERN dErr dMeshCreate(MPI_Comm comm,dMesh *inm);
EXTERN dErr dMeshOrientFacets(dMesh m);
EXTERN dErr dMeshDestroy(dMesh);
EXTERN dErr dMeshView(dMesh,PetscViewer);
EXTERN dErr dMeshRegisterAll(const char path[]);
EXTERN dErr dMeshGetEntSetName(dMesh m,dMeshESH set,char **str);
EXTERN dErr dMeshCreateRuleTagIsotropic(dMesh,dMeshESH,dJacobi,const char*,dInt,dMeshTag*);
EXTERN dErr dMeshDestroyRuleTag(dMesh,dMeshTag);
EXTERN dErr dMeshGetInstance(dMesh,iMesh_Instance*);

PETSC_EXTERN_CXX_END
#endif
