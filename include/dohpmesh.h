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
#define ICHKERRQ(mesh,err)                                      \
  if (err) {                                                    \
    dErr _l_ret = err;                                          \
    char _l_desc[512];                                          \
    iMesh_getDescription(mesh,_l_desc,&err,512);dCHK(err);      \
    dERROR(1,"%s: %s",iBase_ErrorString[_l_ret],_l_desc);       \
  }

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
#define MLREF(m) &(m).v,&(m).a,&(m).s

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

extern dErr dMeshListIntView(MeshListInt*,const char*);
extern dErr dMeshListEHView(MeshListEH*,const char*);

extern dErr dMeshOrientLoopBounds_Quad(dGeomOrient orient, const dInt *size, DohpLoopBounds *l);
extern dErr dMeshOrientLoopBounds_Line(dGeomOrient orient, const dInt *size, DohpLoopBounds *l);

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
