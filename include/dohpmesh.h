#ifndef _DOHPMESH_H
#define _DOHPMESH_H

#include "petsc.h"
#include <stdlib.h>
#include <string.h>
#include "dohptype.h"
#include "dohpjacobi.h"

PETSC_EXTERN_CXX_BEGIN

typedef struct _p_dMeshPacker *dMeshPacker;
typedef struct _p_dMesh *dMesh;

extern PetscCookie dMESH_COOKIE;

EXTERN const char *const iBase_ErrorString[];
EXTERN const char *const iMesh_TopologyName[];
EXTERN const char *const iBase_TypeName[];
EXTERN const char *const iBase_TagValueTypeName[];
EXTERN const int iMesh_TypeFromTopology[];

/* Unfortunately the explicit `mesh' is necessary to get a useful error string */
#define dICHK(mesh,err)                                         \
  if (err) {                                                    \
    dErr _l_ret = err;                                          \
    char _l_desc[512] = "Description not available";            \
    iMesh_getDescription(mesh,_l_desc,&err,sizeof(_l_desc));    \
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

#define dMeshType char *
#define dMESHPACK   "pack"
#define dMESHSERIAL "serial"

typedef unsigned char dEntStatus;
#define dSTATUS_UNOWNED   (dEntStatus)0x1
#define dSTATUS_SHARED    (dEntStatus)0x2
#define dSTATUS_INTERFACE (dEntStatus)0x4
#define dSTATUS_GHOST     (dEntStatus)0x8

extern dErr dMeshListIntView(MeshListInt*,const char*);
extern dErr dMeshListEHView(MeshListEH*,const char*);

extern dErr dMeshOrientLoopBounds_Quad(dInt orient, const dInt *size, DohpLoopBounds *l);
extern dErr dMeshOrientLoopBounds_Line(dInt orient, const dInt *size, DohpLoopBounds *l);
EXTERN dErr dMeshLoopBounds_Quad(const dInt *size, dInt edge, DohpLoopBounds *l);
EXTERN dErr dMeshLoopBounds_Hex(const dInt *size, dInt face, DohpLoopBounds *l);

EXTERN dErr dMeshGetLocalNodeNumbering(dMesh,dInt,dInt*,dInt*);
EXTERN dErr dMeshGetTagName(dMesh m,dMeshTag tag,char **name);
EXTERN dErr dMeshLoad(dMesh m);
EXTERN dErr dMeshSetInFile(dMesh,const char fname[],const char opt[]);
EXTERN dErr dMeshCreate(MPI_Comm comm,dMesh *inm);
EXTERN dErr dMeshOrientFacets(dMesh m);
EXTERN dErr dMeshDestroy(dMesh);
EXTERN dErr dMeshView(dMesh,PetscViewer);
EXTERN dErr dMeshRegisterAll(const char path[]);
#define dMeshRegisterDynamic(a,b,c,d) dMeshRegister(a,b,c,d)
EXTERN dErr dMeshRegister(const char[],const char[],const char[],dErr(*)(dMesh));
EXTERN dErr dMeshSetType(dMesh,const dMeshType);
EXTERN dErr dMeshInitializePackage(const char[]);
EXTERN dErr dMeshGetEntSetName(dMesh m,dMeshESH set,char **str);
EXTERN dErr dMeshCreateRuleTagIsotropic(dMesh,dMeshESH,dJacobi,const char*,dInt,dMeshTag*);
EXTERN dErr dMeshDestroyRuleTag(dMesh,dMeshTag);
EXTERN dErr dMeshGetInstance(dMesh,iMesh_Instance*);
EXTERN dErr dMeshGetNumEnts(dMesh,dMeshESH,dEntType,dEntTopology,dInt*);
EXTERN dErr dMeshGetEnts(dMesh,dMeshESH,dEntType,dEntTopology,dMeshEH[],dInt,dInt*);
EXTERN dErr dMeshGetAdjIndex(dMesh,const dMeshEH[],dInt,const dMeshEH[],dInt,dInt[],dInt*);

EXTERN dErr dMeshGetTag(dMesh mesh,const char name[],dMeshTag *intag);
EXTERN dErr dMeshTagDestroy(dMesh mesh,dMeshTag tag);
EXTERN dErr dMeshTagCreate(dMesh mesh,const char[],dInt count,dDataType type,dMeshTag *intag);
EXTERN dErr dMeshTagCreateTemp(dMesh mesh,const char[],dInt count,dDataType type,dMeshTag *intag);
EXTERN dErr dMeshTagSetData(dMesh mesh,dMeshTag tag,const dMeshEH ents[],dInt ecount,const void *data,dInt count,dDataType type);
EXTERN dErr dMeshTagGetData(dMesh mesh,dMeshTag tag,const dMeshEH ents[],dInt ecount,void *data,dInt count,dDataType type);
EXTERN dErr dMeshSetFromOptions(dMesh);
EXTERN dErr dMeshTagBcast(dMesh mesh,dMeshTag tag);
EXTERN dErr dMeshGetStatus(dMesh,dInt,const dMeshEH[],dEntStatus[]);
EXTERN dErr dMeshGetTopo(dMesh,dInt,const dMeshEH[],dEntTopology[]);

PETSC_EXTERN_CXX_END
#endif
