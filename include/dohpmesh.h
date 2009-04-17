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
    dERROR(1,"iMesh(%d) %s: %s",_l_ret,iBase_ErrorString[_l_ret],_l_desc); \
  }
#define dIGCHK(geom,err)                                        \
  if (err) {                                                    \
    dErr _l_ret = err;                                          \
    char _l_desc[512] = "Description not available";            \
    iGeom_getDescription(geom,_l_desc,&err,sizeof(_l_desc));    \
    dERROR(1,"iGeom(%d) %s: %s",_l_ret,iBase_ErrorString[_l_ret],_l_desc); \
  }
#define dIRCHK(assoc,err)                                               \
  if (err) {                                                            \
    dERROR(1,"iRel(%d) %s: %s",err,iBase_ErrorString[iRel_LAST_ERROR.error_type],iRel_LAST_ERROR.description); \
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

/* name tag on manifold sets, NEUMANN_SET is the default when coming from Cubit */
#define dTAG_MANIFOLD_ID      "NEUMANN_SET"
#define dTAG_SENSE            "SENSE"
#define dTAG_EMPTYSET         "EMPTYSET"

#define dMeshType char *
#define dMESHPACK   "pack"
#define dMESHSERIAL "serial"

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
EXTERN dErr dMeshGetRoot(dMesh mesh,dMeshESH *inroot);
EXTERN dErr dMeshCreate(MPI_Comm comm,dMesh *inm);
EXTERN dErr dMeshDestroy(dMesh);
EXTERN dErr dMeshView(dMesh,PetscViewer);
EXTERN dErr dMeshRegisterAll(const char path[]);
#define dMeshRegisterDynamic(a,b,c,d) dMeshRegister(a,b,c,d)
EXTERN dErr dMeshRegister(const char[],const char[],const char[],dErr(*)(dMesh));
EXTERN dErr dMeshSetType(dMesh,const dMeshType);
EXTERN dErr dMeshInitializePackage(const char[]);
EXTERN dErr dMeshCreateRuleTagIsotropic(dMesh,dMeshESH,dJacobi,const char*,dInt,dMeshTag*);
EXTERN dErr dMeshDestroyRuleTag(dMesh,dMeshTag);
EXTERN dErr dMeshGetInstance(dMesh,iMesh_Instance*);
EXTERN dErr dMeshGetNumEnts(dMesh,dMeshESH,dEntType,dEntTopology,dInt*);
EXTERN dErr dMeshGetEnts(dMesh,dMeshESH,dEntType,dEntTopology,dMeshEH[],dInt,dInt*);
EXTERN dErr dMeshGetNumSubsets(dMesh,dMeshESH,dInt,dInt*);
EXTERN dErr dMeshGetSubsets(dMesh,dMeshESH,dInt,dMeshESH[],dInt,dInt*);
EXTERN dErr dMeshGetEntsOff(dMesh,dMeshESH,dInt*,dMeshEH**);
EXTERN dErr dMeshGetAdjIndex(dMesh,const dMeshEH[],dInt,const dMeshEH[],dInt,dInt[],dInt*);

EXTERN dErr dMeshGetTag(dMesh mesh,const char name[],dMeshTag *intag);
EXTERN dErr dMeshTagDestroy(dMesh mesh,dMeshTag tag);
EXTERN dErr dMeshTagCreate(dMesh mesh,const char[],dInt count,dDataType type,dMeshTag *intag);
EXTERN dErr dMeshTagCreateTemp(dMesh mesh,const char[],dInt count,dDataType type,dMeshTag *intag);
EXTERN dErr dMeshTagSetData(dMesh mesh,dMeshTag tag,const dMeshEH ents[],dInt ecount,const void *data,dInt count,dDataType type);
EXTERN dErr dMeshTagGetData(dMesh mesh,dMeshTag tag,const dMeshEH ents[],dInt ecount,void *data,dInt count,dDataType type);
EXTERN dErr dMeshTagSSetData(dMesh mesh,dMeshTag tag,const dMeshESH esets[],dInt ecount,const void *data,dInt count,dDataType type);
EXTERN dErr dMeshTagSGetData(dMesh mesh,dMeshTag tag,const dMeshESH esets[],dInt ecount,void *data,dInt count,dDataType type);
EXTERN dErr dMeshGetTaggedSet(dMesh,dMeshTag,const void*,dMeshESH*);
EXTERN dErr dMeshSetFromOptions(dMesh);
EXTERN dErr dMeshTagBcast(dMesh mesh,dMeshTag tag);
EXTERN dErr dMeshSetCreate(dMesh,dMeshESH*);
EXTERN dErr dMeshGetStatus(dMesh,const dMeshEH[],dInt,dEntStatus[]);
EXTERN dErr dMeshGetTopo(dMesh,dInt,const dMeshEH[],dEntTopology[]);
EXTERN dErr dMeshGetAdjacency(dMesh,dMeshESH,dMeshAdjacency*);
EXTERN dErr dMeshRestoreAdjacency(dMesh,dMeshESH,dMeshAdjacency*);
EXTERN dErr dMeshGetVertexCoords(dMesh,dInt,const dMeshEH[],dInt**,dReal(**)[3]);
EXTERN dErr dMeshRestoreVertexCoords(dMesh,dInt,const dMeshEH[],dInt**,dReal(**)[3]);
EXTERN dErr dMeshPartitionOnOwnership(dMesh,dMeshEH[],dInt,dInt*);

PETSC_EXTERN_CXX_END
#endif
