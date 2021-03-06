#include <dohpmeshimpl.h>
#include <dohpgeom.h>
#include <dohpviewer.h>
#include <iMesh_extensions.h>
#include <MBParallelConventions.h>
#include <errno.h>

const char *const iBase_ErrorString[] = {
  [iBase_SUCCESS] = "iBase_SUCCESS",
  [iBase_MESH_ALREADY_LOADED] = "iBase_MESH_ALREADY_LOADED",
  [iBase_FILE_NOT_FOUND] = "iBase_FILE_NOT_FOUND",
  [iBase_FILE_WRITE_ERROR] = "iBase_FILE_WRITE_ERROR",
  [iBase_NIL_ARRAY] = "iBase_NIL_ARRAY",
  [iBase_BAD_ARRAY_SIZE] = "iBase_BAD_ARRAY_SIZE",
  [iBase_BAD_ARRAY_DIMENSION] = "iBase_BAD_ARRAY_DIMENSION",
  [iBase_INVALID_ENTITY_HANDLE] = "iBase_INVALID_ENTITY_HANDLE",
  [iBase_INVALID_ENTITY_COUNT] = "iBase_INVALID_ENTITY_COUNT",
  [iBase_INVALID_ENTITY_TYPE] = "iBase_INVALID_ENTITY_TYPE",
  [iBase_INVALID_ENTITY_TOPOLOGY] = "iBase_INVALID_ENTITY_TOPOLOGY",
  [iBase_BAD_TYPE_AND_TOPO] = "iBase_BAD_TYPE_AND_TOPO",
  [iBase_ENTITY_CREATION_ERROR] = "iBase_ENTITY_CREATION_ERROR",
  [iBase_INVALID_TAG_HANDLE] = "iBase_INVALID_TAG_HANDLE",
  [iBase_TAG_NOT_FOUND] = "iBase_TAG_NOT_FOUND",
  [iBase_TAG_ALREADY_EXISTS] = "iBase_TAG_ALREADY_EXISTS",
  [iBase_TAG_IN_USE] = "iBase_TAG_IN_USE",
  [iBase_INVALID_ENTITYSET_HANDLE] = "iBase_INVALID_ENTITYSET_HANDLE",
  [iBase_INVALID_ITERATOR_HANDLE] = "iBase_INVALID_ITERATOR_HANDLE",
  [iBase_INVALID_ARGUMENT] = "iBase_INVALID_ARGUMENT",
  [iBase_MEMORY_ALLOCATION_FAILED] = "iBase_MEMORY_ALLOCATION_FAILED",
  [iBase_NOT_SUPPORTED] = "iBase_NOT_SUPPORTED",
  [iBase_FAILURE] = "iBase_FAILURE",
};

static dInt iBase_SizeFromType(dDataType type) {
  static const dInt sizes [dDATA_UB] = {
    [dDATA_BYTE] = sizeof(char),
    [dDATA_INT]  = sizeof(int),
    [dDATA_REAL] = sizeof(dReal),
    [dDATA_EH]   = sizeof(dMeshEH),
    [dDATA_ESH]  = sizeof(dMeshESH),
  };
  return (type < sizeof(sizes)/sizeof(sizes[0]))
    ? sizes[type] : -1;
}

const char *dMeshEntTopologyName(dEntTopology topo) {
  static const char *const iMesh_TopologyName[12] = {
    "iMesh_POINT",
    "iMesh_LINE_SEGMENT",
    "iMesh_POLYGON",
    "iMesh_TRIANGLE",
    "iMesh_QUADRILATERAL",
    "iMesh_POLYHEDRON",
    "iMesh_TETRAHEDRON",
    "iMesh_HEXAHEDRON",
    "iMesh_PRISM",
    "iMesh_PYRAMID",
    "iMesh_SEPTAHEDRON",
    "iMesh_ALL_TOPOLOGIES"
  };
  return (topo <= sizeof(iMesh_TopologyName)/sizeof(iMesh_TopologyName[0]))
    ? iMesh_TopologyName[topo]
    : "iMesh_TOPOLOGY_OUT_OF_RANGE";
}

dEntType dMeshEntTypeFromTopology(dEntTopology topo) {
  static const int iMesh_TypeFromTopology[12] = {
    [dTOPO_POINT]       = dTYPE_VERTEX,
    [dTOPO_LINE]        = dTYPE_EDGE,
    [dTOPO_POLYGON]     = dTYPE_FACE,
    [dTOPO_TRIANGLE]    = dTYPE_FACE,
    [dTOPO_QUAD]        = dTYPE_FACE,
    [dTOPO_POLYHEDRON]  = dTYPE_REGION,
    [dTOPO_TET]         = dTYPE_REGION,
    [dTOPO_HEX]         = dTYPE_REGION,
    [dTOPO_PRISM]       = dTYPE_REGION,
    [dTOPO_PYRAMID]     = dTYPE_REGION,
    [dTOPO_SEPTAHEDRON] = dTYPE_REGION,
    [dTOPO_ALL]         = dTYPE_ALL,
  };
  return (topo <= dTOPO_ALL)
    ? iMesh_TypeFromTopology[topo]
    : dTYPE_ALL;
}

const char *dMeshEntTypeName(dEntType type) {
  static const char *const iBase_TypeName[] = {
    [dTYPE_VERTEX] = "iBase_VERTEX",
    [dTYPE_EDGE]   = "iBase_EDGE",
    [dTYPE_FACE]   = "iBase_FACE",
    [dTYPE_REGION] = "iBase_REGION",
    [dTYPE_ALL]    = "iBase_ALL_TYPES",
  };
  return (type <= sizeof(iBase_TypeName)/sizeof(iBase_TypeName[0]))
    ? iBase_TypeName[type]
    : "iBase_TYPE_OUT_OF_RANGE";
}

const char *const iBase_TagValueTypeName[] = {
  [iBase_BYTES]             = "iBase_BYTES",
  [iBase_INTEGER]           = "iBase_INTEGER",
  [iBase_DOUBLE]            = "iBase_DOUBLE",
  [iBase_ENTITY_HANDLE]     = "iBase_ENTITY_HANDLE",
  [iBase_ENTITY_SET_HANDLE] = "iBase_ENTITY_SET_HANDLE",
};

dErr dMeshListIntView(MeshListInt *ml,const char *name)
{
  dErr err;

  dFunctionBegin;
  err = dPrintf(PETSC_COMM_SELF,"# %s [%d]\n", name, ml->s);dCHK(err);
  err = PetscIntView(ml->s,ml->v,PETSC_VIEWER_STDOUT_SELF);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshListEHView(MeshListEH *ml,const char *name)
{
  dInt n=ml->s/20,p=ml->s%20;
  dErr err;

  dFunctionBegin;
  err = dPrintf(PETSC_COMM_SELF,"# %s [%d]\n",name,ml->s);dCHK(err);
  for (dInt i=0; i<n; i++) {
    err = dPrintf(PETSC_COMM_SELF,"%D:",i*20);dCHK(err);
    for (dInt j=0; j<20; j++) {
      err = dPrintf(PETSC_COMM_SELF," %#4x",0xffffffff & (long)ml->v[i*20+j]);dCHK(err);
    }
    err = dPrintf(PETSC_COMM_SELF,"\n");dCHK(err);
  }
  if (p) {
    err = dPrintf(PETSC_COMM_SELF,"%D:",n*20);dCHK(err);
    for (dInt i=0; i<p; i++) {
      err = dPrintf(PETSC_COMM_SELF," %#4x",0xffffffff & (long)ml->v[n*20+i]);dCHK(err);
    }
    err = dPrintf(PETSC_COMM_SELF,"\n");dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dMeshOrientLoopBounds_Quad(dInt orient, const dInt *size, DohpLoopBounds *l)
{
  const dInt ox=size[0], oy=size[1];

  dFunctionBegin;
  switch (orient) {
    case 0: {
      l[0].start = 0;         l[0].stride = oy;  l[0].end = ox*oy;
      l[1].start = 0;         l[1].stride = 1;   l[1].end = oy;
    } break;
    case 1: {
      l[0].start = 0;         l[0].stride = 1;   l[0].end = oy;
      l[1].start = (ox-1)*oy; l[1].stride = -oy; l[1].end = -oy;
    } break;
    case 2: {
      l[0].start = (ox-1)*oy; l[0].stride = -oy; l[0].end = -oy;
      l[1].start = oy-1;      l[1].stride = -1;  l[1].end = -1;
    } break;
    case 3: {
      l[0].start = oy-1;      l[0].stride = -1;  l[0].end = -1;
      l[1].start = 0;         l[1].stride = oy;  l[1].end = ox*oy;
    } break;
    case 4: {
      l[0].start = 0;         l[0].stride = 1;   l[0].end = oy;
      l[1].start = 0;         l[1].stride = oy;  l[1].end = ox*oy;
    } break;
    case 5: {
      l[0].start = 0;         l[0].stride = oy;  l[0].end = ox*oy;
      l[1].start = oy-1;      l[1].stride = -1;  l[1].end = -1;
    } break;
    case 6: {
      l[0].start = oy-1;      l[0].stride = -1;  l[0].end = -1;
      l[1].start = (ox-1)*oy; l[1].stride = -oy; l[1].end = -oy;
    } break;
    case 7: {
      l[0].start = (ox-1)*oy; l[0].stride = -oy;   l[0].end = -oy;
      l[1].start = 0;         l[1].stride = 1;     l[1].end = oy;
    } break;
    default:
      dERROR(PETSC_COMM_SELF,1,"Orientation not supported.");
  }
  dFunctionReturn(0);
}

dErr dMeshOrientLoopBounds_Line(dInt orient, const dInt *size, DohpLoopBounds *l)
{

  dFunctionBegin;
  switch (orient) {
    case 0: l->start = 0;         l->stride = 1;  l->end = size[0]; break;
    case 1: l->start = size[0]-1; l->stride = -1; l->end = -1;       break;
    default: dERROR(PETSC_COMM_SELF,1,"Orientation not supported.");
  }
  dFunctionReturn(0);
}

/* On each face, we need a loop which traverses the face (indicated by
* DohpHexQuad[][]) in the positive order.  The ordering of degrees of freedom on the
* Hex is [i][j][k] (C-style ordering). */
dErr dMeshLoopBounds_Hex(const dInt *size, dInt face, DohpLoopBounds *l)
{
  const dInt ox=size[0], oy=size[1], oz=size[2];
  dFunctionBegin;
  switch (face) {
    case 0: {                   /* 0,1,5,4 */
      l[0].start = 0;            l[0].stride = oy*oz;  l[0].end = ox*oy*oz;
      l[1].start = 0;            l[1].stride = 1;      l[1].end = oz;
    } break;
    case 1: {                   /* 1,2,6,5 */
      l[0].start = (ox-1)*oy*oz; l[0].stride = oz;     l[0].end = ox*oy*oz;
      l[1].start = 0;            l[1].stride = 1;      l[1].end = oz;
    } break;
    case 2: {                   /* 2,3,7,6 */
      l[0].start = (ox*oy-1)*oz; l[0].stride = -oy*oz; l[0].end = -oz;
      l[1].start = 0;            l[1].stride = 1;      l[1].end = oz;
    } break;
    case 3: {                   /* 3,0,4,7 */
      l[0].start = (oy-1)*oz;    l[0].stride = -oz;    l[0].end = -oz;
      l[1].start = 0;            l[1].stride = 1;      l[1].end = oz;
    } break;
    case 4: {                   /* 0,3,2,1 */
      l[0].start = 0;            l[0].stride = oz;     l[0].end = oy*oz;
      l[1].start = 0;            l[1].stride = oy*oz;  l[1].end = ox*oy*oz;
    } break;
    case 5: {                   /* 4,5,6,7 */
      l[0].start = oz-1;         l[0].stride = oy*oz;  l[0].end = ox*oy*oz;
      l[1].start = 0;            l[1].stride = oz;     l[1].end = oy*oz;
    } break;
    default:
      dERROR(PETSC_COMM_SELF,1,"Face number not recognized.");
  }
  dFunctionReturn(0);
}

dErr dMeshLoopBounds_Quad(const dInt *size, dInt edge, DohpLoopBounds *l)
{
  const dInt ox=size[0], oy=size[1];

  dFunctionBegin;
  switch (edge) {
    case 0:                     /* 0,1 */
      l->start = 0;         l->stride = oy;  l->end = ox*oy; break;
    case 1:                     /* 1,2 */
      l->start = (ox-1)*oy; l->stride = 1;   l->end = ox*oy; break;
    case 2:                     /* 2,3 */
      l->start = ox*oy-1;   l->stride = -oy; l->end = -1; break;
    case 3:                     /* 3,0 */
      l->start = oy-1;      l->stride = -1;  l->end = -1; break;
    default:
      dERROR(PETSC_COMM_SELF,1,"Edge number not recognized.");
  }
  dFunctionReturn(0);
}

// Causes dMeshLoad() to read this file name.
dErr dMeshSetInFile(dMesh mesh,const char *fname,const char *options)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (fname) {
    err = PetscFree(mesh->infile);dCHK(err);
    err = PetscStrallocpy(fname,&mesh->infile);dCHK(err);
  }
  if (options) {
    err = PetscFree(mesh->inoptions);dCHK(err);
    err = PetscStrallocpy(options,&mesh->inoptions);dCHK(err);
  }
  dFunctionReturn(0);
}

// Set the mesh generation type to be used when dMeshLoad() is called.
dErr dMeshSetGenerate(dMesh mesh,dMeshGenType gtype)
{
  dFunctionBegin;
  mesh->gentype = gtype;
  dFunctionReturn(0);
}

dErr dMeshGetRoot(dMesh mesh,dMeshESH *inroot)
{
  iMesh_Instance mi = mesh->mi;
  dIInt ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(inroot,2);
  *inroot = 0;
  if (!mesh->root) {
    iMesh_getRootSet(mi,&mesh->root,&ierr);dICHK(mi,ierr);
  }
  *inroot = mesh->root;
  dFunctionReturn(0);
}

/* Makes a new set containing all entities in original, but not subsets */
dErr dMeshSetDuplicateEntsOnly(dMesh mesh,dMeshESH set,dMeshESH *copy)
{
  dInt    size;
  dMeshEH *ents;
  dMeshSetOrdering ordering;
  dErr    err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(copy,3);
  err = dMeshGetNumEnts(mesh,set,dTYPE_ALL,dTOPO_ALL,&size);dCHK(err);
  err = dMallocA(size,&ents);dCHK(err);
  err = dMeshGetEnts(mesh,set,dTYPE_ALL,dTOPO_ALL,ents,size,NULL);dCHK(err);
  err = dMeshSetGetOrdering(mesh,set,&ordering);dCHK(err);
  err = dMeshSetCreate(mesh,ordering,copy);dCHK(err);
  err = dMeshSetAddEnts(mesh,*copy,ents,size);dCHK(err);
  err = dFree(ents);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshSetGetOrdering(dMesh mesh,dMeshESH set,dMeshSetOrdering *ordering)
{
  dIInt ierr,islist;

  dFunctionBegin;
  iMesh_isList(mesh->mi,set,&islist,&ierr);dICHK(mesh->mi,ierr);
  *ordering = islist ? dMESHSET_ORDERED : dMESHSET_UNORDERED;
  dFunctionReturn(0);
}

/**
* This function allocates memory for the tag name.  It should be freed with PetscFree()
*/
dErr dMeshGetTagName(dMesh mesh,dMeshTag tag,char **name)
{
  iMesh_Instance mi = mesh->mi;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(mesh,dMESH_CLASSID,1);
  dValidPointer(name,2);
  err = PetscMalloc(dNAME_LEN,name);dCHK(err);
  iMesh_getTagName(mi,tag,*name,&err,dNAME_LEN);dICHK(mi,err);
  dFunctionReturn(0);
}

dErr dMeshGetTag(dMesh mesh,const char name[],dMeshTag *intag)
{
  dIInt ierr;
  size_t namelen;
  dMeshTag tag;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidCharPointer(name,2);
  dValidPointer(intag,3);
  *intag = 0;
  err = dStrlen(name,&namelen);dCHK(err);
  iMesh_getTagHandle(mesh->mi,name,&tag,&ierr,(int)namelen);dICHK(mesh->mi,ierr);
  *intag = tag;
  dFunctionReturn(0);
}

dErr dMeshTagCreateTemp(dMesh mesh,const char template[],dInt count,dDataType type,dMeshTag *intag)
{
  static dInt unique_id = 0;
  char name[dNAME_LEN];
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(intag,4);
  err = PetscSNPrintf(name,sizeof(name),"__DOHP_%s_%d",template,unique_id++);dCHK(err);
  err = dMeshTagCreate(mesh,name,count,type,intag);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshTagCreate(dMesh mesh,const char name[],dInt count,dDataType type,dMeshTag *intag)
{
  dMeshTag       tag;
  iMesh_Instance mi;
  dIInt          ierr,itype;
  size_t         namelen;
  dErr           err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(intag,4);
  *intag = 0;
  mi = mesh->mi;
  if (count <= 0) dERROR(PETSC_COMM_SELF,1,"Tags must have positive count");
  err = dDataTypeToITAPS(type,&itype);dCHK(err);
  err = dStrlen(name,&namelen);dCHK(err);
  iMesh_getTagHandle(mi,name,&tag,&ierr,(dIInt)namelen);
  if (ierr == iBase_TAG_NOT_FOUND) {
    iMesh_createTag(mi,name,count,itype,&tag,&ierr,(dIInt)namelen);dICHK(mi,ierr);
    {
      dIInt bytes;
      iMesh_getTagSizeBytes(mi,tag,&bytes,&ierr);dICHK(mi,ierr);
    }
  } else {
    dIInt existing_itype,existing_count;
    dICHK(mi,ierr);
    iMesh_getTagType(mi,tag,&existing_itype,&ierr);dICHK(mi,ierr);
    if (itype != existing_itype) dERROR(PETSC_COMM_SELF,1,"tag exists with different type");
    iMesh_getTagSizeValues(mi,tag,&existing_count,&ierr);dICHK(mi,ierr);
    if (count != existing_count) dERROR(PETSC_COMM_SELF,1,"tag exists with same type but different count");
  }
  *intag = tag;
  dFunctionReturn(0);
}


dErr dMeshTagDestroy(dMesh mesh,dMeshTag tag)
{
  dIInt ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  iMesh_destroyTag(mesh->mi,tag,1,&ierr);dICHK(mesh->mi,ierr);
  dFunctionReturn(0);
}

dErr dMeshTagSetData(dMesh mesh,dMeshTag tag,const dMeshEH ents[],dInt ecount,const void *data,dInt count,dDataType type)
{
  iMesh_Instance mi = mesh->mi;
  const char *dptr = data;
  dIInt size,ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(ents,3);
  dValidPointer(data,5);
  size = count * iBase_SizeFromType(type);
  if (1) {
    dIInt bytes;
    iMesh_getTagSizeBytes(mi,tag,&bytes,&ierr);dICHK(mi,ierr);
    if (size != ecount * (dInt)bytes)
      dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Trying to tag %D entities with %D byte tags, but only requested %D bytes",ecount,bytes,size);
  }
  iMesh_setArrData(mi,ents,ecount,tag,dptr,size,&ierr);dICHK(mi,ierr);
  dFunctionReturn(0);
}

dErr dMeshTagGetData(dMesh mesh,dMeshTag tag,const dMeshEH ents[],dInt ecount,void *data,dInt count,dDataType type)
{
  iMesh_Instance mi = mesh->mi;
  dIInt size,alloc,ierr;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (!ecount) dFunctionReturn(0);
  dValidPointer(ents,3);
  dValidPointer(data,5);
  alloc = count * iBase_SizeFromType(type);
  size = 0;                     /* protect against degenerate case */
  if (1) {
    dIInt bytes;
    char *name;
    err = dMeshGetTagName(mesh,tag,&name);dCHK(err);
    iMesh_getTagSizeBytes(mi,tag,&bytes,&ierr); if (ierr == iBase_TAG_NOT_FOUND) dERROR(PETSC_COMM_SELF,PETSC_ERR_LIB,"Cannot find tag '%s'\n",name); else dICHK(mi,ierr);
    if (ecount && bytes != iBase_SizeFromType(type)*count/ecount) {
      dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Trying ot retrieve tags of %D bytes, but caller expects %D",bytes,iBase_SizeFromType(type)*count/ecount);
    }
    err = dFree(name);dCHK(err);
  }
  iMesh_getArrData(mi,ents,ecount,tag,&data,&alloc,&size,&ierr);dICHK(mi,ierr);
  if (alloc != count * iBase_SizeFromType(type))
    dERROR(PETSC_COMM_SELF,1,"Looks like an iMesh inconsistency, the library shouldn't be messing with this");
  if (size > alloc) dERROR(PETSC_COMM_SELF,1,"Insufficient allocation, iMesh should have thrown an error already");
  dFunctionReturn(0);
}

dErr dMeshTagSSetData(dMesh mesh,dMeshTag tag,const dMeshESH esets[],dInt ecount,const void *data,dInt count,dDataType type)
{
  iMesh_Instance mi = mesh->mi;
  const char *dptr = data;
  dIInt size,ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(esets,3);
  dValidPointer(data,5);
  size = count * iBase_SizeFromType(type);
  if (1) {
    dIInt bytes;
    iMesh_getTagSizeBytes(mi,tag,&bytes,&ierr);dICHK(mi,ierr);
    if (size != ecount * (dInt)bytes)
      dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"Trying to tag %D entities with %D byte tags, but only requested %D bytes",ecount,bytes,size);
  }
#if 0
  iMesh_setArrData(mi,(const iBase_EntityHandle*)esets,ecount,tag,dptr,size,&ierr);dICHK(mi,ierr);
#else
  for (dInt i=0; i<ecount; i++) {
    iMesh_setEntSetData(mi,esets[i],tag,dptr+i*iBase_SizeFromType(type),iBase_SizeFromType(type),&ierr);dICHK(mi,ierr);
  }
#endif
  dFunctionReturn(0);
}

dErr dMeshTagSGetData(dMesh mesh,dMeshTag tag,const dMeshESH esets[],dInt ecount,void *data,dInt count,dDataType type)
{
  iMesh_Instance mi = mesh->mi;
  dIInt size,alloc,peritem,ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (!ecount) dFunctionReturn(0);
  dValidPointer(esets,3);
  dValidPointer(data,5);
  alloc = count * iBase_SizeFromType(type);
  peritem = alloc / ecount;
  size = 0;                     /* protect against degenerate case */
#if 0
  iMesh_getArrData(mi,(const iBase_EntityHandle*)esets,ecount,tag,&dptr,&alloc,&size,&ierr);dICHK(mi,ierr);
  if (dptr != (char*)data || alloc != count * iBase_SizeFromType[type])
    dERROR(PETSC_COMM_SELF,1,"Looks like an iMesh inconsistency, the library shouldn't be messing with this");
#else
  for (dInt i=0; i<ecount; i++) {
    void *ptr = (char*)data + i*iBase_SizeFromType(type);
    dIInt al = peritem,sz = 0;
    iMesh_getEntSetData(mi,esets[i],tag,&ptr,&al,&sz,&ierr);dICHK(mi,ierr);
  }
#endif
  if (size > alloc) dERROR(PETSC_COMM_SELF,1,"Insufficient allocation, iMesh should have thrown an error already");
  dFunctionReturn(0);
}

dErr dMeshGetTaggedSet(dMesh mesh,dMeshTag tag,const void *value,dMeshESH *set)
{
  dIInt ierr,alloc=1,size;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(set,4);
  *set = 0;
  iMesh_getEntSetsByTagsRec(mesh->mi,mesh->root,&tag,value?(const char*const*)&value:NULL,1,0,&set,&alloc,&size,&ierr);dICHK(mesh->mi,ierr);
  if (size != 1) {
    char name[dSTR_LEN],val[dSTR_LEN] = "unknown";
    dIInt ttype;
    iMesh_getTagName(mesh->mi,tag,name,&ierr,sizeof name);dICHK(mesh->mi,ierr);
    iMesh_getTagType(mesh->mi,tag,&ttype,&ierr);dICHK(mesh->mi,ierr);
    switch (ttype) {
    case iBase_INTEGER: snprintf(val,sizeof val,"INTEGER %d",*(int*)value); break;
    case iBase_DOUBLE: snprintf(val,sizeof val,"DOUBLE %e",*(double*)value); break;
    default:
      snprintf(val,sizeof val,"UNKNOWN");
    }
    dERROR(PETSC_COMM_SELF,PETSC_ERR_LIB,"Got %d sets when searching for tag %s, but expected exactly one",size,val);
  }
  dFunctionReturn(0);
}

/* free sets with dFree() */
dErr dMeshGetTaggedSets(dMesh mesh,dMeshTag tag,dInt nvalues,const void *values,dInt *nsets,dMeshESH **sets)
{
  const char *const*pvals = values ? (const char*const*)&values : NULL;
  dIInt ierr,alloc=0,size;
  dMeshESH *isets = NULL;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (nvalues > 0) dValidPointer(values,4);
  if (nvalues > 1) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"Can only search for one value or all values");
  dValidPointer(nsets,5);
  dValidPointer(sets,6);
  *nsets = -1;
  *sets = NULL;
  iMesh_getEntSetsByTagsRec(mesh->mi,mesh->root,&tag,nvalues?pvals:NULL,1,0,&isets,&alloc,&size,&ierr);dICHK(mesh->mi,ierr);
  *nsets = size;
  err = dMallocA(*nsets,sets);dCHK(err);
  err = dMemcpy(*sets,isets,*nsets*sizeof isets[0]);dCHK(err);
  free(isets);
  dFunctionReturn(0);
}

dErr dMeshGetNumEnts(dMesh mesh,dMeshESH set,dEntType type,dEntTopology topo,dInt *num)
{
  dIInt ierr,n;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(num,5);
  if (topo == dTOPO_ALL) {
    iMesh_getNumOfType(mesh->mi,set,type,&n,&ierr);dICHK(mesh->mi,ierr);
  } else {
    iMesh_getNumOfTopo(mesh->mi,set,topo,&n,&ierr);dICHK(mesh->mi,ierr);
  }
  *num = n;
  dFunctionReturn(0);
}

dErr dMeshGetEnts(dMesh mesh,dMeshESH set,dEntType type,dEntTopology topo,dMeshEH ents[],dInt esize,dInt *nents)
{
  iBase_EntityHandle *e;
  dIInt ea,es,ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (esize) dValidPointer(ents,5);
  else {
#if defined dUSE_DEBUG // Count the number of entities to confirm that zero really was the right number to expect
    dErr err;
    dInt mycount;
    err = dMeshGetNumEnts(mesh,set,type,topo,&mycount);dCHK(err);
    if (mycount) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"The set has %D qualifying entities, but caller expected none");
#endif
    if (nents) *nents = 0;
    dFunctionReturn(0);
  }
  e = ents; ea = esize;
  iMesh_getEntities(mesh->mi,set,type,topo,&e,&ea,&es,&ierr);dICHK(mesh->mi,ierr);
  if (e != ents || ea != esize) dERROR(PETSC_COMM_SELF,1,"should not happen");
  if (nents) *nents = es;
  dFunctionReturn(0);
}

dErr dMeshGetNumSubsets(dMesh mesh,dMeshESH set,dInt hops,dInt *nsubsets)
{
  dIInt ierr,n;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(nsubsets,4);
  *nsubsets = -1;
  iMesh_getNumEntSets(mesh->mi,set,hops,&n,&ierr);dICHK(mesh->mi,ierr);
  *nsubsets = n;
  dFunctionReturn(0);
}

dErr dMeshGetSubsets(dMesh mesh,dMeshESH set,dInt hops,dMeshESH subsets[],dInt alloc,dInt *size)
{
  dIInt ierr,s;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (alloc) dValidPointer(subsets,4);
  s = 0;
  iMesh_getEntSets(mesh->mi,set,hops,&subsets,&alloc,&s,&ierr);dICHK(mesh->mi,ierr);
  if (size) *size = s;
  dFunctionReturn(0);
}

/** Get entities of every type in a set.
* @param mesh mesh object
* @param set set to get entities from
* @param toff array of length 5, offset of entities of each type
* @param inents pointer to array of entities, free with dFree()
*/
dErr dMeshGetEntsOff(dMesh mesh,dMeshESH set,dInt toff[],dMeshEH **inents)
{
  dMeshEH *ents;
  dInt     n;
  dErr     err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(toff,3);
  dValidPointer(inents,4);
  *inents = 0;
  err = dMeshGetNumEnts(mesh,set,dTYPE_ALL,dTOPO_ALL,&n);dCHK(err);
  err = dMallocA(n,&ents);dCHK(err);
  toff[0] = 0;
  for (dEntType type=dTYPE_VERTEX; type<dTYPE_ALL; type++) {
    dInt used;
    err = dMeshGetEnts(mesh,set,type,dTOPO_ALL,ents+toff[type],n-toff[type],&used);dCHK(err);
    toff[type+1] = toff[type] + used;
  }
  *inents = ents;
  dFunctionReturn(0);
}

/* Create a subset of an existing set, consisting of entities matching the type+topo condition. */
dErr dMeshSetFilterEnts(dMesh mesh,dMeshESH set,dEntType etype,dEntTopology etopo,dMeshESH *subset)
{
  dErr err;
  dMeshSetOrdering ordering;
  dInt nents;
  dMeshEH *ents;

  dFunctionBegin;
  err = dMeshGetNumEnts(mesh,set,etype,etopo,&nents);dCHK(err);
  err = dMallocA(nents,&ents);dCHK(err);
  err = dMeshGetEnts(mesh,set,etype,etopo,ents,nents,NULL);dCHK(err);
  err = dMeshSetGetOrdering(mesh,set,&ordering);dCHK(err);
  err = dMeshSetCreate(mesh,ordering,subset);dCHK(err);
  err = dMeshSetAddEnts(mesh,set,ents,nents);dCHK(err);
  err = dFree(ents);dCHK(err);
  dFunctionReturn(0);
}

/**
* Get the parallel status for an array of entities.  In serial, returns that all are interior, in parallel, sets bits
* for dSTATUS_XXX.
*
* @param mesh mesh containing entities
* @param ents entities
* @param count number of entities
* @param status status array
*
* @return err
*/
dErr dMeshGetStatus(dMesh mesh,const dMeshEH ents[],dInt count,dEntStatus status[])
{
  dMPIInt size;
  dMeshTag tag;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(ents,2);
  dValidPointer(status,4);
  err = MPI_Comm_size(((dObject)mesh)->comm,&size);dCHK(err);
  if (size == 1) {
    for (dInt i=0; i<count; i++) status[i] = 0;
  } else {
    err = dMeshGetTag(mesh,PARALLEL_STATUS_TAG_NAME,&tag);dCHK(err);
    err = dMeshTagGetData(mesh,tag,ents,count,status,count*sizeof(status[0]),dDATA_BYTE);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dMeshGetTopo(dMesh mesh,dInt count,const dMeshEH ents[],dEntTopology topo[])
{
  dIInt ierr,talloc,tsize,*ttopo;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (!count) dFunctionReturn(0);
  dValidPointer(ents,2);
  dValidPointer(topo,4);
  ttopo = (int*)topo; talloc = count;
  iMesh_getEntArrTopo(mesh->mi,ents,count,&ttopo,&talloc,&tsize,&ierr);dICHK(mesh->mi,ierr);
  if (tsize != count) dERROR(PETSC_COMM_SELF,1,"Wrong number of topologies returned");
  dFunctionReturn(0);
}

dErr dMeshTagBcast(dMesh m,dMeshTag tag)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(m,dMESH_CLASSID,1);
  if (m->ops->tagbcast) {
    err = (*m->ops->tagbcast)(m,tag);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dMeshSetCreate(dMesh mesh,dMeshSetOrdering ordering,dMeshESH *inset)
{
  dIInt ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(inset,3);
  *inset = 0;
  iMesh_createEntSet(mesh->mi,ordering,inset,&ierr);dICHK(mesh->mi,ierr);
  dFunctionReturn(0);
}

dErr dMeshSetDestroy(dMesh mesh,dMeshESH set)
{
  dIInt ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  iMesh_destroyEntSet(mesh->mi,set,&ierr);dICHK(mesh->mi,ierr);
  dFunctionReturn(0);
}

dErr dMeshSetAddEnts(dMesh mesh,dMeshESH set,const dMeshEH *ents,dInt nents)
{
  dIInt ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (!nents) dFunctionReturn(0);
  dValidPointer(ents,3);
  iMesh_addEntArrToSet(mesh->mi,ents,nents,set,&ierr);dICHK(mesh->mi,ierr);
  dFunctionReturn(0);
}

dErr dMeshEntClassifyExclusive(dMesh mesh,dMeshEH ent,dInt nsets,const dMeshESH sets[],dInt *member)
{

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(sets,4);
  dValidPointer(member,5);
  *member = -1;
  for (dInt i=0; i<nsets; i++) {
    dIInt ismember,ierr;
    iMesh_isEntContained(mesh->mi,sets[i],ent,&ismember,&ierr);dICHK(mesh->mi,ierr);
    if (ismember) {
      if (*member >= 0) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"sets are not disjoint");
      *member = i;
    }
  }
  if (*member < 0) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"entity missing from all sets");
  dFunctionReturn(0);
}

/*
 * @param nents number of entities
 * @param ents entities to classify
 * @param inodes value associated with each entity
 * @param nsets number of sets
 * @param sets sets to classify entities into
 * @param counts array of length \a nsets to hold the sum of the values in each set
 * @param rstarts relative start for each count
 */
dErr dMeshClassifyCountInt(dMesh mesh,dInt nents,const dMeshEH ents[],const dInt vals[],dInt nsets,const dMeshESH sets[],dInt counts[])
{
  dErr err;

  dFunctionBegin;
  err = dMemzero(counts,nsets*sizeof(counts[0]));dCHK(err);
  for (dInt i=0; i<nents; i++) {
    dInt member;
    err = dMeshEntClassifyExclusive(mesh,ents[i],3,sets,&member);dCHK(err);
    counts[member] += vals[i];
  }
  dFunctionReturn(0);
}

dErr dMeshLoad(dMesh mesh)
{
  iMesh_Instance mi = mesh->mi;
  //dMeshTag arf,afe,orf,ofe;
  dMeshTag emptysettag;
  dMeshESH root;
  //MeshListInt off=MLZ;
  PetscBool flg;
  dIInt ierr;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  err = PetscLogEventBegin(dLOG_MeshLoad,mesh,0,0,0);dCHK(err);
  err = dMeshGetRoot(mesh,&root);dCHK(err);
  switch (mesh->gentype) {
  case dMESHGEN_NONE: {
    PetscMPIInt rank;
    FILE *file;
    err = MPI_Comm_rank(((dObject)mesh)->comm,&rank);dCHK(err);
    if (!rank) {
      file = fopen(mesh->infile,"r");
      if (!file) dERROR(PETSC_COMM_SELF,1,"Could not open %s for reading: %s",mesh->infile,strerror(errno));
      if (fclose(file)) dERROR(PETSC_COMM_SELF,1,"Error closing %s",mesh->infile);
    }
    if (mesh->ops->load) {
      err = (*mesh->ops->load)(mesh);dCHK(err);
    } else {
      dERROR(PETSC_COMM_SELF,1,"No load function set");
    }
  } break;
  case dMESHGEN_BLOCK: {
    dBool do_geom = dFALSE; // @todo make optional
    err = dMeshGenerateBlock(mesh,root,do_geom);dCHK(err);
  } break;
  }

  /* Make sure the sense tag exists, to protect against degenerate cases */
  iMesh_getTagHandle(mi,dTAG_SENSE,&mesh->senseTag,&ierr,sizeof(dTAG_SENSE));
  if (ierr == iBase_TAG_NOT_FOUND) {
    iMesh_createTag(mi,dTAG_SENSE,1,iBase_INTEGER,&mesh->senseTag,&ierr,sizeof(dTAG_SENSE));dICHK(mi,ierr);
  } else dICHK(mi,ierr);

  /* Create an empty set if it doesn't already exist, useful to prevent degenerate cases */
  iMesh_getTagHandle(mi,dTAG_EMPTYSET,&emptysettag,&ierr,sizeof(dTAG_EMPTYSET));
  if (ierr == iBase_TAG_NOT_FOUND) {
    dMeshESH eset;
    iMesh_createTag(mi,dTAG_EMPTYSET,1,iBase_ENTITY_HANDLE,&emptysettag,&ierr,sizeof(dTAG_EMPTYSET));dICHK(mi,ierr);
    iMesh_createEntSet(mi,0,&eset,&ierr);dICHK(mi,ierr);
    iMesh_setEntSetEHData(mi,root,emptysettag,(iBase_EntityHandle)eset,&ierr);dICHK(mi,ierr);
  } else dICHK(mi,ierr);
  iMesh_getEntSetEHData(mi,root,emptysettag,(iBase_EntityHandle*)&mesh->emptyset,&ierr);dICHK(mi,ierr);

  /* Get all entities of each type. */
  iMesh_getEntities(mi,root,iBase_REGION,iMesh_ALL_TOPOLOGIES,&mesh->r.v,&mesh->r.a,&mesh->r.s,&ierr);dICHK(mi,ierr);
  iMesh_getEntities(mi,root,iBase_FACE,iMesh_ALL_TOPOLOGIES,&mesh->f.v,&mesh->f.a,&mesh->f.s,&ierr);dICHK(mi,ierr);
  iMesh_getEntities(mi,root,iBase_EDGE,iMesh_ALL_TOPOLOGIES,&mesh->e.v,&mesh->e.a,&mesh->e.s,&ierr);dICHK(mi,ierr);
  iMesh_getEntities(mi,root,iBase_VERTEX,iMesh_ALL_TOPOLOGIES,&mesh->v.v,&mesh->v.a,&mesh->v.s,&ierr);dICHK(mi,ierr);
#if 0
  /* Get tags for custom adjacencies, needed since our meshes are nonconforming with respect to the adjacent lower dim entity */
  iMesh_getTagHandle(mi,dTAG_ADJ_REGION_FACE,&arf,&err,strlen(dTAG_ADJ_REGION_FACE));dICHK(mi,err);
  iMesh_getTagHandle(mi,dTAG_ADJ_FACE_EDGE,&afe,&err,strlen(dTAG_ADJ_FACE_EDGE));dICHK(mi,err);
  iMesh_getTagHandle(mi,dTAG_ORIENT_REGION_FACE,&orf,&err,strlen(dTAG_ORIENT_REGION_FACE));dICHK(mi,err);
  iMesh_getTagHandle(mi,dTAG_ORIENT_FACE_EDGE,&ofe,&err,strlen(dTAG_ORIENT_FACE_EDGE));dICHK(mi,err);
  /* Get full adjacencies */
  iMesh_getEHArrData(mi,mesh->r.v,mesh->r.s,arf,&mesh->arf.v,&mesh->arf.a,&mesh->arf.s,&err);dICHK(mi,err); /* region -> face */
  iMesh_getEHArrData(mi,mesh->f.v,mesh->f.s,afe,&mesh->afe.v,&mesh->afe.a,&mesh->afe.s,&err);dICHK(mi,err); /* face -> edge */
  iMesh_getEntArrAdj(mi,mesh->e.v,mesh->e.s,iBase_VERTEX,&mesh->aev.v,&mesh->aev.a,&mesh->aev.s,&off.v,&off.a,&off.s,&err);dICHK(mi,err); /* edge -> vertex */
  MeshListFree(off);      /* We don't use the offsets because we know there are always exactly two vertices per edge. */
  /* Get orientation of lower dimensional entities, we don't need vertex orientation */
  iMesh_getArrData(mi,mesh->r.v,mesh->r.s,orf,&mesh->orf.v,&mesh->orf.a,&mesh->orf.s,&err);dICHK(mi,err); /* region[face] */
  iMesh_getArrData(mi,mesh->f.v,mesh->f.s,ofe,&mesh->ofe.v,&mesh->ofe.a,&mesh->ofe.s,&err);dICHK(mi,err); /* face[edge] */
#endif

  err = PetscLogEventEnd(dLOG_MeshLoad,mesh,0,0,0);dCHK(err);

  /* View if requested */
  err = PetscOptionsHasName(((dObject)mesh)->prefix,"-dmesh_view",&flg);dCHK(err);
  if (flg) {
    dViewer viewer;
    err = PetscViewerASCIIGetStdout(((dObject)mesh)->comm,&viewer);dCHK(err);
    err = dMeshView(mesh,viewer);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dMeshInferIntermediateAdjacencies(dMesh mesh)
{
  dIInt ierr,adjtable_arr[dTYPE_ALL*dTYPE_ALL],*table = adjtable_arr,table_a = dTYPE_ALL*dTYPE_ALL,table_s;

  dFunctionBegin;
  iMesh_getAdjTable(mesh->mi,&table,&table_a,&table_s,&ierr);dICHK(mesh->mi,ierr);
  if (table_a != table_s) dERROR(PETSC_COMM_SELF,PETSC_ERR_LIB,"iMesh changed the size of the adjacency table");
  table[dTYPE_EDGE*dTYPE_ALL + dTYPE_EDGE] = iBase_AVAILABLE;
  table[dTYPE_FACE*dTYPE_ALL + dTYPE_FACE] = iBase_AVAILABLE;
  iMesh_setAdjTable(mesh->mi,table,table_s,&ierr);dICHK(mesh->mi,ierr);
  dFunctionReturn(0);
}

dErr dMeshGetInstance(dMesh m,iMesh_Instance *mi)
{

  dFunctionBegin;
  dValidHeader(m,dMESH_CLASSID,1);
  dValidPointer(mi,2);
  *mi = m->mi;
  dFunctionReturn(0);
}

dErr dMeshDestroy(dMesh *meshp)
{
  dMesh m = *meshp;
  dErr err;
  dIInt ierr;

  dFunctionBegin;
  if (!m) dFunctionReturn(0);
  PetscValidHeaderSpecific(m,dMESH_CLASSID,1);
  if (--((PetscObject)m)->refct > 0) dFunctionReturn(0);
  if (m->ops->destroy) {
    err = (*m->ops->destroy)(m);dCHK(err);
  }
  MeshListFree(m->v);   MeshListFree(m->e);   MeshListFree(m->f);   MeshListFree(m->r);
  MeshListFree(m->arf); MeshListFree(m->afe); MeshListFree(m->aev);
  MeshListFree(m->orf); MeshListFree(m->ofe);
  MeshListFree(m->x);
  iMesh_dtor(m->mi,&ierr);dICHK(m->mi,ierr);
  err = PetscFree(m->infile);dCHK(err);
  err = PetscFree(m->inoptions);dCHK(err);
#ifdef dHAVE_ITAPS_REL
  if (m->irel)  {iRel_destroy(m->irel,&ierr);dIRCHK(m->irel,ierr);}
  if (m->igeom) {iGeom_dtor(m->igeom,&ierr);dIGCHK(m->igeom,ierr);}
#endif
  err = PetscHeaderDestroy(meshp);dCHK(err);
  dFunctionReturn(0);
}


/**
* Creates a \p dRule tag over all non-vertex topological entities in the mesh.  Also tags the root entity set with the
* given tag name and value equal to the base pointer for the rule storage.  This tag should be removed using
* dMeshDestroyRuleTag().
*
* @param mesh mesh
* @param set entity set handle on which to add the tag, tags all non-vertex entities in \a set
* @param name unique identifier for the tag (mostly for debugging)
* @param degree polynomial degree which should be integrated exactly when the element has an affine map
* @param[out] inrtag tag
*
* @return err
*/
dErr dMeshCreateRuleTagIsotropic(dMesh mesh,dMeshESH set,const char name[],dInt degree,dMeshTag *inrtag)
{
  dMeshTag rtag;
  dMeshEH *ents;
  dInt nents,toff[dTYPE_ALL+1];
  dPolynomialOrder *rdeg;
  dEntTopology *topo;
  dEntType firstType;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(inrtag,2);
  *inrtag = 0;
  //err = dMeshTagCreate(mesh,name,sizeof(dRule),dDATA_BYTE,&rtag);dCHK(err);
  err = dMeshTagCreate(mesh,name,1,dDATA_INT,&rtag);dCHK(err);

  firstType = dTYPE_VERTEX;  /* Hack: Get vertices as well (we should not need to label them) */
  toff[dTYPE_VERTEX] = toff[dTYPE_EDGE] = 0;
  for (dEntType type=firstType; type<dTYPE_ALL; type++) {
    dInt t;
    err = dMeshGetNumEnts(mesh,set,type,dTOPO_ALL,&t);dCHK(err);
    toff[type+1] = toff[type] + t;
  }
  nents = toff[dTYPE_ALL];
  err = dMallocA3(nents,&ents,nents,&topo,nents,&rdeg);dCHK(err);
  for (dEntType type=firstType; type<dTYPE_ALL; type++) {
    err = dMeshGetEnts(mesh,set,type,dTOPO_ALL,ents+toff[type],toff[type+1]-toff[type],NULL);dCHK(err);
  }
  err = dMeshGetTopo(mesh,nents,ents,topo);dCHK(err);
  for (dInt i=0; i<nents; i++) {
    switch (topo[i]) {
      case dTOPO_POINT: rdeg[i] = dPolynomialOrderCreate(0,0,0,0); break;
      case dTOPO_LINE: rdeg[i] = dPolynomialOrderCreate(0,degree,0,0); break;
      case dTOPO_QUAD: rdeg[i] = dPolynomialOrderCreate(0,degree,degree,0); break;
      case dTOPO_HEX:  rdeg[i] = dPolynomialOrderCreate(0,degree,degree,degree); break;
      default: dERROR(PETSC_COMM_SELF,1,"Topology %d not supported",topo[i]);
    }
  }
  //err = dJacobiGetRule(jac,nents,topo,rdeg,rules);dCHK(err);
  //err = dMeshTagSetData(mesh,rtag,ents,nents,rules,nents*(dInt)sizeof(s_dRule),dDATA_BYTE);dCHK(err);
  err = dMeshTagSetData(mesh,rtag,ents,nents,rdeg,nents,dDATA_INT);dCHK(err);
  err = dFree3(ents,topo,rdeg);dCHK(err);
  *inrtag = rtag;
  dFunctionReturn(0);
}

dErr dMeshDestroyRuleTag(dMesh mesh,dMeshTag rtag)
{
  iMesh_Instance mi = mesh->mi;
  void *base;
  int bsize,balloc=sizeof(void*);
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  iMesh_getEntSetData(mi,mesh->root,rtag,&base,&balloc,&bsize,&err);dICHK(mi,err);
  err = dFree(base);dCHK(err);
  dFunctionReturn(0);
}

static dErr dMeshAdjacencyPermutations_Private(dMeshAdjacency ma,const dInt connoff[],const dMeshEH conn[])
{
  dInt i,e,ai,aind;
  dErr err;

  dFunctionBegin;
  for (e=ma->toff[dTYPE_EDGE]; e<ma->toff[dTYPE_ALL]; e++) { /* All non-vertex entities */
    switch (ma->topo[e]) {
      case dTOPO_HEX:
        for (i=0; i<6; i++) {
          ai = ma->adjoff[e]+i; aind = ma->adjind[ai];
          err = dGeomOrientFindPerm_HexQuad(conn+connoff[e],conn+connoff[aind],i,&ma->adjperm[ai]);dCHK(err);
        }
        break;
      case dTOPO_QUAD:
        for (i=0; i<4; i++) {
          ai = ma->adjoff[e]+i; aind = ma->adjind[ai];
          err = dGeomOrientFindPerm_QuadLine(conn+connoff[e],conn+connoff[aind],i,&ma->adjperm[ai]);dCHK(err);
        }
        break;
      case dTOPO_LINE:
        for (i=0; i<2; i++) { /* Both endpoints are vertices, they always have degree 1 */
          ai = ma->adjoff[e]+i; aind = ma->adjind[ai];
          ma->adjperm[ai] = 0;  /* Vertices cannot be permuted */
        }
        break;
      default: dERROR(PETSC_COMM_SELF,1,"Topology %s not supported",dMeshEntTopologyName(ma->topo[e]));
    }
  }
  dFunctionReturn(0);
}

/**
@note Not collective
*/
dErr dMeshGetAdjacency(dMesh mesh,dMeshESH set,dMeshAdjacency *inadj)
{
  iMesh_Instance mi = mesh->mi;
  struct _p_dMeshAdjacency ma;
  dMeshEH *adj,*conn;
  dEntType type;
  dInt i,cnt,nadj,*connoff,tnents,*eind;
  dIInt ierr;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(inadj,3);
  *inadj = 0;

  ma.set = set;                 /* Store set used to check out adjacencies */
  /*
  * Step 1: number all entities in \a set
  */

  err = dMeshGetNumEnts(mesh,set,dTYPE_ALL,dTOPO_ALL,&ma.nents);dCHK(err);
  /* These arrays are persistant for the life of dMeshAdjacency, the arrays are populated in the remainder of this function */
  err = dMallocA3(ma.nents,&ma.ents,ma.nents+1,&ma.adjoff,ma.nents,&ma.topo);dCHK(err);

  /* Get counts and populate \c ma.ents, guaranteed to be in type-order */
  ma.toff[dTYPE_VERTEX] = 0; cnt = 0;
  for (type=dTYPE_VERTEX; type<dTYPE_ALL; type++) {
    err = dMeshGetEnts(mesh,set,type,iMesh_ALL_TOPOLOGIES,ma.ents+cnt,ma.nents-cnt,&tnents);dCHK(err);
    ma.toff[type] = cnt;
    cnt += tnents;
  }
  ma.toff[dTYPE_ALL] = cnt;
  if (cnt != ma.nents) {         /* Check whether the iMesh implementation is sane */
    dMeshEH *allents;
    err = dMallocA(ma.nents,&allents);dCHK(err);
    err = dMeshGetEnts(mesh,set,dTYPE_ALL,dTOPO_ALL,allents,ma.nents,NULL);dCHK(err);
    for (i=0; i<ma.nents; i++) {
      if (ma.ents[i] != allents[i]) {
        dERROR(PETSC_COMM_SELF,1,"mismatch: ents[%d]=%p  allents[%d]=%p\n",i,ma.ents[i],i,allents[i]);
      }
    }
    dERROR(PETSC_COMM_SELF,1,"count by type %d does not agree with total count %d",cnt,ma.nents);
  }

  /* Populate \c ma.topo */
  err = dMeshGetTopo(mesh,ma.nents,ma.ents,ma.topo);dCHK(err);

  /* Set indices (into \c ma.ents[]) for all entities */
  err = dMallocA(ma.nents,&eind);dCHK(err);
  for (i=0; i<ma.nents; i++) {
    eind[i] = i;
  }

  /* Create the tag to hold indices and set it with strictly increasing values */
  err = dMeshTagCreateTemp(mesh,"index",1,dDATA_INT,&ma.indexTag);dCHK(err);
  err = dMeshTagSetData(mesh,ma.indexTag,ma.ents,ma.nents,eind,ma.nents,dDATA_INT);dCHK(err);
  err = dFree(eind);dCHK(err);

  /*
  * Step 2: use connectivity and indices for all adjacent entities
  */

  { /* connectivity for all entities, vertices have null connectivity */
    MeshListEH  ml_conn=MLZ;
    MeshListInt ml_connoff=MLZ;
    iMesh_getEntArrAdj(mi,ma.ents+ma.toff[dTYPE_EDGE],ma.nents-ma.toff[dTYPE_EDGE],iBase_VERTEX,MLREF(ml_conn),MLREF(ml_connoff),&ierr);dICHK(mi,ierr);
    err = dMallocA2(ml_connoff.s+ma.toff[dTYPE_EDGE],&connoff,ml_conn.s,&conn);dCHK(err);
    for (i=0; i<ma.toff[dTYPE_EDGE]; i++) { connoff[i] = 0; } /* empty connectivity for vertices */
    err = dMemcpy(conn,ml_conn.v,ml_conn.s*sizeof(*conn));dCHK(err);
    err = dMemcpy(connoff+ma.toff[dTYPE_EDGE],ml_connoff.v,ml_connoff.s*sizeof(*connoff));dCHK(err);
    MeshListFree(ml_conn); MeshListFree(ml_connoff);dCHK(err);
  }

  { /* Downward adjacency for all entities, vertices have no downward adjacencies */
    MeshListEH  ml_adj[4]={MLZ,MLZ,MLZ,MLZ};
    MeshListInt ml_adjoff[4]={MLZ,MLZ,MLZ,MLZ};
    dInt adjcnt;
    for (type=dTYPE_EDGE; type<dTYPE_ALL; type++) {
      iMesh_getEntArrAdj(mi,ma.ents+ma.toff[type],ma.toff[type+1]-ma.toff[type],type-1,MLREF(ml_adj[type]),MLREF(ml_adjoff[type]),&ierr);dICHK(mi,ierr);
    }
    for (type=dTYPE_EDGE; type<dTYPE_ALL; type++) {
      if (ml_adjoff[type].s != ma.toff[type+1]-ma.toff[type]+1) dERROR(PETSC_COMM_SELF,1,"unexpected number of adjacent offsets");
    }
    nadj = ml_adj[1].s+ml_adj[2].s+ml_adj[3].s; /* total number of adjacent entities */
    if (!nadj) dERROR(PETSC_COMM_SELF,1,"No adjacent entities, seems like a deficient mesh");
    err = dMallocA(nadj,&adj);                  /* Freed after getting index tag values */
    err = dMallocA2(nadj,&ma.adjind,nadj,&ma.adjperm);dCHK(err);      /* Persistant */
    for (i=0; i<ma.toff[dTYPE_EDGE]+1; i++) { ma.adjoff[i] = 0; }     /* vertices have no adjacent entities */
    adjcnt = 0;                                                       /* Finger in the adjacent entities array */
    for (type=dTYPE_EDGE; type<dTYPE_ALL; type++) {
      for (i=ma.toff[type]; i<ma.toff[type+1]+1; i++) {
        ma.adjoff[i] = ma.adjoff[ma.toff[type]] + ml_adjoff[type].v[i-ma.toff[type]];
      }
      for (i=0; i<ml_adj[type].s; i++) { /* Pack the adjacent entities */
        adj[adjcnt++] = ml_adj[type].v[i];
      }
    }
    if (adjcnt != nadj) dERROR(PETSC_COMM_SELF,1,"unexpected adjacent entity count");
    for (i=0; i<4; i++) {
      MeshListFree(ml_adj[i]);
      MeshListFree(ml_adjoff[i]);
    }
  }

  err = dMeshTagGetData(mesh,ma.indexTag,adj,nadj,ma.adjind,nadj,dDATA_INT);dCHK(err); /* get indices of adj ents */
  err = dFree(adj);dCHK(err);

  /* Determine permutation of adjacent entities */
  err = dMeshAdjacencyPermutations_Private(&ma,connoff,conn);dCHK(err);
#if defined(dMESHADJACENCY_HAS_CONNECTIVITY)
  ma.connoff = connoff; ma.conn = conn;
#else
  err = dFree2(connoff,conn);dCHK(err);
#endif
  err = dMalloc(sizeof(ma),inadj);dCHK(err);
  err = dMemcpy(*inadj,&ma,sizeof(ma));dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshRestoreAdjacency(dMesh dUNUSED mesh,dMeshESH set,dMeshAdjacency *inma)
{
  dErr err;
  dMeshAdjacency ma = *inma;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(inma,3);
  ma = *inma;
  if (set != ma->set) dERROR(PETSC_COMM_SELF,1,"Adjacency for the wrong set");
  err = dMeshTagDestroy(mesh,ma->indexTag);dCHK(err);
  err = dFree3(ma->ents,ma->adjoff,ma->topo);dCHK(err);
  err = dFree2(ma->adjind,ma->adjperm);dCHK(err);
#if defined(dMESHADJACENCY_HAS_CONNECTIVITY)
  err = dFree2(ma->connoff,ma->conn);dCHK(err);
#endif
  err = dFree(*inma);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshGetVertexCoords(dMesh mesh,dInt nvtx,const dMeshEH vtx[],const dReal **incoords)
{
  dIInt ierr;
  MeshListReal vcoords=MLZ;
  dReal *coords;
  dErr err;

  dFunctionBegin;
  iMesh_getVtxArrCoords(mesh->mi,vtx,nvtx,iBase_INTERLEAVED,MLREF(vcoords),&ierr);dICHK(mesh->mi,ierr);
  err = dMallocA(vcoords.s,&coords);dCHK(err);
  for (dInt i=0; i<vcoords.s; i++) coords[i] = (dReal)vcoords.v[i];
  *incoords = coords;
  dFunctionReturn(0);
}

dErr dMeshRestoreVertexCoords(dMesh dUNUSED mesh,dInt dUNUSED nvtx,const dMeshEH dUNUSED vtx[],const dReal **incoords)
{
  dErr err;

  dFunctionBegin;
  err = dFree(*incoords);dCHK(err);
  dFunctionReturn(0);
}

/** Get vertex coordinates for the vertices representing the connectivity of the entities.
*
* @example <c>x+xoff[i] ... x+xoff[i+1]</c> are the coordinates of the vertices in the connectivity of \c ents[i]
*
* @param mesh Mesh object
* @param n number of entities
* @param ents entity handles
* @param inxoff address of offsets of vertices for each entity, NOT the offset of the first vertex coordinate because \a x has array type
* @param inx address of vertex values
*/
dErr dMeshGetAdjVertexCoords(dMesh mesh,dInt n,const dMeshEH ents[],const dInt **inxoff,const dReal **inx)
{
  iMesh_Instance mi      = mesh->mi;
  MeshListEH     conn    = MLZ;
  MeshListInt    connoff = MLZ;
  MeshListReal   vtx     = MLZ;
  dIInt          ierr;
  dReal          *x;
  dInt           *xoff;
  dErr           err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(ents,3);
  dValidPointer(inxoff,4);
  dValidPointer(inx,5);
  *inxoff = NULL; *inx = NULL;
  iMesh_getEntArrAdj(mi,ents,n,iBase_VERTEX,MLREF(conn),MLREF(connoff),&ierr);dICHK(mi,ierr);
  iMesh_getVtxArrCoords(mi,conn.v,conn.s,iBase_INTERLEAVED,MLREF(vtx),&ierr);dICHK(mi,ierr);

  /* Allocate for coordinates \a x first to be sure they are suitable aligned (dReal has stricter alignment requirements
  * than dInt) */
  err = dMallocA2(connoff.v[n]*3,&x,n+1,&xoff);dCHK(err);
  for (dInt i=0; i<n; i++) {
    xoff[i] = connoff.v[i];     /* Leave the offset as number of vertices */
    for (dInt j=3*connoff.v[i]; j<3*connoff.v[i+1]; j++) {
      x[j] = vtx.v[j];
    }
  }
  xoff[n] = connoff.v[n];       /* Leave the offset as number of vertices, not first coordinate of each vtx */
  MeshListFree(conn); MeshListFree(connoff); MeshListFree(vtx);
  *inxoff = xoff;
  *inx = x;
  dFunctionReturn(0);
}

/*
* Since the vertex coords are often persistent for the life of the dFS, it's common that \a ents will not be available
* when this function is called, hence we accept a NULL argument.  This is a hack to preserve a symmetric interface.
**/
dErr dMeshRestoreAdjVertexCoords(dMesh dUNUSED mesh,dUNUSED dInt n,const dUNUSED dMeshEH ents[],const dInt **inxoff,const dReal **inx)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidPointer(inxoff,4);
  dValidPointer(inx,5);
  err = dFree2(*inx,*inxoff);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshPartitionOnOwnership(dMesh mesh,dMeshEH ents[],dInt n,dInt *ghstart)
{
  dEntStatus *status;
  dMeshEH    *pents;
  dInt        owned;
  dErr        err;

  dFunctionBegin;
  err = dMallocA2(n,&status,n,&pents);dCHK(err);
  err = dMeshGetStatus(mesh,ents,n,status);dCHK(err);
  owned = 0;
  for (dInt i=0; i<n; i++) owned += !!(~status[i] & dSTATUS_UNOWNED);
  *ghstart = owned;
  for (dInt i=0,o=0,g=owned; i<n; i++) {
    if (status[i] & dSTATUS_UNOWNED) pents[g++] = ents[i];
    else                             pents[o++] = ents[i];
  }
  err = dMemcpy(ents,pents,n*sizeof ents[0]);dCHK(err);
  err = dFree2(status,pents);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshMorph(dMesh mesh,void (*morph)(void*,double*),void *ctx)
{
  iMesh_Instance mi;
  dMeshEH       *ents=0;
  double        *coords=0;
  dIInt          ierr,ents_a=0,ents_s,coords_a=0,coords_s;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  dValidFunctionPointer(morph,2);   /* It's actually a function pointer */
  dValidPointer(ctx,3);
  mi = mesh->mi;
  iMesh_getEntities(mi,mesh->root,iBase_VERTEX,iMesh_ALL_TOPOLOGIES,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
  iMesh_getVtxArrCoords(mi,ents,ents_s,iBase_INTERLEAVED,&coords,&coords_a,&coords_s,&ierr);dICHK(mi,ierr);
  for (dInt i=0; i<coords_s; i+=3) morph(ctx,coords+i);
  iMesh_setVtxArrCoords(mi,ents,ents_s,iBase_INTERLEAVED,coords,coords_s,&ierr);dICHK(mi,ierr);
  free(ents);
  free(coords);
  dFunctionReturn(0);
}

// Mutates a set so that it is its own closure
dErr dMeshSetClosure(dMesh mesh,dMeshESH set)
{
  dErr err;
  dIInt ierr,ents_a,ents_s,adj_a,adj_s,*off,off_s,off_a;
  dMeshEH *ents,*adj;
  iMesh_Instance mi;
  dMeshSetOrdering ordering;

  dFunctionBegin;
  err = dMeshSetGetOrdering(mesh,set,&ordering);dCHK(err);
  if (ordering != dMESHSET_UNORDERED) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"This function only works for unordered sets");
  err = dMeshGetInstance(mesh,&mi);dCHK(err);

  ents_a = adj_a = off_a = 0;
  iMesh_getEntitiesRec(mi,set,dTYPE_FACE,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr); // Recursive closure because this method is sort of specialized for boundary sets :-(
  iMesh_getEntArrAdj(mi,ents,ents_s,dTYPE_EDGE,&adj,&adj_a,&adj_s,&off,&off_a,&off_s,&ierr);dICHK(mi,ierr);
  err = dMeshSetAddEnts(mesh,set,adj,adj_s);dCHK(err);
  free(adj); free(off); free(ents);
  adj = NULL; off = NULL; ents = NULL;

  ents_a = adj_a = off_a = 0;
  iMesh_getEntitiesRec(mi,set,dTYPE_EDGE,dTOPO_ALL,1,&ents,&ents_a,&ents_s,&ierr);dICHK(mi,ierr);
  iMesh_getEntArrAdj(mi,ents,ents_s,dTYPE_VERTEX,&adj,&adj_a,&adj_s,&off,&off_a,&off_s,&ierr);dICHK(mi,ierr);
  err = dMeshSetAddEnts(mesh,set,adj,adj_s);dCHK(err);
  free(adj); free(off); free(ents);
  adj = NULL; off = NULL; ents = NULL;
  ents_a = adj_a = off_a = 0;
  dFunctionReturn(0);
}

#ifdef dHAVE_ITAPS_REL
// The mesh takes ownership with this call. ITAPS interfaces do not have reference counting.
dErr dMeshSetGeometryRelation(dMesh mesh,iGeom_Instance geom,iRel_Instance rel)
{
  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (rel) {
    if (mesh->irel && rel != mesh->irel) dERROR(((dObject)mesh)->comm,PETSC_ERR_ARG_WRONGSTATE,"Can only set relation once");
    mesh->irel = rel;
  }
  if (geom) {
    if (mesh->igeom && geom != mesh->igeom) dERROR(((dObject)mesh)->comm,PETSC_ERR_ARG_WRONGSTATE,"Can only set geometry once");
    mesh->igeom = geom;
  }
  dFunctionReturn(0);
}
dErr dMeshGetGeometryRelation(dMesh mesh,iGeom_Instance *geom,iRel_Instance *rel)
{
  dFunctionBegin;
  dValidHeader(mesh,dMESH_CLASSID,1);
  if (rel)  *rel = mesh->irel;
  if (geom) *geom = mesh->igeom;
  dFunctionReturn(0);
}
#endif
