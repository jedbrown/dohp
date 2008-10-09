#include "private/dmeshimpl.h"
#include <ctype.h>              /* needed for isprint() */

static dErr dMeshView_EntSet(dMesh m,dMeshESH root,PetscViewer viewer);

static const dInt iBase_SizeFromType[4] = {sizeof(int),sizeof(double),sizeof(void*),sizeof(char)};

const char *const iBase_ErrorString[] = {
  "iBase_SUCCESS",
  "iBase_MESH_ALREADY_LOADED",
  "iBase_NO_MESH_DATA",
  "iBase_FILE_NOT_FOUND",
  "iBase_FILE_WRITE_ERROR",
  "iBase_NIL_ARRAY",
  "iBase_BAD_ARRAY_SIZE",
  "iBase_BAD_ARRAY_DIMENSION",
  "iBase_INVALID_ENTITY_HANDLE",
  "iBase_INVALID_ENTITY_COUNT",
  "iBase_INVALID_ENTITY_TYPE",
  "iBase_INVALID_ENTITY_TOPOLOGY",
  "iBase_BAD_TYPE_AND_TOPO",
  "iBase_ENTITY_CREATION_ERROR",
  "iBase_INVALID_TAG_HANDLE",
  "iBase_TAG_NOT_FOUND",
  "iBase_TAG_ALREADY_EXISTS",
  "iBase_TAG_IN_USE",
  "iBase_INVALID_ENTITYSET_HANDLE",
  "iBase_INVALID_ITERATOR_HANDLE",
  "iBase_INVALID_ARGUMENT",
  "iBase_MEMORY_ALLOCATION_FAILED",
  "iBase_NOT_SUPPORTED",
  "iBase_FAILURE"
};

const char *const iMesh_TopologyName[12] = {
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

const int iMesh_TypeFromTopology[12] = {
  iBase_VERTEX,                 /* POINT */
  iBase_EDGE,                   /* LINE_SEGMENT */
  iBase_FACE,                   /* POLYGON */
  iBase_FACE,                   /* TRIANGLE */
  iBase_FACE,                   /* QUADRILATERAL */
  iBase_REGION,                 /* POLYHEDRON */
  iBase_REGION,                 /* TETRAHEDRON */
  iBase_REGION,                 /* HEXAHEDRON */
  iBase_REGION,                 /* PRISM */
  iBase_REGION,                 /* PYRAMID */
  iBase_REGION,                 /* SEPTAHEDRON */
  iBase_ALL_TYPES,              /* ALL_TOPOLOGIES */
};

const char *const iBase_TypeName[] = {
  "iBase_VERTEX",
  "iBase_EDGE",
  "iBase_FACE",
  "iBase_REGION",
  "iBase_ALL_TYPES"
};

const char *const iBase_TagValueTypeName[] = {
  "iBase_INTEGER",
  "iBase_DOUBLE",
  "iBase_ENTITY_HANDLE",
  "iBase_BYTES"
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
      dERROR(1,"Orientation not supported.");
  }
  dFunctionReturn(0);
}

dErr dMeshOrientLoopBounds_Line(dInt orient, const dInt *size, DohpLoopBounds *l)
{

  dFunctionBegin;
  switch (orient) {
    case 0: l->start = 0;         l->stride = 1;  l->end = size[0]; break;
    case 1: l->start = size[0]-1; l->stride = -1; l->end = -1;       break;
    default: dERROR(1,"Orientation not supported.");
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
      dERROR(1,"Face number not recognized.");
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
      dERROR(1,"Edge number not recognized.");
  }
  dFunctionReturn(0);
}

#if 0
/* Maps facet degrees of freedom to element degrees of freedom, adding
* contributions.  This function is actually an optimization for conforming
* elements since it does not need to do interpolation. */
static dErr EFSFacetToElem_HexQuad_Conforming(dInt dof,const dInt rsize[],const dInt fsize[],dInt fnum,dInt forient,const dScalar fvals[],dScalar rvals[])
{
  dInt ri,rj,fi,fj,k;
  DohpLoopBounds rl[2],fl[2];
  dErr err;

  dFunctionBegin;
  err = dMeshLoopBounds_Hex(rsize,fnum,rl);dCHK(err);
  err = dMeshOrientLoopBounds_Quad(forient,fsize,fl);dCHK(err);
  for (ri=rl[0].start,fi=fl[0].start; ri!=rl[0].end && fi!=fl[0].end; ri+=rl[0].stride,fi+=fl[0].stride) {
    for (rj=rl[1].start,fj=fl[1].start; rj!=rl[1].end && fj!=fl[1].end; rj+=rl[1].stride,fj+=fl[1].stride) {
      for (k=0; k<dof; k++) {
        rvals[(ri+rj)*dof+k] += fvals[(fi+fj)*dof+k];
      }
    }
    if (!(rj==rl[1].end && fj==fl[1].end)) {
      dERROR(1,"Inner loop bounds do not agree.  Is this relation conforming?");
    }
  }
  if (!(ri==rl[0].end && fi==fl[0].end)) {
    dERROR(1,"Outer loop bounds do not agree.  Is this relation conforming?");
  }
  dFunctionReturn(0);
}

static dErr EFSFacetToElem_QuadLine_Conforming(dInt dof,const dInt fsize[],const dInt esize[],dInt en,dInt eorient,const dScalar evals[],dScalar fvals[])
{
  dInt fi,ei,j;
  DohpLoopBounds fl,el;
  dErr err;

  dFunctionBegin;
  err = dMeshLoopBounds_Quad(fsize,en,&fl);dCHK(err);
  err = dMeshOrientLoopBounds_Line(eorient,esize,&el);dCHK(err);
  for (fi=fl.start,ei=el.start; fi!=fl.end && ei!=el.end; fi+=fl.stride,ei+=el.stride) {
    for (j=0; j<dof; j++) {
      fvals[fi*dof+j] += evals[ei*dof+j];
    }
  }
  if (!(fi==fl.end && ei==el.end)) {
    dERROR(1,"Loop bounds do not agree.  Is this relation conforming?");
  }
  dFunctionReturn(0);
}

#endif

dErr dMeshSetInFile(dMesh mesh,const char *fname,const char *options)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  if (fname) {
    err = PetscStrfree(mesh->infile);dCHK(err);
    err = PetscStrallocpy(fname,&mesh->infile);dCHK(err);
  }
  if (options) {
    err = PetscStrfree(mesh->inoptions);dCHK(err);
    err = PetscStrallocpy(options,&mesh->inoptions);dCHK(err);
  }
  dFunctionReturn(0);
}


/**
* This function allocates memory for the tag.  It should be freed with PetscFree()
* 
*/
dErr dMeshGetTagName(dMesh mesh,dMeshTag tag,char **name)
{
  iMesh_Instance mi = mesh->mi;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(mesh,dMESH_COOKIE,1);
  dValidPointer(name,2);
  err = PetscMalloc(dNAME_LEN,name);dCHK(err);
  iMesh_getTagName(mi,tag,*name,&err,dNAME_LEN);dICHK(mi,err);
  dFunctionReturn(0);
}

dErr dMeshTagCreateTemp(dMesh mesh,const char template[],dInt count,dDataType type,dMeshTag *intag)
{
  static dInt unique_id = 0;
  char name[dNAME_LEN];
  dMeshTag tag;
  dIInt ierr,itype;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  dValidPointer(intag,4);
  *intag = 0;
  if (count > 0) {
    err = PetscSNPrintf(name,sizeof(name),"TEMP_%s_%d",template,unique_id++);dCHK(err);
    err = dDataTypeToITAPS(type,&itype);dCHK(err);
    iMesh_createTag(mesh->mi,name,count,itype,&tag,&ierr,sizeof(name));dICHK(mesh->mi,ierr);
    *intag = tag;
  }
  dFunctionReturn(0);
}

dErr dMeshTagDestroy(dMesh mesh,dMeshTag tag)
{
  dIInt ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  iMesh_destroyTag(mesh->mi,tag,1,&ierr);dICHK(mesh->mi,ierr);
  dFunctionReturn(0);
}

dErr dMeshTagSetData(dMesh mesh,dMeshTag tag,const dMeshEH ents[],dInt ecount,const void *data,dInt count,dDataType type)
{
  iMesh_Instance mi = mesh->mi;
  dMeshEH *ments;
  const char *dptr = data;
  dIInt size,ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  dValidPointer(ents,3);
  dValidPointer(data,5);
  ments = (dMeshEH*)ents;       /* Cast away const, pretty sure iMesh will never modify the handles */
  size = count * iBase_SizeFromType[type];
  iMesh_setArrData(mi,ments,ecount,tag,dptr,size,&ierr);dICHK(mi,ierr);
  dFunctionReturn(0);
}

dErr dMeshTagGetData(dMesh mesh,dMeshTag tag,const dMeshEH ents[],dInt ecount,void *data,dInt count,dDataType type)
{
  iMesh_Instance mi = mesh->mi;
  char *dptr = data;
  dIInt size,alloc,ierr;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  dValidPointer(ents,3);
  dValidPointer(data,5);
  alloc = count * iBase_SizeFromType[type];
  iMesh_getArrData(mi,ents,ecount,tag,&dptr,&alloc,&size,&ierr);dICHK(mi,ierr);
  if (dptr != (char*)data || alloc != count * iBase_SizeFromType[type] || size != alloc)
    dERROR(1,"Looks like an iMesh inconsistency, the library shouldn't be messing with this");
  dFunctionReturn(0);
}  

dErr dMeshGetNumEnts(dMesh mesh,dMeshESH set,dEntType type,dEntTopology topo,dInt *num)
{
  dIInt ierr,n;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  dValidPointer(num,5);
  if (topo == iMesh_ALL_TOPOLOGIES) {
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
  dValidHeader(mesh,dMESH_COOKIE,1);
  dValidPointer(ents,5);
  e = ents; ea = esize;
  iMesh_getEntities(mesh->mi,set,type,topo,&e,&ea,&es,&ierr);dICHK(mesh->mi,ierr);
  if (e != ents || ea != esize) dERROR(1,"should not happen");
  if (nents) *nents = es;
  dFunctionReturn(0);
}

dErr dMeshGetTopo(dMesh mesh,const dMeshEH ents[],dInt count,dEntTopology topo[])
{
  dIInt ierr,talloc,tsize,*ttopo;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  if (!count) dFunctionReturn(0);
  dValidPointer(ents,2);
  dValidPointer(topo,4);
  ttopo = (int*)topo; talloc = count;
  iMesh_getEntArrTopo(mesh->mi,ents,count,&ttopo,&talloc,&tsize,&ierr);dICHK(mesh->mi,ierr);
  if (tsize != count) dERROR(1,"Wrong number of topologies returned");
  dFunctionReturn(0);
}

dErr dMeshTagBcast(dMesh m,dMeshTag tag)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(m,dMESH_COOKIE,1);
  if (m->ops->tagbcast) {
    err = (*m->ops->tagbcast)(m,tag);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dMeshLoad(dMesh mesh)
{
  iMesh_Instance mi = mesh->mi;
  //dMeshTag arf,afe,orf,ofe;
  dMeshESH root;
  //MeshListInt off=MLZ;
  dBool flg;
  dErr err;

  dFunctionBegin;
  if (mesh->ops->load) {
    err = (*mesh->ops->load)(mesh);dCHK(err);
  } else {
    dERROR(1,"No load function set");
  }
  iMesh_getRootSet(mi,&root,&err);dICHK(mi,err);
  mesh->root = root;

  /* Get all entities of each type. */
  iMesh_getEntities(mi,root,iBase_REGION,iMesh_ALL_TOPOLOGIES,&mesh->r.v,&mesh->r.a,&mesh->r.s,&err);ICHKERRQ(mi,err);
  iMesh_getEntities(mi,root,iBase_FACE,iMesh_ALL_TOPOLOGIES,&mesh->f.v,&mesh->f.a,&mesh->f.s,&err);ICHKERRQ(mi,err);
  iMesh_getEntities(mi,root,iBase_EDGE,iMesh_ALL_TOPOLOGIES,&mesh->e.v,&mesh->e.a,&mesh->e.s,&err);ICHKERRQ(mi,err);
  iMesh_getEntities(mi,root,iBase_VERTEX,iMesh_ALL_TOPOLOGIES,&mesh->v.v,&mesh->v.a,&mesh->v.s,&err);ICHKERRQ(mi,err);
#if 0
  /* Get tags for custom adjacencies, needed since our meshes are nonconforming with respect to the adjacent lower dim entity */
  iMesh_getTagHandle(mi,dTAG_ADJ_REGION_FACE,&arf,&err,strlen(dTAG_ADJ_REGION_FACE));ICHKERRQ(mi,err);
  iMesh_getTagHandle(mi,dTAG_ADJ_FACE_EDGE,&afe,&err,strlen(dTAG_ADJ_FACE_EDGE));ICHKERRQ(mi,err);
  iMesh_getTagHandle(mi,dTAG_ORIENT_REGION_FACE,&orf,&err,strlen(dTAG_ORIENT_REGION_FACE));ICHKERRQ(mi,err);
  iMesh_getTagHandle(mi,dTAG_ORIENT_FACE_EDGE,&ofe,&err,strlen(dTAG_ORIENT_FACE_EDGE));ICHKERRQ(mi,err);
  /* Get full adjacencies */
  iMesh_getEHArrData(mi,mesh->r.v,mesh->r.s,arf,&mesh->arf.v,&mesh->arf.a,&mesh->arf.s,&err);ICHKERRQ(mi,err); /* region -> face */
  iMesh_getEHArrData(mi,mesh->f.v,mesh->f.s,afe,&mesh->afe.v,&mesh->afe.a,&mesh->afe.s,&err);ICHKERRQ(mi,err); /* face -> edge */
  iMesh_getEntArrAdj(mi,mesh->e.v,mesh->e.s,iBase_VERTEX,&mesh->aev.v,&mesh->aev.a,&mesh->aev.s,&off.v,&off.a,&off.s,&err);ICHKERRQ(mi,err); /* edge -> vertex */
  MeshListFree(off);      /* We don't use the offsets because we know there are always exactly two vertices per edge. */
  /* Get orientation of lower dimensional entities, we don't need vertex orientation */
  iMesh_getArrData(mi,mesh->r.v,mesh->r.s,orf,&mesh->orf.v,&mesh->orf.a,&mesh->orf.s,&err);ICHKERRQ(mi,err); /* region[face] */
  iMesh_getArrData(mi,mesh->f.v,mesh->f.s,ofe,&mesh->ofe.v,&mesh->ofe.a,&mesh->ofe.s,&err);ICHKERRQ(mi,err); /* face[edge] */
#endif

  /* View if requested */
  err = PetscOptionsHasName(((dObject)mesh)->prefix,"-dmesh_view",&flg);dCHK(err);
  if (flg) {
    dViewer viewer;
    err = PetscViewerASCIIGetStdout(((dObject)mesh)->comm,&viewer);dCHK(err);
    err = dMeshView(mesh,viewer);dCHK(err);
  }
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dMeshOrientFacets"
/*@
   dMeshOrientFacets -

@*/
dErr dMeshOrientFacets(dMesh m)
{
  //dErr err;

  dFunctionBegin;
  if (m || !m) dERROR(1,"not implemented");
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dMeshView"
/*@
   dMeshView -

@*/
dErr dMeshView(dMesh m,PetscViewer viewer)
{
  const char *type;
  iMesh_Instance mi;
  dBool iascii;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(m,dMESH_COOKIE,1);
  mi = m->mi;
  if (!viewer) {
    printf("Changing Viewer.");
    err = PetscViewerASCIIGetStdout(((PetscObject)m)->comm,&viewer);dCHK(err);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_COOKIE,2);
  PetscCheckSameComm(m,1,viewer,2);
  err = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);dCHK(err);
  if (iascii) {
    err = PetscObjectGetType((PetscObject)m,&type);dCHK(err);
    if (((PetscObject)m)->prefix) {
      err = PetscViewerASCIIPrintf(viewer,"dMesh object:(%s)\n",((PetscObject)m)->prefix);dCHK(err);
    } else {
      err = PetscViewerASCIIPrintf(viewer,"dMesh object:\n");dCHK(err);
    }
    err = PetscViewerASCIIPrintf(viewer,"Mesh type: %s\n",(type ? type : "not yet set"));dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Internal count by type: V=%d E=%d F=%d R=%d\n",m->v.s,m->e.s,m->f.s,m->r.s);dCHK(err);
    err = dMeshView_EntSet(m,m->root,viewer);dCHK(err);
    if (m->ops->view) {
      err = PetscViewerASCIIPushTab(viewer);dCHK(err);
      err = (*m->ops->view)(m,viewer);dCHK(err);
      err = PetscViewerASCIIPopTab(viewer);dCHK(err);
    }
  } else {
    if (m->ops->view) {
      err = (*m->ops->view)(m,viewer);dCHK(err);
    }
  }
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dMeshView_EntSet"
/*@
   dMeshView_EntSet -

@*/
static dErr dMeshView_EntSet(dMesh m,dMeshESH root,PetscViewer viewer)
{
  size_t valuesLen = 256;
  char values[256];
  iMesh_Instance mi = m->mi;
  char *tagname,*name,*z;
  int tagtype,tagsize,intdata;
  double dbldata;
  dMeshEH ehdata;
  MeshListTag tag=MLZ;
  MeshListData data=MLZ;
  MeshListESH esh=MLZ;
  dInt i,j,ntopo;
  dBool canprint;
  dErr err;

  dFunctionBegin;
  err = dMeshGetEntSetName(m,root,&name);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Entity Set %10p : %s\n",root,name);dCHK(err);
  err = PetscStrfree(name);dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  {
    for (i=iMesh_POINT; i<iMesh_ALL_TOPOLOGIES; i++) {
    iMesh_getNumOfTopo(mi,root,i,&ntopo,&err);ICHKERRQ(mi,err);
      if (ntopo) {
        err = PetscViewerASCIIPrintf(viewer,"%20s : %d\n",iMesh_TopologyName[i],ntopo);dCHK(err);
      }
    }
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);

  iMesh_getAllEntSetTags(mi,root,&tag.v,&tag.a,&tag.s,&err);ICHKERRQ(mi,err);
  err = PetscViewerASCIIPrintf(viewer,"Number of tags %d\n",tag.s);dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  {
    for (i=0; i<tag.s; i++) {
      err = dMeshGetTagName(m,tag.v[i],&tagname);dCHK(err);
      iMesh_getTagType(mi,tag.v[i],&tagtype,&err);ICHKERRQ(mi,err);
      iMesh_getTagSizeValues(mi,tag.v[i],&tagsize,&err);ICHKERRQ(mi,err);
      switch (tagtype) {        /* this needs a refactor */
        case iBase_INTEGER:
          iMesh_getEntSetIntData(mi,root,tag.v[i],&intdata,&err);ICHKERRQ(mi,err);
          err = PetscSNPrintf(values,valuesLen,"%d",intdata);dCHK(err);
          break;
        case iBase_DOUBLE:
          iMesh_getEntSetDblData(mi,root,tag.v[i],&dbldata,&err);ICHKERRQ(mi,err);
          err = PetscSNPrintf(values,valuesLen,"%f",dbldata);dCHK(err);
          break;
        case iBase_ENTITY_HANDLE:
          iMesh_getEntSetEHData(mi,root,tag.v[i],&ehdata,&err);ICHKERRQ(mi,err);
          err = PetscSNPrintf(values,valuesLen,"%p",ehdata);dCHK(err);
          break;
        case iBase_BYTES:
          iMesh_getEntSetData(mi,root,tag.v[i],&data.v,&data.a,&data.s,&err);ICHKERRQ(mi,err);
          canprint = PETSC_TRUE;
          for (j=0; j<data.s && data.v[j]; j++) {
            if (!isprint(data.v[i])) canprint = PETSC_FALSE;
          }
          if (canprint && false) {
            err = PetscSNPrintf(values,(size_t)data.s,"%s",data.v);dCHK(err); /* Just a copy, but ensures a NULL byte */
          } else {
            z = values;
            for (j=0; j<data.s && data.v[j] && (size_t)(z-values) < valuesLen-5; j++) {
              err = PetscSNPrintf(z,3,"%02uhhx ",data.v[j]);dCHK(err);
              z += 3;
              if (j%4 == 0) {
                *(z++) = ' ';
              }
              *(z++) = '\0';       /* Terminate the string */
            }
          }
          err = MeshListFree(data);dCHK(err);
          break;
        default: dERROR(1,"Invalid tag type, iMesh probably corrupt");
      }
      err = PetscViewerASCIIPrintf(viewer,"Tag: %30s : %20s [%3d] = %s\n",tagname,iBase_TagValueTypeName[tagtype],tagsize,values);dCHK(err);
      err = PetscFree(tagname);dCHK(err);
    }
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  err = MeshListFree(tag);dCHK(err);

  iMesh_getEntSets(mi,root,1,&esh.v,&esh.a,&esh.s,&err);ICHKERRQ(mi,err);
  err = PetscViewerASCIIPrintf(viewer,"Number of contained Entity Sets: %d\n",esh.s);dCHK(err);

  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  for (i=0; i<esh.s; i++) {
    err = PetscViewerASCIIPrintf(viewer,"Contained set %d/%d:\n",i+1,esh.s);dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = dMeshView_EntSet(m,esh.v[i],viewer);dCHK(err);
    err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  err = MeshListFree(esh);dCHK(err);

  iMesh_getChldn(mi,root,1,&esh.v,&esh.a,&esh.s,&err);ICHKERRQ(mi,err);
  err = PetscViewerASCIIPrintf(viewer,"Number of child Entity Sets: %d\n",esh.s);dCHK(err);

  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  for (i=0; i<esh.s; i++) {
    err = PetscViewerASCIIPrintf(viewer,"Child %d/%d:\n",i+1,esh.s);dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    err = dMeshView_EntSet(m,esh.v[i],viewer);dCHK(err);
    err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  }
  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  err = MeshListFree(esh);dCHK(err);
  dFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "dMeshGetEntSetName"
/*@
   dMeshGetEntSetName -

@*/
dErr dMeshGetEntSetName(dMesh m,dMeshESH set,char **str)
{
  MeshListData buf=MLZ;
  dMeshTag tag;
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(m,dMESH_COOKIE,1);
  dValidPointer(str,2);
  iMesh_getTagHandle(m->mi,dENT_SET_NAME,&tag,&err,strlen(dENT_SET_NAME));
  if (!err) {
    iMesh_getEntSetData(m->mi,set,tag,&buf.v,&buf.a,&buf.s,&err);
  }
  if (!err) {
    err = PetscStrallocpy(buf.v,str);dCHK(err);
    err = MeshListFree(buf);dCHK(err);
  } else if (err == iBase_TAG_NOT_FOUND) {
    err = PetscStrallocpy("NO_NAME",str);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dMeshGetInstance(dMesh m,iMesh_Instance *mi)
{

  dFunctionBegin;
  dValidHeader(m,dMESH_COOKIE,1);
  dValidPointer(mi,2);
  *mi = m->mi;
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dMeshDestroy"
/*@
   dMeshDestroy -

@*/
dErr dMeshDestroy(dMesh m)
{
  dErr err;

  dFunctionBegin;
  PetscValidHeaderSpecific(m,dMESH_COOKIE,1);
  if (m->ops->destroy) {
    err = (*m->ops->destroy)(m);dCHK(err);
  }
  MeshListFree(m->v); MeshListFree(m->e); MeshListFree(m->f); MeshListFree(m->r);
  MeshListFree(m->arf); MeshListFree(m->afe); MeshListFree(m->aev);
  MeshListFree(m->orf); MeshListFree(m->ofe);
  MeshListFree(m->x);
  iMesh_dtor(m->mi,&err);ICHKERRQ(m->mi,err);
  err = PetscStrfree(m->infile);dCHK(err);
  err = PetscStrfree(m->inoptions);dCHK(err);
  err = PetscHeaderDestroy(m);dCHK(err);
  dFunctionReturn(0);
}


/**
* Creates a \p dRule tag over all non-vertex topological entities in the mesh.  Also tags the root entity set with the
* given tag name and value equal to the base pointer for the rule storage.  This tag should be removed using
* dMeshDestroyRuleTag().
*
* @param mesh mesh
* @param esh entity set handle on which to add the tag, tags all non-vertex entities in \a esh
* @param jac dJacobi to use when generating the tags
* @param name unique identifier for the tag (mostly for debugging)
* @param degree polynomial degree which should be integrated exactly when the element has an affine map
* @param[out] inrtag tag
*
* @return err
*/
dErr dMeshCreateRuleTagIsotropic(dMesh mesh,dMeshESH esh,dJacobi jac,const char *name,dInt degree,dMeshTag *inrtag)
{
  static const int topologies[] = {0,iMesh_LINE_SEGMENT,iMesh_QUADRILATERAL,iMesh_HEXAHEDRON};
  static const int ntopologies[] = {0,1,1,1};
  const dInt rdeg[3] = {degree,degree,degree};
  iMesh_Instance mi = mesh->mi;
  dMeshTag rtag;
  dMeshESH root;
  MeshListEH ent=MLZ;
  void **base;
  dInt rulesize,needed,nents,maxnents,index;
  dRule *rules;
  size_t namelen;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  dValidPointer(inrtag,2);
  *inrtag = 0;
  err = PetscStrlen(name,&namelen);dCHK(err);
  iMesh_createTag(mi,name,sizeof(dRule),iBase_BYTES,&rtag,&err,(int)namelen);dICHK(mi,err);

  base = NULL; needed = 0; maxnents = 0;
  for (dInt pass=1; pass<=2; pass++) { /* determine preallocation, then allocate and fill */
    if (pass == 2) {
      err = dMalloc(needed*sizeof(*base),&base);dCHK(err);
      err = dMalloc(maxnents*sizeof(dMeshEH),&ent.v);dCHK(err);
      ent.a = maxnents; ent.s = 0;
      err = dMalloc(maxnents*sizeof(dRule),&rules);dCHK(err);
      index = 0;
    }
    for (int type=iBase_EDGE; type<iBase_ALL_TYPES; type++) {
      for (int topo=topologies[type]; topo<topologies[type]+ntopologies[type]; topo++) {
        if (pass == 1) {
          rulesize = 0;
          err = dJacobiGetRule(jac,topo,rdeg,NULL,NULL,&rulesize);dCHK(err); /* determine the amount of space needed for this topology */
          iMesh_getNumOfTopo(mi,esh,topo,&nents,&err);dICHK(mi,err);
          needed += nents*rulesize;
          maxnents = dMax(maxnents,nents);
        } else {
          iMesh_getEntities(mi,esh,type,topo,&ent.v,&ent.a,&ent.s,&err);dICHK(mi,err);
          for (dInt i=0; i<ent.s; i++) {
            err = dJacobiGetRule(jac,topo,rdeg,rules,&base[index],&index);dCHK(err);
          }
          iMesh_setArrData(mi,ent.v,ent.s,rtag,(char*)rules,ent.s*(int)sizeof(dRule),&err);dICHK(mi,err);
        }
      }
    }
  }
  err = dFree(ent.v);dCHK(err);
  err = dFree(rules);dCHK(err);
  iMesh_getRootSet(mi,&root,&err);dICHK(mi,err);
  if (sizeof(base) != sizeof(dRule)) dERROR(1,"Very strange, pointer sizes don't agree");
  iMesh_setEntSetData(mi,root,rtag,(char*)&base,sizeof(base),&err);dICHK(mi,err);
  *inrtag = rtag;
  dFunctionReturn(0);
}

dErr dMeshDestroyRuleTag(dMesh mesh,dMeshTag rtag)
{
  iMesh_Instance mi = mesh->mi;
  void *base;
  int bsize,balloc=sizeof(void*);
  dMeshESH root;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  iMesh_getRootSet(mi,&root,&err);dICHK(mi,err);
  iMesh_getEntSetData(mi,root,rtag,(char**)&base,&balloc,&bsize,&err);dICHK(mi,err);
  err = dFree(base);dCHK(err);
  dFunctionReturn(0);
}
