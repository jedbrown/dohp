#include "private/dohpimpl.h"
#include <ctype.h>              /* needed for isprint() */

static PetscErrorCode DohpMeshView_EntSet(DohpMesh m,DohpESH root,PetscViewer viewer);

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

const char *const iMesh_TopologyName[] = {
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

const char *const iBase_TagValueTypeName[] = {
  "iBase_INTEGER",
  "iBase_DOUBLE",
  "iBase_ENTITY_HANDLE",
  "iBase_BYTES"
};

#undef __FUNCT__
#define __FUNCT__ "MeshListIntView"
PetscErrorCode MeshListIntView(MeshListInt *m, const char *name) {
  CHKERRQ(PetscPrintf(PETSC_COMM_SELF,"# %s [%d]\n", name, m->s));
  CHKERRQ(PetscIntView(m->s,m->v,PETSC_VIEWER_STDOUT_SELF));
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpOrientFindPerm_HexQuad"
PetscErrorCode DohpOrientFindPerm_HexQuad(const iBase_EntityHandle *rv, const iBase_EntityHandle *fv,
                                          int fnum, DohpOrient *orient)
{
  PetscInt i;
  static const PetscInt perm[8][4] = {{0,1,2,3},{1,2,3,0},{2,3,0,1},{3,0,1,2},
                                      {0,3,2,1},{3,2,1,0},{2,1,0,3},{1,0,3,2}};
  static const DohpOrient permorient[8] = {0,1,2,3,4,5,6,7};

  PetscFunctionBegin;
#if 0
  printf("# rv:");              /* The vertices of face fnum on this region in canonical order. */
  for (i=0; i<4; i++) { printf(" %ld",(long)rv[DohpHexQuad[fnum][i]]); }
  printf("\n# fv:");
  /* The order of vertices on the face.  We want to find a permutation of these
  vertices to matches the canonical order. */
  for (i=0; i<4; i++) { printf(" %ld",(long)fv[i]); }
  printf("\n");
#endif
  for (i=0; i<8; i++) { /* loop over permutations */
    if (fv[perm[i][0]] == rv[DohpHexQuad[fnum][0]] && fv[perm[i][1]] == rv[DohpHexQuad[fnum][1]]) {
      /* we have found a match, as an extra check we can check that the other vertices match */
      if (fv[perm[i][2]] != rv[DohpHexQuad[fnum][2]] || fv[perm[i][3]] != rv[DohpHexQuad[fnum][3]]) {
        SETERRQ(1,"Faces cannot be matched.");
      }
      *orient = permorient[i];
      break;
    }
  }
  if (i==8) SETERRQ(1,"Faces cannot be matched.");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpOrientFindPerm_QuadLine"
PetscErrorCode DohpOrientFindPerm_QuadLine(const iBase_EntityHandle *fv, const iBase_EntityHandle *ev,
                                          int en, DohpOrient *orient)
{
  static const PetscInt perm[2][2] = {{0,1},{1,0}};
  static const DohpOrient permorient[4] = {0,1,2,3};
  PetscInt i;

  PetscFunctionBegin;
#if 0
  printf("# fv:");              /* The vertices of edge en on this face in canonical order. */
  for (i=0; i<2; i++) { printf(" %ld",(long)fv[DohpQuadLine[en][i]]); }
  printf("\n# ev:");
  /* The order of vertices on the edge.  We want to find a permutation of these
  vertices to matches the canonical order. */
  for (i=0; i<2; i++) { printf(" %ld",(long)ev[i]); }
  printf("\n");
#endif
  for (i=0; i<2; i++) { /* loop over permutations */
    if (ev[perm[i][0]] == fv[DohpQuadLine[en][0]] && ev[perm[i][1]] == fv[DohpQuadLine[en][1]]) {
      *orient = permorient[i];
      break;
    }
  }
  if (i==4) SETERRQ(1,"Edges cannot be matched.");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpOrientLoopBounds_Quad"
PetscErrorCode DohpOrientLoopBounds_Quad(DohpOrient orient, const PetscInt *size, DohpLoopBounds *l)
{
  const PetscInt ox=size[0], oy=size[1];

  PetscFunctionBegin;
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
      SETERRQ(1,"Orientation not supported.");
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpOrientLoopBounds_Line"
PetscErrorCode DohpOrientLoopBounds_Line(DohpOrient orient, const PetscInt *size, DohpLoopBounds *l)
{

  PetscFunctionBegin;
  switch (orient) {
    case 0: l->start = 0;         l->stride = 1;  l->end = size[0]; break;
    case 1: l->start = size[0]-1; l->stride = -1; l->end = -1;       break;
    default: SETERRQ(1,"Orientation not supported.");
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpLoopBounds_Hex"
/* On each face, we need a loop which traverses the face (indicated by
* DohpHexQuad[][]) in the positive order.  The ordering of degrees of freedom on the
* Hex is [i][j][k] (C-style ordering). */
PetscErrorCode DohpLoopBounds_Hex(const PetscInt *size, PetscInt face, DohpLoopBounds *l)
{
  const PetscInt ox=size[0], oy=size[1], oz=size[2];
  PetscFunctionBegin;
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
      SETERRQ(1,"Face number not recognized.");
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpLoopBounds_Quad"
PetscErrorCode DohpLoopBounds_Quad(const PetscInt *size, PetscInt edge, DohpLoopBounds *l)
{
  const PetscInt ox=size[0], oy=size[1];

  PetscFunctionBegin;
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
      SETERRQ(1,"Edge number not recognized.");
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EFSFacetToElem_HexQuad_Conforming"
/* Maps facet degrees of freedom to element degrees of freedom, adding
* contributions.  This function is actually an optimization for conforming
* elements since it does not need to do interpolation. */
PetscErrorCode EFSFacetToElem_HexQuad_Conforming(PetscInt dof,const PetscInt rsize[],const PetscInt fsize[],PetscInt fnum,DohpOrient forient,const PetscScalar fvals[],PetscScalar rvals[])
{
  PetscInt ri,rj,fi,fj,k;
  DohpLoopBounds rl[2],fl[2];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DohpLoopBounds_Hex(rsize,fnum,rl);CHKERRQ(ierr);
  ierr = DohpOrientLoopBounds_Quad(forient,fsize,fl);CHKERRQ(ierr);
  for (ri=rl[0].start,fi=fl[0].start; ri!=rl[0].end && fi!=fl[0].end; ri+=rl[0].stride,fi+=fl[0].stride) {
    for (rj=rl[1].start,fj=fl[1].start; rj!=rl[1].end && fj!=fl[1].end; rj+=rl[1].stride,fj+=fl[1].stride) {
      for (k=0; k<dof; k++) {
        rvals[(ri+rj)*dof+k] += fvals[(fi+fj)*dof+k];
      }
    }
    if (!(rj==rl[1].end && fj==fl[1].end)) {
      SETERRQ(1,"Inner loop bounds do not agree.  Is this relation conforming?");
    }
  }
  if (!(ri==rl[0].end && fi==fl[0].end)) {
    SETERRQ(1,"Outer loop bounds do not agree.  Is this relation conforming?");
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EFSFacetToElem_QuadLine_Conforming"
PetscErrorCode EFSFacetToElem_QuadLine_Conforming(PetscInt dof,const PetscInt fsize[],const PetscInt esize[],PetscInt en,DohpOrient eorient,const PetscScalar evals[],PetscScalar fvals[])
{
  PetscInt fi,ei,j;
  DohpLoopBounds fl,el;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DohpLoopBounds_Quad(fsize,en,&fl);CHKERRQ(ierr);
  ierr = DohpOrientLoopBounds_Line(eorient,esize,&el);CHKERRQ(ierr);
  for (fi=fl.start,ei=el.start; fi!=fl.end && ei!=el.end; fi+=fl.stride,ei+=el.stride) {
    for (j=0; j<dof; j++) {
      fvals[fi*dof+j] += evals[ei*dof+j];
    }
  }
  if (!(fi==fl.end && ei==el.end)) {
    SETERRQ(1,"Loop bounds do not agree.  Is this relation conforming?");
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpMeshGetTagName"
/*@
   DohpMeshGetTagName -

   This function allocates memory for the tag.  It should be freed with PetscFree()

@*/
PetscErrorCode DohpMeshGetTagName(DohpMesh m,DohpTag tag,char **name)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(m,DOHP_MESH_COOKIE,1);
  PetscValidPointer(name,2);
  ierr = PetscMalloc(DOHP_NAME_LEN,name);CHKERRQ(ierr);
  iMesh_getTagName(m->mi,tag,*name,&ierr,DOHP_NAME_LEN);ICHKERRQ(m->mi,ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DohpMeshCreate"
PetscErrorCode DohpMeshCreate(MPI_Comm comm,DohpMesh *inm)
{
  static const struct _DohpMeshOps dfltOps = { .orientfacets = DohpMeshOrientFacets, .view = 0, .destroy = 0 };
  DohpMesh m;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidPointer(inm,2);
#if !defined(PETSC_USE_DYNAMIC_LIBRARIES)
  ierr = DohpQuotientInitializePackage(PETSC_NULL);CHKERRQ(ierr);
#endif
  *inm = 0;
  ierr = PetscHeaderCreate(m,_p_DohpMesh,struct _DohpMeshOps,DOHP_MESH_COOKIE,0,"DohpMesh",comm,DohpMeshDestroy,DohpMeshView);CHKERRQ(ierr);
  ierr = PetscObjectChangeTypeName((PetscObject)m,"iMesh");CHKERRQ(ierr);
  iMesh_newMesh("",&m->mi,&ierr,0);ICHKERRQ(m->mi,ierr);
  ierr = PetscMemcpy(m->ops,&dfltOps,sizeof(dfltOps));CHKERRQ(ierr);
  *inm = m;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpMeshLoad"
PetscErrorCode DohpMeshLoad(DohpMesh m,const char fname[],const char opt[])
{
  iBase_TagHandle arf,afe,orf,ofe;
  MeshListInt off=MLZ;
  iMesh_Instance mi;
  size_t fnamelen,optlen;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscStrlen(fname,&fnamelen);CHKERRQ(ierr);
  ierr = PetscStrlen(opt,&optlen);CHKERRQ(ierr);
  mi = m->mi;
  iMesh_load(mi,0,fname,opt,&ierr,(int)fnamelen,(int)optlen);ICHKERRQ(mi,ierr);
  iMesh_getRootSet(mi,&m->root,&ierr);ICHKERRQ(mi,ierr);
  /* Get all entities of each type. */
  iMesh_getEntities(mi,m->root,iBase_REGION,iMesh_ALL_TOPOLOGIES,&m->r.v,&m->r.a,&m->r.s,&ierr);ICHKERRQ(mi,ierr);
  iMesh_getEntities(mi,m->root,iBase_FACE,iMesh_ALL_TOPOLOGIES,&m->f.v,&m->f.a,&m->f.s,&ierr);ICHKERRQ(mi,ierr);
  iMesh_getEntities(mi,m->root,iBase_EDGE,iMesh_ALL_TOPOLOGIES,&m->e.v,&m->e.a,&m->e.s,&ierr);ICHKERRQ(mi,ierr);
  iMesh_getEntities(mi,m->root,iBase_VERTEX,iMesh_ALL_TOPOLOGIES,&m->v.v,&m->v.a,&m->v.s,&ierr);ICHKERRQ(mi,ierr);
  /* Get tags for custom adjacencies, needed since our meshes are nonconforming with respect to the adjacent lower dim entity */
  iMesh_getTagHandle(mi,DOHP_TAG_ADJ_REGION_FACE,&arf,&ierr,strlen(DOHP_TAG_ADJ_REGION_FACE));ICHKERRQ(mi,ierr);
  iMesh_getTagHandle(mi,DOHP_TAG_ADJ_FACE_EDGE,&afe,&ierr,strlen(DOHP_TAG_ADJ_FACE_EDGE));ICHKERRQ(mi,ierr);
  iMesh_getTagHandle(mi,DOHP_TAG_ORIENT_REGION_FACE,&orf,&ierr,strlen(DOHP_TAG_ORIENT_REGION_FACE));ICHKERRQ(mi,ierr);
  iMesh_getTagHandle(mi,DOHP_TAG_ORIENT_FACE_EDGE,&ofe,&ierr,strlen(DOHP_TAG_ORIENT_FACE_EDGE));ICHKERRQ(mi,ierr);
  /* Get full adjacencies */
  iMesh_getEHArrData(mi,m->r.v,m->r.s,arf,&m->arf.v,&m->arf.a,&m->arf.s,&ierr);ICHKERRQ(mi,ierr); /* region -> face */
  iMesh_getEHArrData(mi,m->f.v,m->f.s,afe,&m->afe.v,&m->afe.a,&m->afe.s,&ierr);ICHKERRQ(mi,ierr); /* face -> edge */
  iMesh_getEntArrAdj(mi,m->e.v,m->e.s,iBase_VERTEX,&m->aev.v,&m->aev.a,&m->aev.s,&off.v,&off.a,&off.s,&ierr);ICHKERRQ(mi,ierr); /* edge -> vertex */
  MeshListFree(off);      /* We don't use the offsets because we know there are always exactly two vertices per edge. */
  /* Get orientation of lower dimensional entities, we don't need vertex orientation */
  iMesh_getArrData(mi,m->r.v,m->r.s,orf,&m->orf.v,&m->orf.a,&m->orf.s,&ierr);ICHKERRQ(mi,ierr); /* region[face] */
  iMesh_getArrData(mi,m->f.v,m->f.s,ofe,&m->ofe.v,&m->ofe.a,&m->ofe.s,&ierr);ICHKERRQ(mi,ierr); /* face[edge] */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpMeshOrientFacets"
/*@
   DohpMeshOrientFacets - 

@*/
PetscErrorCode DohpMeshOrientFacets(DohpMesh m)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (ierr || m || !m) SETERRQ(1,"not implemented");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpMeshView"
/*@
   DohpMeshView - 

@*/
PetscErrorCode DohpMeshView(DohpMesh m,PetscViewer viewer)
{
  const char *type;
  iMesh_Instance mi;
  PetscTruth iascii;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(m,DOHP_MESH_COOKIE,1);
  mi = m->mi;
  if (!viewer) {
    printf("Changing Viewer.");
    ierr = PetscViewerASCIIGetStdout(((PetscObject)m)->comm,&viewer);CHKERRQ(ierr);
  }
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_COOKIE,2);
  PetscCheckSameComm(m,1,viewer,2);
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscObjectGetType((PetscObject)m,&type);CHKERRQ(ierr);
    if (((PetscObject)m)->prefix) {
      ierr = PetscViewerASCIIPrintf(viewer,"DohpMesh object:(%s)\n",((PetscObject)m)->prefix);CHKERRQ(ierr);
    } else {
      ierr = PetscViewerASCIIPrintf(viewer,"DohpMesh object:\n");CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"Mesh type: %s\n",(type ? type : "not yet set"));CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"Internal count by type: V=%d E=%d F=%d R=%d\n",m->v.s,m->e.s,m->f.s,m->r.s);CHKERRQ(ierr);
    ierr = DohpMeshView_EntSet(m,m->root,viewer);CHKERRQ(ierr);
    if (m->ops->view) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = (*m->ops->view)(m,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
  } else {
    if (m->ops->view) {
      ierr = (*m->ops->view)(m,viewer);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpMeshView_EntSet"
/*@
   DohpMeshView_EntSet - 

@*/
PetscErrorCode DohpMeshView_EntSet(DohpMesh m,DohpESH root,PetscViewer viewer)
{
  size_t valuesLen = 256;
  char values[256];
  iMesh_Instance mi = m->mi;
  char *tagname,*name,*z;
  int tagtype,tagsize,intdata;
  double dbldata;
  DohpEH ehdata;
  MeshListTag tag=MLZ;
  MeshListData data=MLZ;
  MeshListESH esh=MLZ;
  PetscInt i,j,ntopo;
  PetscTruth canprint;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DohpMeshGetEntSetName(m,root,&name);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Entity Set %10p : %s\n",root,name);CHKERRQ(ierr);
  if (name) {
    ierr = PetscFree(name);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  {
    for (i=iMesh_POINT; i<iMesh_ALL_TOPOLOGIES; i++) {
    iMesh_getNumOfTopo(mi,root,i,&ntopo,&ierr);ICHKERRQ(mi,ierr);
      if (ntopo) {
        ierr = PetscViewerASCIIPrintf(viewer,"%20s : %d\n",iMesh_TopologyName[i],ntopo);CHKERRQ(ierr);
      }
    }
  }
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);

  iMesh_getAllEntSetTags(mi,root,&tag.v,&tag.a,&tag.s,&ierr);ICHKERRQ(mi,ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Number of tags %d\n",tag.s);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  {
    for (i=0; i<tag.s; i++) {
      ierr = DohpMeshGetTagName(m,tag.v[i],&tagname);CHKERRQ(ierr);
      iMesh_getTagType(mi,tag.v[i],&tagtype,&ierr);ICHKERRQ(mi,ierr);
      iMesh_getTagSizeValues(mi,tag.v[i],&tagsize,&ierr);ICHKERRQ(mi,ierr);
      switch (tagtype) {        /* this needs a refactor */
        case iBase_INTEGER:
          iMesh_getEntSetIntData(mi,root,tag.v[i],&intdata,&ierr);ICHKERRQ(mi,ierr);
          ierr = PetscSNPrintf(values,valuesLen,"%d",intdata);CHKERRQ(ierr);
          break;
        case iBase_DOUBLE:
          iMesh_getEntSetDblData(mi,root,tag.v[i],&dbldata,&ierr);ICHKERRQ(mi,ierr);
          ierr = PetscSNPrintf(values,valuesLen,"%f",dbldata);CHKERRQ(ierr);
          break;
        case iBase_ENTITY_HANDLE:
          iMesh_getEntSetEHData(mi,root,tag.v[i],&ehdata,&ierr);ICHKERRQ(mi,ierr);
          ierr = PetscSNPrintf(values,valuesLen,"%p",ehdata);CHKERRQ(ierr);
          break;
        case iBase_BYTES:
          iMesh_getEntSetData(mi,root,tag.v[i],&data.v,&data.a,&data.s,&ierr);ICHKERRQ(mi,ierr);
          canprint = PETSC_TRUE;
          for (j=0; j<data.s && data.v[j]; j++) {
            if (!isprint(data.v[i])) canprint = PETSC_FALSE;
          }
          if (canprint) {
            ierr = PetscSNPrintf(values,(size_t)data.s,"%s",data.v);CHKERRQ(ierr); /* Just a copy, but ensures a NULL byte */
          } else {
            z = values;
            for (j=0; j<data.s && data.v[j] && (size_t)(z-values) < valuesLen-5; j++) {
              ierr = PetscSNPrintf(z,3,"%02x ",data.v[j]);CHKERRQ(ierr);
              z += 3;
              if (j%4 == 0) {
                *(z++) = ' ';
              }
              *(z++) = '\0';       /* Terminate the string */
            }
          }
          ierr = MeshListFree(data);CHKERRQ(ierr);
          break;
        default: SETERRQ(1,"Invalid tag type, iMesh probably corrupt");
      }
      ierr = PetscViewerASCIIPrintf(viewer,"Tag: %30s : %20s [%3d] = %s\n",tagname,iBase_TagValueTypeName[tagtype],tagsize,values);CHKERRQ(ierr);
      ierr = PetscFree(tagname);CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  ierr = MeshListFree(tag);CHKERRQ(ierr);

  iMesh_getEntSets(mi,root,1,&esh.v,&esh.a,&esh.s,&ierr);ICHKERRQ(mi,ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Number of contained Entity Sets: %d\n",esh.s);CHKERRQ(ierr);

  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  for (i=0; i<esh.s; i++) {
    ierr = PetscViewerASCIIPrintf(viewer,"Contained set %d/%d:\n",i+1,esh.s);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = DohpMeshView_EntSet(m,esh.v[i],viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  ierr = MeshListFree(esh);CHKERRQ(ierr);

  iMesh_getChldn(mi,root,1,&esh.v,&esh.a,&esh.s,&ierr);ICHKERRQ(mi,ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Number of child Entity Sets: %d\n",esh.s);CHKERRQ(ierr);

  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  for (i=0; i<esh.s; i++) {
    ierr = PetscViewerASCIIPrintf(viewer,"Child %d/%d:\n",i+1,esh.s);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = DohpMeshView_EntSet(m,esh.v[i],viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  ierr = MeshListFree(esh);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DohpMeshGetEntSetName"
/*@
   DohpMeshGetEntSetName - 

@*/
PetscErrorCode DohpMeshGetEntSetName(DohpMesh m,DohpESH set,char **str)
{
  MeshListData buf=MLZ;
  DohpTag tag;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(m,DOHP_MESH_COOKIE,1);
  PetscValidPointer(str,2);
  iMesh_getTagHandle(m->mi,DOHP_ENT_SET_NAME,&tag,&ierr,strlen(DOHP_ENT_SET_NAME));ICHKERRQ(m->mi,ierr);
  iMesh_getEntSetData(m->mi,set,tag,&buf.v,&buf.a,&buf.s,&ierr);
  if (!ierr) {
    ierr = PetscStrallocpy(buf.v,str);CHKERRQ(ierr);
    ierr = MeshListFree(buf);CHKERRQ(ierr);
  } else if (ierr == iBase_TAG_NOT_FOUND) {
    ierr = PetscStrallocpy("NO_NAME",str);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpMeshDestroy"
/*@
   DohpMeshDestroy - 

@*/
PetscErrorCode DohpMeshDestroy(DohpMesh m)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(m,DOHP_MESH_COOKIE,1);
  if (m->ops->destroy) {
    ierr = (*m->ops->destroy)(m);CHKERRQ(ierr);
  }
  MeshListFree(m->v); MeshListFree(m->e); MeshListFree(m->f); MeshListFree(m->r);
  MeshListFree(m->arf); MeshListFree(m->afe); MeshListFree(m->aev);
  MeshListFree(m->orf); MeshListFree(m->ofe);
  MeshListFree(m->x);
  iMesh_dtor(m->mi,&ierr);ICHKERRQ(m->mi,ierr);
  ierr = PetscHeaderDestroy(m);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpMeshRegisterAll"
/*@
   DohpMeshRegisterAll - 

@*/
PetscErrorCode DohpMeshRegisterAll(const char path[])
{
  static PetscTruth registered = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (registered) PetscFunctionReturn(0);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"DohpMeshRegisterAll: %s (nothing to do)\n",path);CHKERRQ(ierr);
  registered = PETSC_TRUE;
  PetscFunctionReturn(0);
}
