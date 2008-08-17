#include "dohp.h"

const char *iBase_ErrorString[] = {
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
