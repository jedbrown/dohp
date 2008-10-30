#ifndef _DOHPBASE_H
#define _DOHPBASE_H

#include "dohptype.h"

PETSC_EXTERN_CXX_BEGIN

#define INLINE static inline

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

INLINE void dGeomVecPlus(const dReal a[],const dReal b[],dReal c[]) {dInt i; for (i = 0; i<3; i++) c[i] = a[i]+b[i];}
INLINE void dGeomVecMinus(const dReal a[],const dReal b[],dReal c[]) {dInt i; for (i = 0; i<3; i++) c[i] = a[i]-b[i];}
INLINE dReal dGeomDotProd(const dReal a[],const dReal b[]) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
/* c = a \cross b */
INLINE void dGeomCrossProd(const dReal a[],const dReal b[],dReal c[])
{c[0]=a[1]*b[2]-a[2]*b[1]; c[1]=a[2]*b[0]-a[0]*b[2]; c[2]=a[0]*b[1]+a[1]*b[0];}
/* coordinates for vertices are interlaced in a, result goes in b */
INLINE void dGeomVecMeanI(dInt n,const dReal a[],dReal b[])
{dInt i,j; b[0]=b[1]=b[2]=0; for (i=0;i<n;i++) {for (j=0;j<3;j++) b[j]+=a[3*i+j]/n;}}
/* coordinates are interlaced in a, result in b */
INLINE void dGeomQuadFaceNormal(const dReal a[],dReal b[]) {
  dReal f[3],g[3],h[3];
  dGeomVecMinus(a+3,a,b); dGeomVecMinus(a+6,a+9,f); dGeomVecPlus(b,f,g); /* g = mean vector in x direction */
  dGeomVecMinus(a+9,a,b); dGeomVecMinus(a+6,a+3,f); dGeomVecPlus(b,f,h); /* h = mean vector in y direction */
  dGeomCrossProd(g,h,b);
}
INLINE dBool dGeomQuadParallel(const dReal a[],const dReal b[]) /* return true if Quad a \dot b > 0 */
{dReal f[3]; dGeomQuadFaceNormal(a,f); return (dGeomDotProd(b,f) > 0);}

INLINE void dGeomConvexComb_2_4(dReal x,dReal y,const dReal v[],const dInt p[],dReal f[])
{
  dInt i;
  for (i=0; i<3; i++) {
    f[i] = 0.25*((-1-x)*((-1-y)*v[p[0]*3+i] + (1+y)*v[p[3]*3+i]) + (1+x)*((-1-y)*v[p[1]*3+i] + (1+y)*v[p[2]*3+i]));
  }
}

/**
* These are macros so we can use them with ints and entity handles.
* 
*/
#define dGeomMatch2(a0,a1,b0,b1) (((a0)==(b0) && (a1)==(b1)) || ((a0)==(b1) && (a1)==(b0)))
#define dGeomMatchQuadLine(q,l)                                         \
  (dGeomMatch2((q)[0],(q)[1],(l)[0],(l)[1]) ? 0                         \
   : (dGeomMatch2((q)[1],(q)[2],(l)[0],(l)[1]) ? 1                      \
      : (dGeomMatch2((q)[2],(q)[3],(l)[0],(l)[1]) ? 2                   \
         : (dGeomMatch2((q)[3],(q)[0],(l)[0],(l)[1]) ? 3 : -1)))) 

/** 
* Finds the permutation of a Quad face so that it fits on a Hex, assumes that you know which number the face is.
* 
* @param rv vertex handles of the region (the Hex), length 8
* @param fv vertex handles of the face (the Quad), lenth 4
* @param fnum face number in canonical face ordering
* @param orient permutation, range {0..7}, error thrown if they don't fit
* 
* @return 
*/
INLINE dErr dGeomOrientFindPerm_HexQuad(const dMeshEH rv[],const dMeshEH fv[],int fnum,dInt *orient)
{
  static const dInt perm[8][4] = {{0,1,2,3},{1,2,3,0},{2,3,0,1},{3,0,1,2},
                                  {0,3,2,1},{3,2,1,0},{2,1,0,3},{1,0,3,2}};
  static const dInt permorient[8] = {0,1,2,3,4,5,6,7};
  dInt i;

  dFunctionBegin;
#if 0
  printf("# rv:");              /* The vertices of face fnum on this region in canonical order. */
  for (i=0; i<4; i++) { printf(" %ld",(long)rv[dMeshConnectHexQuad[fnum][i]]); }
  printf("\n# fv:");
  /* The order of vertices on the face.  We want to find a permutation of these
  vertices to matches the canonical order. */
  for (i=0; i<4; i++) { printf(" %ld",(long)fv[i]); }
  printf("\n");
#endif
  for (i=0; i<8; i++) { /* loop over permutations */
    if (fv[perm[i][0]] == rv[dMeshConnectHexQuad[fnum][0]] && fv[perm[i][1]] == rv[dMeshConnectHexQuad[fnum][1]]) {
      /* we have found a match, as an extra check we can check that the other vertices match */
      if (fv[perm[i][2]] != rv[dMeshConnectHexQuad[fnum][2]] || fv[perm[i][3]] != rv[dMeshConnectHexQuad[fnum][3]]) {
        dERROR(1,"Faces cannot be matched, but part matches, perhaps adjacencies are corrupt");
      }
      *orient = permorient[i];
      dFunctionReturn(0);
    }
  }
  dERROR(1,"Faces cannot be matched");
  dFunctionReturn(0);
}

INLINE dErr dGeomOrientFindPerm_QuadLine(const dMeshEH fv[],const dMeshEH ev[],int en,dInt *orient)
{
  static const dInt perm[2][2] = {{0,1},{1,0}};
  static const dInt permorient[4] = {0,1};
  dInt i;

  dFunctionBegin;
#if 0
  printf("# fv:");              /* The vertices of edge en on this face in canonical order. */
  for (i=0; i<2; i++) { printf(" %ld",(long)fv[dMeshConnectQuadLine[en][i]]); }
  printf("\n# ev:");
  /* The order of vertices on the edge.  We want to find a permutation of these
  vertices to matches the canonical order. */
  for (i=0; i<2; i++) { printf(" %ld",(long)ev[i]); }
  printf("\n");
#endif
  for (i=0; i<2; i++) { /* loop over permutations */
    if (ev[perm[i][0]] == fv[dMeshConnectQuadLine[en][0]] && ev[perm[i][1]] == fv[dMeshConnectQuadLine[en][1]]) {
      *orient = permorient[i];
      dFunctionReturn(0);
    }
  }
  dERROR(1,"Edges cannot be matched.");
  dFunctionReturn(0);
}



PETSC_EXTERN_CXX_END
#endif  /* _DOHPBASE_H */
