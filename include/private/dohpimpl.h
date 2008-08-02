#ifndef _DOHPIMPL_H
#define _DOHPIMPL_H

#include "dohp.h"

struct _DohpMeshOps {
  PetscErrorCode(*orientfacets)(DohpMesh);
};

struct _p_DohpMesh {
  PETSCHEADER(struct _DohpMeshOps);
  iMesh_Instance m;
  iBase_EntitySetHandle root;
  MeshListEH v,e,f,r;
  MeshListReal x;
};

struct _DohpDMOps {
};

struct _p_DohpDM {
  PETSCHEADER(struct _DohpDMOps);

  /* Distributed object */
  DohpMesh mesh;

  /* Defined on the locally owned portion of the mesh, can be reconstructed from information stored in the mesh. */
  DohpRule_Hex *ruleR;
  DohpRule_Quad *ruleF;
  DohpEMap_Hex *emapR;
  DohpEMap_Quad *emapF;
};

struct _EFSOps {
  PetscErrorCode (*deriv)(void *,const PetscScalar *,PetscScalar *); /* evaluate derivatives at quadrature points */
  PetscErrorCode (*derivt)(void *,const PetscScalar *,PetscScalar *); /* weak derivatives of test functions */
  PetscErrorCode (*getnodes)(void *,PetscReal **);                   /* coordinates of the nodal basis */
  PetscErrorCode (*getqweights)(void *,PetscReal **);                /* weights at quadrature points  */
  PetscErrorCode (*getqnodes)(void *,PetscReal *);                   /* nodes */
  PetscErrorCode (*loctoint)(void *,const PetscScalar *,PetscScalar *);
  PetscErrorCode (*inttoloc)(void *,const PetscScalar *,PetscScalar *);
  PetscErrorCode (*facettoelem)(void *,void *,const PetscScalar *,PetscScalar *);
  PetscErrorCode (*elemtofacet)(void *,void *,const PetscScalar *,PetscScalar *);
};


// Held once for every Function Space with support on this element.
typedef struct {
  // Normally a constant pointer to the operations for a basis type, but could point to unrolled versions.
  struct _EFSOps ops;
  // Private storage to define the basis operations.  For Hex elements, this
  // would be three pointers to the line contexts in the tensor product.
  // Used by EFSOps for operations in the star space. Must agree with quadrature order.
  void *elem;
  // Private storage, defines how to project to/from facets.
  // An example implementation would have pointers to projection matrices.
  void *facet;
} DohpEFS;

// For every Field, we need only the number of degrees of freedom and the local offset.
typedef struct {
  PetscInt localoffset;
  // PetscInt ndof;
} EField;

typedef struct {
  PetscInt *localoffset;
  PetscInt ndof;
} MField;

struct _DohpMFSOps {};

struct _p_DohpMFS {
  PETSCHEADER(struct _MFSOps);
  DohpEFS *efs;
};

// The facets of a line are 2 vertices.
// The projection matrix is the 1x1 identity so vertices are indexed straight out of the local vector.
typedef struct {
  PetscInt localoffset[2];
} ElemFacets_Line;

// The facets of a quad are 4 lines.
// F^* -> R^*: c[i] = SUM(j=0:size[fnum]) projrf[fnum][i*size[fnum]+j] * f[j]
// R^* -> F^*: f[i] = SUM(j=0:csize)      projfr[fnum][i*csize+j] * c[j]
typedef struct {
  PetscInt     size[4];         // Encode orientation and offset in the high-order bits of this value
  PetscScalar *projef[4];       // Element to facet projection
  PetscScalar *projfe[4];       // Facet to element projection
  PetscInt     facetoffset[4];  // Edge degrees of freedom are indexed out of the facet vector.
} ElemFacets_Quad;

// The facets of a hex are 6 quads.
typedef struct {
  PetscInt     size[6];         // How do we encode both x- and y-orientation in the high bits?
  PetscScalar *projef[6];       // Element to facet projection
  PetscScalar *projfe[6];       // Facet to element projection
  PetscInt     facetoffset[6];  // Face degrees of freedom are indexed out of the facet vector.
} ElemFacets_Hex;


/* The Mesh Quotient map.  A quadrature rule and element coordinate mapping
 defines a quotient map on the Hilbert space $H^1(\Omega)$ which induces a
 finite topology.  The same quotient map can be used for many different function
 spaces $X \subset H^1$. */

// Since element quadrature and mapping contexts may have different sizes in
// different elements, we cannot simply use the PETSc data structures to handle
// local to global communication.  For now, we can just use flat arrays of
// pointers on the local process, but this will require some thought when we
// want to do dynamic load balancing.
typedef struct _p_MQuot *MQuot;
struct _p_MQuot {
  PETSCHEADER(struct _MQuotOps);
  void **equad; // locally owned element quadrature context
  void **emap; // locally owned element map context
};

#endif /* _DOHPIMPL_H */
