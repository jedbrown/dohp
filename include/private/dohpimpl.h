#ifndef _DOHPIMPL_H
#define _DOHPIMPL_H

#include "dohp.h"
#include "src/dm/dmimpl.h"

PetscCookie DOHP_DM_COOKIE;
PetscLogEvent DOHP_MatMult, DOHP_FunctionEval, DOHP_JacobianEval;

struct _DohpMeshOps {
  PetscErrorCode(*orientfacets)(DohpMesh);
  PetscErrorCode(*view)(DohpMesh,PetscViewer);
  PetscErrorCode(*destroy)(DohpMesh);
};

struct _p_DohpMesh {
  PETSCHEADER(struct _DohpMeshOps);
  iMesh_Instance mi;
  iBase_EntitySetHandle root;
  MeshListEH v,e,f,r;           /* vertices, edges, faces, vertices */
  MeshListEH arf,afe,aev;       /* adjacencies region -> face -> edge -> vertex */
  MeshListData orf,ofe;
  MeshListReal x;
};

typedef struct {
  PetscErrorCode (*dofspernode)(/* what goes here? */);
  /* number of fields, number of nodes, normal vector, boundary values, element nodal values (in/out), output vector (out) */
  PetscErrorCode (*apply)(PetscInt,PetscTruth,const PetscScalar*,const PetscScalar*,PetscScalar*,PetscScalar*);
  PetscErrorCode (*unapply)(PetscInt,PetscInt,const PetscScalar*,const PetscScalar*,PetscScalar*,PetscScalar*);
} DohpBCOps;

typedef unsigned char DohpUInt8;

struct _DohpMFSOps {};

struct _DohpDMOps {
  DMOps(DohpDM);
};

struct _p_DohpDM {
  PETSCHEADER(struct _DohpDMOps);
  DMHEADER
  /* Distributed object */
  DohpMesh mesh;

  /* Defined on the locally owned portion of the mesh, can be reconstructed from information stored in the mesh. */
  DohpRule_Hex  *ruleR;
  DohpRule_Quad *ruleF;
  DohpEMap_Hex  *emapR;
  DohpEMap_Quad *emapF;

  /* Approximation space */
  PetscInt   nmfs;              /* Number of function spaces in this DM */
  DohpMFS   *mfs;               /* handles of the spaces */
  PetscInt  *nfields;           /* [nmfs], number of fields on each function space */
  char ***fieldnames;           /* names of the fields in each function space */
  PetscInt   nelems;            /* total number of elements */
  PetscInt  *ind;               /* [nelems], starting index of interior dofs in local vector */
  PetscInt  *find;              /* index of face dofs in local vector */
  PetscInt  *foff;              /* [nelems], offset into find for each element */
  void     **edata;             /* private data (Jacobian storage) for each element */

  /* Blocks: during preconditioning, the */
  /* Preconditioning: this information is used for matrix-free application of blocks and assembly of low-order blocks */
  PetscInt **fieldoff;          /* [nfields][nelems], offset to add to ind[r] to get interior dofs in local vector */
  PetscInt **sind;              /* [nfields][nelems], starting index of interior dofs in split local vectors */
  PetscInt  *sfind;             /* index of face dofs in split local vector */
  PetscInt  *sfoff;             /* [nelems], offset into sfind for each element */
};

struct _DohpBlockOps {
};

struct _p_DohpBlock {
  PETSCHEADER(struct _DohpBlockOps);
  DohpDM dm;
  Sliced sliced;
  VecScatter unifiedtosplit;    /* Extract the relevant fields from a unified vector */
  VecScatter localtoglobal;     /* Scatter the local vector to the global split vector */

  PetscInt lmfs,rmfs;           /* Index of left and right MFS in  */
  PetscInt nlf,nrf;             /* Number of left fields, number of right fields */
  PetscInt *lfields,*rfields;   /* Indices of left and right fields in dm->[] */

  PetscInt nelems;              /* Total number of elements active in this block */
  PetscInt *ind;                /* [nelems] starting index of interior dofs in local vector */
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
struct _p_DohpEFS {
  // Normally a constant pointer to the operations for a basis type, but could point to unrolled versions.
  struct _EFSOps *ops;
  // Private storage to define the basis operations.  For Hex elements, this
  // would be three pointers to the line contexts in the tensor product.
  // Used by EFSOps for operations in the star space. Must agree with quadrature order.
  void *elem;
  // Private storage, defines how to project to/from facets.
  // An example implementation would have pointers to projection matrices.
  void *facet;
};

struct _p_DohpMFS {
  PETSCHEADER(struct _DohpMFSOps);
  DohpMesh  mesh;
  char     *name;               /* name of the domain, used to find the corresponding entity set in the mesh */
  char     *sizetagname;        /* a tag with this name is stored in the mesh to indicate the approximation order */
  PetscInt  ndof;               /* Number of degrees of freedom per node */

  /* Boundary conditions */
  PetscInt          nbc;        /* Number of boundary conditions */
  char            **bcname;     /* Names of boundary conditions (used to get entity sets from mesh) */
  PetscInt         *bcfaces;    /* Index of face on which to apply boundary condition */
  char             *bcfacesflip; /* 0 indicates the face normal coincides with region normal, 1 indicates flip */
  DohpBCOps **bcops;      /* Pointers to function tables defining how to apply the conditions */

  /* Element operations */
  PetscInt nelems;              /* Number of elements in this function space */
  DohpEFS **efs;                /* Pointers to the function space operations */
  DohpUInt8 *efssize;           /* Stored size at each element, stride is topological dimension of the element */
};

// For every Field, we need only the number of degrees of freedom and the local offset.
typedef struct {
  PetscInt localoffset;
  // PetscInt ndof;
} EField;

typedef struct {
  PetscInt *localoffset;
  PetscInt ndof;
} MField;

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


/* The Quotient map.  A quadrature rule and element coordinate mapping
 defines a quotient map on the Hilbert space $H^1(\Omega)$ which induces a
 finite topology.  The same quotient map can be used for many different function
 spaces $X \subset H^1$. */

// Since element quadrature and mapping contexts may have different sizes in
// different elements, we cannot simply use the PETSc data structures to handle
// local to global communication.  For now, we can just use flat arrays of
// pointers on the local process, but this will require some thought when we
// want to do dynamic load balancing.
struct _DohpQuotientOps {
  PetscErrorCode (*update)(DohpQuotient);
  PetscErrorCode (*setup)(DohpQuotient);
  PetscErrorCode (*setfromoptions)(DohpQuotient);
  PetscErrorCode (*view)(DohpQuotient,PetscViewer);
  PetscErrorCode (*destroy)(DohpQuotient);
};
struct _p_DohpQuotient {
  PETSCHEADER(struct _DohpQuotientOps);
  DohpMesh   mesh;
  DohpTag    qsizetag;
  DohpESH    loc;
  PetscInt   nelems;
  void     **quad;              // element quadrature context for locally owned elements
  void     **map;               // locally owned element map context, compatible with equad
  PetscInt   setupcalled;
};


PETSC_EXTERN_CXX_BEGIN
EXTERN PetscErrorCode DohpQuotientCreate_Gauss(DohpQuotient);
PETSC_EXTERN_CXX_END

#endif /* _DOHPIMPL_H */
