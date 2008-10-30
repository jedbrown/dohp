#ifndef _DOHPIMPL_H
#define _DOHPIMPL_H

#include "dohp.h"
#include "src/dm/dmimpl.h"
// #include "private/jacimpl.h"

PetscCookie dDM_COOKIE,dMFS_COOKIE;
PetscLogEvent dLOG_MatMult, dLOG_FunctionEval, dLOG_JacobianEval;

typedef struct {
  dErr (*dofspernode)(void /* what goes here? */);
  /* number of fields, number of nodes, normal vector, boundary values, element nodal values (in/out), output vector (out) */
  dErr (*apply)(dInt,dBool,const dScalar*,const dScalar*,dScalar*,dScalar*);
  dErr (*unapply)(dInt,dInt,const dScalar*,const dScalar*,dScalar*,dScalar*);
} DohpBCOps;

struct _DohpMFSOps {
  dErr (*lots)(void /* FIXME: what here? */);
};

struct _DohpDMOps {
  DMOps(DohpDM);
};

struct p_dohpDM {
  PETSCHEADER(struct _DohpDMOps);
  DMHEADER
  /* Distributed object */
  dMesh mesh;

  /* Defined on the locally owned portion of the mesh, can be reconstructed from information stored in the mesh. */
  dRule *ruleR;
  dRule *ruleF;
  dEFS *efsR;
  dEFS *efsF;
  DohpEMap_Hex  *emapR;
  DohpEMap_Quad *emapF;

  /* Approximation space */
  dInt   nmfs;              /* Number of function spaces in this DM */
  DohpMFS   *mfs;               /* handles of the spaces */
  dInt  *nfields;           /* [nmfs], number of fields on each function space */
  char ***fieldnames;           /* names of the fields in each function space */
  dInt   nelems;            /* total number of elements */
  dInt  *ind;               /* [nelems], starting index of interior dofs in local vector */
  dInt  *find;              /* index of face dofs in local vector */
  dInt  *foff;              /* [nelems], offset into find for each element */
  void     **edata;             /* private data (Jacobian storage) for each element */

  /* Blocks: during preconditioning, the */
  /* Preconditioning: this information is used for matrix-free application of blocks and assembly of low-order blocks */
  dInt **fieldoff;          /* [nfields][nelems], offset to add to ind[r] to get interior dofs in local vector */
  dInt **sind;              /* [nfields][nelems], starting index of interior dofs in split local vectors */
  dInt  *sfind;             /* index of face dofs in split local vector */
  dInt  *sfoff;             /* [nelems], offset into sfind for each element */
};

struct _DohpBlockOps {
  dErr (*stuffhere)(void /* FIXME: what here? */);
};

struct p_dBlock {
  PETSCHEADER(struct _DohpBlockOps);
  DohpDM dm;
  Sliced sliced;
  VecScatter unifiedtosplit;    /* Extract the relevant fields from a unified vector */
  VecScatter localtoglobal;     /* Scatter the local vector to the global split vector */

  dInt lmfs,rmfs;           /* Index of left and right MFS in  */
  dInt nlf,nrf;             /* Number of left fields, number of right fields */
  dInt *lfields,*rfields;   /* Indices of left and right fields in dm->[] */

  dInt nelems;              /* Total number of elements active in this block */
  dInt *ind;                /* [nelems] starting index of interior dofs in local vector */
};

struct p_dMFS {
  PETSCHEADER(struct _DohpMFSOps);
  dMesh  mesh;
  char     *name;               /* name of the domain, used to find the corresponding entity set in the mesh */
  char     *sizetagname;        /* a tag with this name is stored in the mesh to indicate the approximation order */
  dInt  ndof;               /* Number of degrees of freedom per node */

  /* Boundary conditions */
  dInt          nbc;        /* Number of boundary conditions */
  char            **bcname;     /* Names of boundary conditions (used to get entity sets from mesh) */
  dInt         *bcfaces;    /* Index of face on which to apply boundary condition */
  char             *bcfacesflip; /* 0 indicates the face normal coincides with region normal, 1 indicates flip */
  DohpBCOps **bcops;      /* Pointers to function tables defining how to apply the conditions */

  /* Element operations */
  dInt nelems;              /* Number of elements in this function space */
  dEFS **efs;                /* Pointers to the function space operations */
  dInt *efssize;           /* Stored size at each element, stride is topological dimension of the element */
};

// For every Field, we need only the number of degrees of freedom and the local offset.
typedef struct {
  dInt localoffset;
  // dInt ndof;
} EField;

typedef struct {
  dInt *localoffset;
  dInt ndof;
} MField;

// The facets of a line are 2 vertices.
// The projection matrix is the 1x1 identity so vertices are indexed straight out of the local vector.
typedef struct {
  dInt localoffset[2];
} ElemFacets_Line;

// The facets of a quad are 4 lines.
// F^* -> R^*: c[i] = SUM(j=0:size[fnum]) projrf[fnum][i*size[fnum]+j] * f[j]
// R^* -> F^*: f[i] = SUM(j=0:csize)      projfr[fnum][i*csize+j] * c[j]
typedef struct {
  dInt     size[4];         // Encode orientation and offset in the high-order bits of this value
  dScalar *projef[4];       // Element to facet projection
  dScalar *projfe[4];       // Facet to element projection
  dInt     facetoffset[4];  // Edge degrees of freedom are indexed out of the facet vector.
} ElemFacets_Quad;

// The facets of a hex are 6 quads.
typedef struct {
  dInt     size[6];         // How do we encode both x- and y-orientation in the high bits?
  dScalar *projef[6];       // Element to facet projection
  dScalar *projfe[6];       // Facet to element projection
  dInt     facetoffset[6];  // Face degrees of freedom are indexed out of the facet vector.
} ElemFacets_Hex;


// Since element quadrature and mapping contexts may have different sizes in
// different elements, we cannot simply use the PETSc data structures to handle
// local to global communication.  For now, we can just use flat arrays of
// pointers on the local process, but this will require some thought when we
// want to do dynamic load balancing.
struct _dQuotientOps {
  dErr (*update)(dQuotient);
  dErr (*setup)(dQuotient);
  dErr (*setfromoptions)(dQuotient);
  dErr (*view)(dQuotient,PetscViewer);
  dErr (*destroy)(dQuotient);
};


/**
* The Quotient map.  A quadrature rule and element coordinate mapping defines a quotient map on the Hilbert space
* $H^1(\Omega)$ which induces a finite topology (on the infinite dimensional space).  The same quotient map can be used
* for many different discrete function spaces $X \subset H^1$.
* */
struct p_dQuotient {
  PETSCHEADER(struct _dQuotientOps);
  dMesh                    mesh;
  dMeshTag                     qsizetag;
  dMeshESH                     loc;
  dInt                    nelems;
  dInt                   *degree;
  void                      **quad; // element quadrature context for locally owned elements
  void                      **emap; // locally owned element map context, compatible with quad
  dInt                    setupcalled;
  dQuotientSetDegreeFunc   setdegreefunc;
  void                       *setdegreectx;
  dBool                  setdegreeset;
};


PETSC_EXTERN_CXX_BEGIN
EXTERN dErr dQuotientCreate_Gauss(dQuotient);
PETSC_EXTERN_CXX_END

#endif /* _DOHPIMPL_H */
