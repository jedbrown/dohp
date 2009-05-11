#ifndef _FSIMPL_H
#define _FSIMPL_H

#include "dohpfs.h"
#include "private/dmimpl.h"

PETSC_EXTERN_CXX_BEGIN

extern PetscLogEvent dLOG_Q1HexComputeQuadrature,dLOG_FSMatSetValuesExpanded;

static inline dInt dFSBStatusStrongCount(dFSBStatus stat) {
  return stat & dFSBSTATUS_MASK;
}
static inline dFSBStatus dFSBStatusSetStrongCount(dFSBStatus stat,dInt count) {
  return (stat & ~dFSBSTATUS_MASK) /* clear lower bits */ & count;
}
static inline dTruth dFSBStatusValid(dFSBStatus stat) {
  /* cannot be Dirichlet and anything else, note that (!a ^ !b) is (a LXOR b) */
  return (!(stat & dFSBSTATUS_DIRICHLET) ^ !((stat & dFSBSTATUS_WEAK) || dFSBStatusStrongCount(stat)));
}

struct dFSConstraintCtx {
  dFSConstraintFunction  cfunc;
  void                   *user;
};

typedef struct {
  dInt status;                  /* 0:unallocated 1:available 2:checked out */
  char name[dNAME_LEN];
  dReal (*q)[3];
  dReal (*jinv)[3][3];
  dReal *jw;
  dScalar *u,*v,*du,*dv;
} s_dFSWorkspace;

struct _dFSOps {
  DMOPS(dFS)
  dErr (*setfromoptions)(dFS);
  dErr (*impldestroy)(dFS);
  dErr (*buildspace)(dFS);
};

#define dFS_MAX_WORKSPACES 64

struct _p_dFS {
  PETSCHEADER(struct _dFSOps);
  DMHEADER
  dMesh        mesh;
  dMeshAdjacency meshAdj;
  dMeshTag     degreetag,ruletag; /**< tags on regions */
  dMeshESH     activeSet;       /**< all entities that will be part of this space, weak forms are evaluated on regions in this set */
  dQuotient    quotient;
  dJacobi      jacobi;
  dMeshESH     boundaries;      /**< Set of boundary conditions to be enforced on this function space */
  dMeshTag     bdyTag;          /**< Tag for each boundary condition */
  char         bdyTagName[dNAME_LEN]; /**< Usually "NEUMANN_SET" */
  dMeshTag     bstatusTag;      /**< Boundary status tag, every NEUMANN_SET=x set will be tagged */
  dMeshTag     bdyConstraintTag; /**< User-defined context for enforcing boundary constraints */
  dTruth       spacebuilt;
  dTruth       assemblefull;    /**< Use full order constraints for assembly */
  dInt         ruleStrength;
  dInt         bs;              /**< Block size (number of dofs per node) */
  dMeshEH     *ents;            /**< All entities in active set */
  dInt         nelem;
  dInt        *off;             /**< Offset of element dofs in expanded vector */
  s_dRule     *rule;            /**< Integration rule */
  s_dEFS      *efs;             /**< Element function space, defined for all entities */
  dInt        *vtxoff;
  dReal       (*vtx)[3];
  dInt         n,nc,ngh;        /**< Vector sizes in blocks: owned, owned closure, ghosts */
  dMeshTag     goffsetTag;      /**< Offset of first node in global vector */
  dMeshTag     gcoffsetTag;     /**< Offset of first node in global closure vector */
  dMeshTag     loffsetTag;      /**< Offset of first node in local vector (split to include both real and Dirichlet vectors, based on dsplit) */
  /**< Vectors have form [explicit,dirichlet,ghost].  Global is first part, closure is first two, local is whole thing.
  * The sets allow us to work with pieces. */
  dMeshESH     explicitSet,dirichletSet,ghostSet;
  dMeshESH     weakFaceSet;     /**< Faces on which weak forms need to be evaluated (I'm not sure this is actually needed) */
  Mat          E;               /**< full-order element assembly matrix (element nodes to local numbering) */
  Mat          Ep;              /**< preconditioning element assembly matrix (element nodes to local numbering, as sparse as possible) */
  Vec          weight;          /**< Vector in global space, used to compensate for overcounting after local to global */
  Vec          gvec;            /**< Global Vec, closure can be obtained with VecDohpGetClosure() */
  Vec          dcache;          /**< All Dirichlet values, this is only a cache so that we can project a vector into the inhomogeneous space */
  VecScatter   dscat;           /**< Scatter from global closure to \a dcache. */
  ISLocalToGlobalMapping bmapping; /**< Block mapping, Dirichlet blocks have negative global index */
  dInt         maxQ;
  dFSRotation  rot;             /**< Rotation for local vector */
  s_dFSWorkspace workspace[dFS_MAX_WORKSPACES];
  char         orderingtype[256];
  void        *data;
};

struct _dFSRotationOps {        /* We might not ever have specialized implementations of these */
  dErr (*apply)(dFSRotation,Vec,dFSRotateMode,dFSHomogeneousMode);
  dErr (*applylocal)(dFSRotation,Vec,dFSRotateMode,dFSHomogeneousMode);
};

struct _p_dFSRotation {
  PETSCHEADER(struct _dFSRotationOps);
  dInt   bs;                    /**< block size */
  dInt   n;                     /**< number of blocks to rotate (size of \a is) */
  IS     is;                    /**< which blocks to rotate */
  dReal *rmat;                  /**< rotation matrices for the blocks in \a is */
  dInt  *nstrong;               /**< number of dofs to enforce strongly, for each block is \a is */
  Vec    strong;                /**< Values for all strongly enforced dofs */
};

PETSC_EXTERN_CXX_END

#endif
