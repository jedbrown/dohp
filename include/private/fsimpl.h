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
  return !((stat & dFSBSTATUS_DIRICHLET) && ((stat & dFSBSTATUS_WEAK) || dFSBStatusStrongCount(stat))); /* cannot be Dirichlet and anything else */
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
  Sliced       slice;           /**< Lower level DM object, manages the global/local aspects with no additional structure */
  Sliced       dslice;          /**< Manages the Dirichlet local/global Dirichlet vectors */
  dInt         bs;              /**< Block size (number of dofs per node) */
  dMeshEH     *ents;            /**< All entities in active set */
  dInt         nelem;
  dInt        *off;             /**< Offset of element dofs in expanded vector */
  s_dRule     *rule;            /**< Integration rule */
  s_dEFS      *efs;             /**< Element function space, defined for all entities */
  dInt        *vtxoff;
  dReal       (*vtx)[3];
  dMeshTag     goffsetTag;      /**< Offset of first node in global vector */
  dMeshTag     gdoffsetTag;     /**< Offset of first node in global Dirichlet vector */
  dMeshTag     gcoffsetTag;     /**< Offset of first node in global closure vector */
  dMeshTag     loffsetTag;      /**< Offset of first node in local vector (split to include both real and Dirichlet vectors, based on dsplit) */
  dMeshESH     ownedExplicitSet,ghostExplicitSet; /**< Set of all entities that should be represented explicitly */
  dMeshESH     ownedDirichletSet,ghostDirichletSet; /**< Set of all entities that have full Dirichlet conditions (removed from global system) */
  dMeshESH     weakFaceSet;     /**< Faces on which weak forms need to be evaluated (I'm not sure this is actually needed) */
  VecScatter   ctod;            /**< Scatter from global closure (includes Dirichlet conditions) to local vectors with Dirichlet values */
  VecScatter   ctog;            /**< Scatter from global closure to global vector */
  Mat          E;               /**< full-order element assembly matrix (element nodes to local numbering) */
  Mat          Ep;              /**< preconditioning element assembly matrix (element nodes to local numbering, as sparse as possible) */
  Mat          Ed;              /**< element assembly matrix for Dirichlet nodes */
  Vec          weight;          /**< Vector in global space, used to compensate for overcounting after local to global */
  Vec          gc;              /**< Global closure vector */
  Vec          d,dl;            /**< Global and local Dirichlet vectors */
  dInt         maxQ;
  s_dFSWorkspace workspace[dFS_MAX_WORKSPACES];
  void        *data;
};

PETSC_EXTERN_CXX_END

#endif
