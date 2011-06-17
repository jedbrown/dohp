#ifndef _DOHPFSIMPL_H
#define _DOHPFSIMPL_H

#include "dohpfs.h"
#include "dohp.h"
#include <private/dmimpl.h>

dEXTERN_C_BEGIN

extern dLogEvent dLOG_Q1HexComputeQuadrature,dLOG_FSMatSetValuesExpanded;

static inline dInt dFSBStatusStrongCount(dFSBStatus stat) {
  return stat & dFSBSTATUS_MASK;
}
static inline dFSBStatus dFSBStatusSetStrongCount(dFSBStatus stat,dInt count) {
  return (stat & ~dFSBSTATUS_MASK) /* clear lower bits */ & count;
}
static inline bool dFSBStatusValid(dFSBStatus stat) {
  /* cannot be Dirichlet and anything else, note that (!a ^ !b) is (a LXOR b) */
  return (!(stat & dFSBSTATUS_DIRICHLET) ^ !((stat & dFSBSTATUS_WEAK) || dFSBStatusStrongCount(stat)));
}

struct dFSConstraintCtx {
  dFSConstraintFunction  cfunc;
  void                   *user;
};

struct dRulesetWorkspaceLink {
  dBool checkedout;
  dInt dof;
  dScalar *u,*v,*du,*dv;
  struct dRulesetWorkspaceLink *next;
};

struct dRulesetWorkspace {
  dScalar *q;
  dScalar *cjac;
  dScalar *cjinv;
  dScalar *jw;
  struct dRulesetWorkspaceLink *link;
};

struct _n_dRuleset {
  dInt refct;
  dMesh mesh;                   /**< Mesh on which the rules are defined */
  dMeshESH set;                 /**< Set containing all entities needing integration (and perhaps others of different type/topology) */
  dEntType type;                /**< Type of entities in the set on which integration is to be done */
  dEntTopology topo;            /**< Topology of entities in the set on which integration is to be done */
  dInt n;                       /**< Number of entities in set */
  dRule *rules;                 /**< Array of rules, one for each patch */
  dInt maxQ;                    /**< Max number of quadrature points in any patch */
  dInt maxnpatches;             /**< Max number of patches in any element */
  dInt maxQelem;                /**< Max number of quadrature points in any element */
  struct dRulesetWorkspace *workspace; /**< Manages workspace storage for evaluations based on dRules in this set */
};

struct _dFSIntegrationLink {
  char *name;
  dQuadrature quad;
  dRule *rule;
  dEFS  *efs;
  struct _dFSIntegrationLink *next;
};

struct _dFSOps {
  dErr (*impldestroy)(dFS);
  dErr (*setfromoptions)(dFS);
  dErr (*view)(dFS,dViewer);
  dErr (*buildspace)(dFS);
  dErr (*getsubelementmeshsize)(dFS,dInt*,dInt*,dInt*);
  dErr (*getsubelementmesh)(dFS,dInt,dInt,dEntTopology[],dInt[],dInt[]);
  dErr (*loadintofs)(dViewer,const char[],dFS);
};

struct _p_dFS {
  struct _p_DM dm;              /**< This must come first so that downcasting to DM works correctly */
  struct _dFSOps *ops;
  dMesh        mesh;
  struct {
    dMeshTag degree;            /**< Effective degree of elements, packed into dPolynomialOrder */
    dMeshTag boundary;          /**< Tag indicating which boundary condition should be enforced */
    dMeshTag bstatus;           /**< Boundary status tag, every NEUMANN_SET=x set will be tagged */
    dMeshTag bdyConstraint;     /**< User-defined context for enforcing boundary constraints */
    dMeshTag goffset;           /**< Offset of first node in global vector */
    dMeshTag gcoffset;          /**< Offset of first node in global closure vector */
    dMeshTag loffset;           /**< Offset of first node in local vector (split to include both real and Dirichlet vectors, based on dsplit) */
    dMeshTag partition;         /**< Part number of activeSet */
    dMeshTag orderedsub;        /**< Part number of ordered set */
  } tag;
  struct {
    dMeshESH active;            /**< All entities that will be part of this space, weak forms are evaluated on regions in this set */
    dMeshESH boundaries;        /**< Entities for which boundary conditions need to be enforced */
    dMeshESH ordered;           /**< All entities in local space, ordered as they are used in the computation */
    dMeshESH explicit,dirichlet,ghost;  /**< Vectors have form [explicit,dirichlet,ghost].  Global is first
                                         * part, closure is first two, local is whole thing.  The sets allow
                                         * us to work with pieces. */
    dMeshESH weakFace;
  } set;
  dJacobi      jacobi;
  char         bdyTagName[dNAME_LEN]; /**< Usually "NEUMANN_SET" */
  dBool        spacebuilt;
  dBool        assemblefull;    /**< Use full order constraints for assembly */
  dBool        assemblereduced; /**< Assemble only diagonal part of blocks, only matters for bs>1 and MATAIJ */
  dInt         nelem;
  dInt        *off;             /**< Offset of element dofs in expanded vector */
  struct _dFSIntegrationLink *integration;
  dInt         n,nc,ngh;        /**< Vector sizes in blocks: owned, owned closure, ghosts */
  Mat          E;               /**< full-order element assembly matrix (element nodes to local numbering) */
  Mat          Ep;              /**< preconditioning element assembly matrix (element nodes to local numbering, as sparse as possible) */
  Vec          weight;          /**< Vector in global space, used to compensate for overcounting after local to global */
  Vec          gvec;            /**< Global Vec, closure can be obtained with VecDohpGetClosure() */
  Vec          dcache;          /**< All Dirichlet values, this is only a cache so that we can project a vector into the inhomogeneous space */
  VecScatter   dscat;           /**< Scatter from global closure to \a dcache. */

  struct {
    Vec expanded;               /**< expanded, used for integration */
    Vec global;                 /**< coordinates at all nodes of geometryfs */
    dFS fs;                     /**< function space for geometry */
  } geometry;
  struct {
    Vec expanded;
    Vec global;
    dFS fs;
  } nodalcoord;

  dFSRotation  rot;             /**< Rotation for local vector */
  char         orderingtype[256];
  char       **fieldname;
  dUnit       *fieldunit;
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

extern dErr dFSCreateLocalToGlobal_Private(dFS fs,dInt n,dInt nc,dInt ngh,dInt *ghidx,dInt rstart);
extern dErr dFSPopulatePartitionedSets_Private(dFS FS,dMeshAdjacency meshadj);
extern dErr dFSBuildSpaceOffsets_Private(dFS fs,dMeshTag indexTag,const dInt inodes[],dInt rstart,dInt crstart,dInt nents,dMeshEH ents[],dInt *ghstart);
extern dErr dFSBuildSpaceVectors_Private(dFS fs,dMeshTag indexTag,const dInt inodes[],dInt rstart,dInt ghents_s,const dMeshEH ghents[]);
extern dErr dFSBuildSpaceWithOrderedSet_Private(dFS fs,dMeshAdjacency meshAdj);
extern dErr DMDestroy_dFS(DM);

dEXTERN_C_END

#endif
