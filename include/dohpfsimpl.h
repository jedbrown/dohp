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

typedef struct {
  dInt status;                  /* 0:unallocated 1:available 2:checked out */
  char name[dNAME_LEN];
  dReal (*q)[3];
  dReal (*jinv)[3][3];
  dReal *jw;
  dScalar *u,*v,*du,*dv;
} s_dFSWorkspace;

struct _n_dRuleSet {
  dMesh mesh;
  dMeshESH set;
  dEntType type;
  dEntTopology topo;
  dInt n;
  dRule *rules;
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

#define dFS_MAX_WORKSPACES 64

struct _p_dFS {
  struct _p_DM dm;              /**< This must come first so that downcasting to DM works correctly */
  struct _dFSOps *ops;
  dMesh        mesh;
  struct {
    dMeshTag degree;            /**< Effective degree of elements, packed into dPolynomialOrder */
    dMeshTag rule;              /**< Degree of rules defined on elements that integration is to be performed on */
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
  dInt         ruleStrength;
  dInt         bs;              /**< Block size (number of dofs per node) */
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
  ISLocalToGlobalMapping bmapping; /**< Block mapping, Dirichlet blocks have negative global index */
  ISLocalToGlobalMapping mapping;  /**< Scalar mapping, mapping[i] = bmapping[i/bs]*bs+i%bs; */

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

  dInt         maxQ;
  dFSRotation  rot;             /**< Rotation for local vector */
  s_dFSWorkspace workspace[dFS_MAX_WORKSPACES];
  char         orderingtype[256];
  char       **fieldname;
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

dEXTERN_C_END

#endif
