#ifndef _FSIMPL_H
#define _FSIMPL_H

#include "dohpfs.h"
#include "private/dmimpl.h"

PETSC_EXTERN_CXX_BEGIN

extern PetscLogEvent dLOG_Q1HexComputeQuadrature,dLOG_FSMatSetValuesExpanded;

struct _p_dFSBoundary {
  dInt                           nDirichlet;
  dInt                           nGlobal;
  dMeshManifold                  manifold;
  dTruth                         flip;
  dFSBoundaryConstraintFunction  cfunc;
  void                          *user;
  struct _p_dFSBoundary         *next;
};

typedef struct {
  dInt Q;
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

struct _p_dFS {
  PETSCHEADER(struct _dFSOps);
  DMHEADER
  dMesh        mesh;
  dMeshTag     degreetag,ruletag; /**< tags on regions */
  dMeshTag     fieldSpecTag;      /**< 2-int tag on every entity in active set */
  dMeshESH     active;          /**< regions that will be part of this space */
  dFSBoundary  bdylist;
  dQuotient    quotient;
  dJacobi      jacobi;
  dTruth       spacebuilt;
  dTruth       assemblefull;    /**< Use full order constraints for assembly */
  dInt         ruleStrength;
  Sliced       sliced;
  dInt         nbdofs;
  dInt         n,N;             /**< length of the owned and global vectors */
  dInt         nlocal;          /**< number of owned+ghost dofs on this process */
  dInt         rstart;          /**< global offset of first owned dof */
  dMeshEH     *ents;            /**< All entities in active set
  dInt         m;               /**< Number of expanded dofs */
  dInt         D;               /**< Number of dofs per (non-boundary) node */
  dInt         nelem;
  dInt        *off;             /**< Offset of element dofs in expanded vector */
  s_dRule     *rule;            /**< Integration rule */
  s_dEFS      *efs;             /**< Element function space, defined for all entities */
  dInt        *vtxoff;
  dReal       (*vtx)[3];
  Mat          C;               /**< full-order constraint matrix (element dofs to local numbering) */
  Mat          Cp;              /**< preconditioning constraint matrix (element dofs to local numbering, as sparse as possible) */
  Mat          Cd;              /**< constraint matrix for Dirichlet values */
  Vec          weight;          /**< Vector in global space, used to compensate for overcounting after local to global */
  s_dFSWorkspace *workspace;
  void        *data;
};

PETSC_EXTERN_CXX_END

#endif
