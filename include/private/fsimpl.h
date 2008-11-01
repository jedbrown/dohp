#ifndef _FSIMPL_H
#define _FSIMPL_H

#include "dohpfs.h"
#include "src/dm/dmimpl.h"

PETSC_EXTERN_CXX_BEGIN

extern PetscLogEvent dLOG_FSConstrain;

struct _p_dFSBoundary {
  PF                  constrain;
  char               *name;
  dMeshESH            entset;
  dMeshTag            orient;
  dBdyType            btype;
  dTruth              fliporient;
  struct _p_dFSBoundary *next;
};

struct _dFSOps {
  DMOPS(dFS)
  dErr (*setfromoptions)(dFS);
  dErr (*impldestroy)(dFS);
  dErr (*buildspace)(dFS);
};

struct _p_dFS {
  PETSCHEADER(struct _dFSOps);
  DMHEADER
  dMesh               mesh;
  dMeshTag            degreetag,ruletag; /**< tags on regions */
  dMeshESH            active;   /**< regions that will be part of this space */
  dFSBoundary         bdylist;
  dQuotient           quotient;
  dJacobi             jacobi;
  dTruth              spacebuilt;
  Sliced              sliced;
  dInt                nbdofs;
  dInt                n,N;      /**< length of the owned and global vectors */
  dInt                nlocal;   /**< number of owned+ghost dofs on this process */
  dInt                rstart;   /**< global offset of first owned dof */
  dInt                m;        /**< Number of expanded dofs */
  dInt                D;        /**< Number of dofs per (non-boundary) node */
  dRule              *rule;     /**< Integration rule */
  dEFS               *efs;      /**< Element function space, defined for all entities */
  Mat                 C;        /**< full-order constraint matrix (element dofs to local numbering) */
  Mat                 Cp;       /**< preconditioning constraint matrix (element dofs to local numbering, as sparse as possible) */
  void               *data;
};

PETSC_EXTERN_CXX_END

#endif
