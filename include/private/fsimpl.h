#ifndef _FSIMPL_H
#define _FSIMPL_H

#include "dohpfs.h"
#include "src/dm/dmimpl.h"

PETSC_EXTERN_CXX_BEGIN

extern PetscLogEvent dLOG_FSConstrain;

struct dFSBoundary {
  PF                  constrain;
  char               *name;
  dMeshESH            facets;
  dMeshTag            orient;
  dBool               fliporient;
  struct dFSBoundary *next;
};

struct dMeshInterface {
  int rank;
  dMeshESH owned,remote;
  UT_hash_handle hh;
};

struct _dFSOps {
  DMOPS(dFS)
  dErr (*setfromoptions)(dFS);
  dErr (*impldestroy)(dFS);
  dErr (*buildspace)(dFS);
};

struct _p_dFS {
  PETSCHEADER(struct _dFSOps);
  dMesh               mesh;
  dMeshTag            degree,ruletag; /* tags on regions */
  dMeshESH            active;         /* regions that will be part of this space */
  struct dMeshInterface  *iface;      /* for each process sharing an interface, sets of owned and remote interface entities */
  struct dFSBoundary *bdylist;
  dQuotient           quotient;
  dJacobi             jacobi;
  dBool               spacebuilt;
  Sliced              sliced;
  MeshListEH          r,f,e,v;  /**< region, face, edge, vertex */
  PetscInt            n,N;      /**< length of the local and global vectors */
  void               *data;
};

PETSC_EXTERN_CXX_END

#endif
