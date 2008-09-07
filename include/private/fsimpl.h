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

struct _dFSOps {
  DMOPS(dFS)
  dErr (*setfromoptions)(dFS);
};

struct _p_dFS {
  PETSCHEADER(struct _dFSOps);
  dMesh               mesh;
  dMeshTag            partition,degree;
  dMeshESH            active;
  struct dFSBoundary *bdy_start;
  dQuotient           quotient;
  dJacobi             jacobi;
  dBool               setupcalled;
  void               *data;
};

PETSC_EXTERN_CXX_END

#endif
