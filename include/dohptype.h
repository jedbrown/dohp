#ifndef _DOHPTYPE_H
#define _DOHPTYPE_H

#define _POSIX_C_SOURCE 199309L

#include <dohpconfig.h>
#include <petscconf.h>
#include <iMesh.h>
#include <mpi.h>
#include <stddef.h>

#if defined __cplusplus
#  define dEXTERN_C_BEGIN extern "C" {
#  define dEXTERN_C_END   }
#  define restrict
#else
#  define dEXTERN_C_BEGIN
#  define dEXTERN_C_END
#  include <stdbool.h>
#endif

dEXTERN_C_BEGIN

/**
* These types all have to be exactly the PETSc versions.  These typedefs are here so that we don't alway have to include
* all the PETSc headers (and to give us shorter names).
*
*/

#if defined(PETSC_USE_64BIT_INDICES)
typedef long long dInt;
#else
typedef int dInt;
#endif

#if defined PETSC_USE_SCALAR_SINGLE
typedef float dReal;
#elif defined PETSC_USE_SCALAR_LONG_DOUBLE
typedef long long dReal;
#else
typedef double dReal;
#endif
#if defined PETSC_USE_COMPLEX
#  include <complex.h>
#  if defined PETSC_USE_SCALAR_SINGLE
typedef float complex dScalar;
#  elif defined PETSC_USE_SCALAR_LONG_DOUBLE
typedef long long complex dScalar;
#  else
typedef double complex dScalar;
#  endif
#else
#  if defined PETSC_USE_SCALAR_SINGLE
typedef float dScalar;
#  elif defined PETSC_USE_SCALAR_LONG_DOUBLE
typedef long long dScalar;
#  else
typedef double dScalar;
#  endif
#endif

/* Do not define a boolean type */
typedef int dErr;               /* PetscErrorCode */
typedef int dClassId;            /* PetscClassId */
typedef int dLogEvent;          /* PetscLogEvent */
typedef struct _n_PetscFList  *dFList;
typedef struct _p_PetscObject *dObject;
typedef struct _p_PetscViewer *dViewer;

typedef int dMPIInt;

typedef int dMeshInt;
typedef double dMeshReal;
typedef iBase_EntityHandle dMeshEH;
typedef iBase_TagHandle dMeshTag;
typedef iBase_EntitySetHandle dMeshESH;

/* ITAPS data types */
typedef int dIInt;
typedef double dIReal;
typedef char dIByte;

typedef enum { dDATA_BYTE, dDATA_INT, dDATA_REAL, dDATA_EH, dDATA_ESH, dDATA_UB } dDataType;
typedef enum { dTOPO_POINT, dTOPO_LINE, dTOPO_POLYGON, dTOPO_TRIANGLE,
               dTOPO_QUAD, dTOPO_POLYHEDRON, dTOPO_TET, dTOPO_HEX, dTOPO_PRISM,
               dTOPO_PYRAMID, dTOPO_SEPTAHEDRON, dTOPO_ALL } dEntTopology;
typedef enum { dTYPE_VERTEX, dTYPE_EDGE, dTYPE_FACE, dTYPE_REGION, dTYPE_ALL } dEntType;

extern const char *dMeshEntTopologyName(dEntTopology);
extern const char *dMeshEntTypeName(dEntType);
extern dEntType dMeshEntTypeFromTopology(dEntTopology);

typedef unsigned char dEntStatus;

#define dUNUSED PETSC_UNUSED

dEXTERN_C_END

#endif  /* _DOHPTYPE_H */
