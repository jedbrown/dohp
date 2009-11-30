#ifndef _DOHPTYPE_H
#define _DOHPTYPE_H

#include "dohpconfig.h"
#include <petsc.h>
#include <iMesh.h>

#if defined __cplusplus
#  define dEXTERN_C_BEGIN extern "C" {
#  define dEXTERN_C_END   }
#else
#  define dEXTERN_C_BEGIN
#  define dEXTERN_C_END
#endif

/**
* These types all have to be exactly the Petsc versions.  These typedefs are here just to shorten the names, not to
* become autonomous.
*
*/

typedef PetscInt       dInt;
typedef PetscReal      dReal;
typedef PetscScalar    dScalar;
typedef PetscTruth     dBool;
typedef PetscTruth     dTruth;
typedef PetscErrorCode dErr;
typedef PetscObject    dObject;
typedef PetscViewer    dViewer;

typedef PetscMPIInt dMPIInt;

typedef int dMeshInt;
typedef double dMeshReal;
typedef iBase_EntityHandle dMeshEH;
typedef iBase_TagHandle dMeshTag;
typedef iBase_EntitySetHandle dMeshESH;

/* ITAPS data types */
typedef int dIInt;
typedef double dIReal;
typedef char dIByte;

typedef enum { dDATA_INT, dDATA_REAL, dDATA_EH, dDATA_BYTE } dDataType;
typedef enum { dTOPO_POINT, dTOPO_LINE, dTOPO_POLYGON, dTOPO_TRIANGLE,
               dTOPO_QUAD, dTOPO_POLYHEDRON, dTOPO_TET, dTOPO_HEX, dTOPO_PRISM,
               dTOPO_PYRAMID, dTOPO_SEPTAHEDRON, dTOPO_ALL } dEntTopology;
typedef enum { dTYPE_VERTEX, dTYPE_EDGE, dTYPE_FACE, dTYPE_REGION, dTYPE_ALL } dEntType;

typedef unsigned char dEntStatus;

#endif  /* _DOHPTYPE_H */
