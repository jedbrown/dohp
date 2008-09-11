#ifndef _DOHPTYPE_H
#define _DOHPTYPE_H

#include "petsc.h"
#include "iMesh.h"


/**
* These types all have to be exactly the Petsc versions.  These typedefs are here just to shorten the names, not to
* become autonomous.
* 
*/

typedef PetscInt       dInt;
typedef PetscReal      dReal;
typedef PetscScalar    dScalar;
typedef PetscTruth     dBool;
typedef PetscErrorCode dErr;
typedef PetscObject    dObject;
typedef PetscViewer    dViewer;

/* #define dEntTopology      enum iMesh_EntityTopology */
typedef enum iMesh_EntityTopology dEntTopology;
typedef enum iBase_EntityType dEntType;

#define dCHK(err) CHKERRQ(err);
#define dERROR(n,...) {return PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,n,1,__VA_ARGS__);}

#define dPrintf PetscPrintf
#define dMemcpy(a,b,c) PetscMemcpy(a,b,c)
#define dMemzero(a,b)  PetscMemzero(a,b)
#define dValidHeader(a,b,c) PetscValidHeaderSpecific(a,b,c)
#define dValidPointer(a,b) PetscValidPointer(a,b)
#define dMalloc(a,b) PetscMalloc(a,b)
#define dNew(a,b) PetscNew(a,b)
#define dNewLog(a,b,c) PetscNewLog(a,b,c)
#define dFree(a) PetscFree(a)
#define dNewM(n,t,p) (dMalloc((n)*sizeof(t),(p)) || dMemzero(*(p),(n)*sizeof(t)))
#define dMallocM(n,t,p) (dMalloc((n)*sizeof(t),(p)))

#define dMax(a,b) PetscMax(a,b)
#define dMin(a,b) PetscMin(a,b)
#define dSqr(a) PetscSqr(a)

#define dGamma(a) tgamma(a) /* This is defined in math.h as of C99. */

#ifndef false
# define false PETSC_FALSE
#endif
#ifndef true
# define true PETSC_TRUE
#endif

#define dMAX_PATH_LEN PETSC_MAX_PATH_LEN
#define dNAME_LEN     256

#define dFunctionBegin \
  {\
   if (petscstack && (petscstack->currentsize < PETSCSTACKSIZE)) {    \
    petscstack->function[petscstack->currentsize]  = __func__; \
    petscstack->file[petscstack->currentsize]      = __FILE__; \
    petscstack->directory[petscstack->currentsize] = __SDIR__; \
    petscstack->line[petscstack->currentsize]      = __LINE__; \
    petscstack->currentsize++; \
  }}

#define dFunctionReturn(a) \
  {\
  PetscStackPop; \
  return(a);}

#define dFunctionReturnVoid() \
  {\
  PetscStackPop; \
  return;}

#endif
