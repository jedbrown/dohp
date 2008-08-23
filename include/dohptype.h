#ifndef _DOHPTYPE_H
#define _DOHPTYPE_H

#include "petsc.h"


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

#define dPrintf PetscPrintf
#define dMemcpy(a,b,c) PetscMemcpy(a,b,c)
#define dMemzero(a,b)  PetscMemzero(a,b)
#define dValidHeader(a,b,c) PetscValidHeaderSpecific(a,b,c)
#define dValidPointer(a,b) PetscValidPointer(a,b)
#define dMalloc(a,b) PetscMalloc(a,b)
#define dCHK(err) CHKERRQ(err);

#ifndef false
# define false PETSC_FALSE
#endif

#ifndef true
# define true PETSC_TRUE
#endif

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

#define dERROR(n,...) {return PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,n,1,__VA_ARGS__);}

#endif
