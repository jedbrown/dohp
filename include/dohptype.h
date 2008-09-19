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
typedef unsigned char dGeomOrient;

typedef int dMeshInt;
typedef double dMeshReal;
typedef iBase_EntityHandle dMeshEH;
typedef iBase_TagHandle dMeshTag;
typedef iBase_EntitySetHandle dMeshESH;

/* ITAPS data types */
typedef int dIInt;
typedef double dIReal;
typedef char dIByte;

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
#define dCalloc(n,p) (dMalloc((n),(p)) || dMemzero(*(p),(n)))

#define dValidPointer2(a,b,c,d) (dValidPointer((a),(b)) || dValidPointer((c),(d)))
#define dValidPointer3(a,b,c,d,e,f) (dValidPointer2((a),(b),(c),(d)) || dValidPointer((e),(f)))
#define dValidPointer4(a,b,c,d,e,f,g,h) (dValidPointer3((a),(b),(c),(d),(e),(f)) || dValidPointer((g),(h)))
#define dValidPointer5(a,b,c,d,e,f,g,h,i,j) (dValidPointer4((a),(b),(c),(d),(e),(f),(g),(h)) || dValidPointer((i),(j)))
#define dValidPointer6(a,b,c,d,e,f,g,h,i,j,k,l) (dValidPointer5((a),(b),(c),(d),(e),(f),(g),(h),(i),(j)) || dValidPointer((k),(l)))
#define dValidPointer7(a,b,c,d,e,f,g,h,i,j,k,l,m,n) (dValidPointer6((a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l)) || dValidPointer((m),(n)))

#define dValidPointerSpecific(p,t,a)                              \
  {if (!p) dERROR(PETSC_ERR_ARG_BADPTR,"Null Pointer: Parameter # %d",(a)); \
    if ((size_t)(p) % sizeof(*(p))) dERROR(PETSC_ERR_ARG_BADPTR,"Insufficient alignment for pointer to %s: Parameter # %d should have %ld alignment",(t),(a),sizeof(*(p)));}
#define dValidPointerSpecific2(p0,t0,a0,p1,t1,a1) {dValidPointerSpecific(p0,t0,a0); dValidPointerSpecific(p1,t1,a1);}
#define dValidPointerSpecific3(p0,t0,a0,p1,t1,a1,p2,t2,a2) \
  {dValidPointerSpecific(p0,t0,a0); dValidPointerSpecific2(p1,t1,a1,p2,t2,a2);}
#define dValidPointerSpecific4(p0,t0,a0,p1,t1,a1,p2,t2,a2,p3,t3,a3) \
  {dValidPointerSpecific(p0,t0,a0); dValidPointerSpecific3(p1,t1,a1,p2,t2,a2,p3,t3,a3);}
#define dValidPointerSpecific5(p0,t0,a0,p1,t1,a1,p2,t2,a2,p3,t3,a3,p4,t4,a4) \
  {dValidPointerSpecific(p0,t0,a0); dValidPointerSpecific4(p1,t1,a1,p2,t2,a2,p3,t3,a3,p4,t4,a4);}
#define dValidPointerSpecific6(p0,t0,a0,p1,t1,a1,p2,t2,a2,p3,t3,a3,p4,t4,a4,p5,t5,a5) \
  {dValidPointerSpecific(p0,t0,a0); dValidPointerSpecific5(p1,t1,a1,p2,t2,a2,p3,t3,a3,p4,t4,a4,p5,t5,a5);}
#define dValidPointerSpecific7(p0,t0,a0,p1,t1,a1,p2,t2,a2,p3,t3,a3,p4,t4,a4,p5,t5,a5,p6,t6,a6) \
  {dValidPointerSpecific(p0,t0,a0); dValidPointerSpecific6(p1,t1,a1,p2,t2,a2,p3,t3,a3,p4,t4,a4,p5,t5,a5,p6,t6,a6);}

#define dValidIntPointer(p,a) PetscValidIntPointer(p,a)
#define dValidHandlePointer(p,a) dValidPointerNamedSpecific(p,void*,"void*",a)
#define dValidScalarPointer(p,a) dValidPointerNamedSpecific(p,dScalar,"dScalar",a)
#define dValidRealPointer(p,a) dValidPointerNamedSpecific(p,dReal,"dReal",a)

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
#define dSTR_LEN      256

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

#endif  /* _DOHPTYPE_H */
