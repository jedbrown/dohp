#if !defined _DOHP_H
#define _DOHP_H

#include "dohptype.h"
#include <stdint.h>
#include <petscsys.h>

#define dSTATUS_UNOWNED   (dEntStatus)0x1
#define dSTATUS_SHARED    (dEntStatus)0x2
#define dSTATUS_INTERFACE (dEntStatus)0x4
#define dSTATUS_GHOST     (dEntStatus)0x8

#define dDataTypeToITAPS(dtype,itype) (*(itype) = (dtype), 0)
#define dEntTopoToITAPS(dtopo,itopo) (*(itopo) = (dtopo), 0)
#define dEntTypeToITAPS(dtype,itype) (*(itype) = (dtype), 0)

#define dCHK(err) do {if (PetscUnlikely(err)) return PetscError(PETSC_COMM_SELF,__LINE__,__func__,__FILE__,__SDIR__,(err),PETSC_ERROR_REPEAT," ");} while (0)
#define dERROR(n,...) return PetscError(PETSC_COMM_SELF,__LINE__,__func__,__FILE__,__SDIR__,(n),PETSC_ERROR_INITIAL,__VA_ARGS__)

#define dPrintf PetscPrintf
#define dMemcpy(a,b,c) PetscMemcpy(a,b,c)
#define dMemzero(a,b)  PetscMemzero(a,b)
#define dValidHeader(a,b,c) PetscValidHeaderSpecific(a,b,c)

#define dValidPointer(a,b) do {                                              \
    if (!(a)) dERROR(PETSC_ERR_ARG_NULL,"Null Pointer: Parameter # %d",(b)); \
    if ((size_t)a & 3) dERROR(PETSC_ERR_ARG_BADPTR,"Invalid Pointer: Parameter # %d",(b)); \
  } while (0)

#define dMalloc(a,b) PetscMalloc(a,b)
#define dNew(a,b) PetscNew(a,b)
#define dNewLog(a,b,c) PetscNewLog(a,b,c)
#define dFree(a) PetscFree(a)
#define dNewM(n,t,p) (dMalloc((n)*sizeof(t),(p)) || dMemzero(*(p),(n)*sizeof(t)))
#define dMallocM(n,t,p) (dMalloc((n)*sizeof(t),(p)))
#define dMallocA(n,p) (dMalloc((n)*sizeof(**(p)),(p)))
#define dCalloc(n,p) (dMalloc((n),(p)) || dMemzero(*(p),(n)))
#define dCallocA(n,p) (dCalloc((n)*sizeof(**(p)),(p)))

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

#define dValidCharPointer(p,a) PetscValidCharPointer((p),(a))
#define dValidIntPointer(p,a) PetscValidIntPointer((p),(a))
#define dValidHandlePointer(p,a) dValidPointerNamedSpecific((p),void*,"void*",(a))
#define dValidScalarPointer(p,a) dValidPointerSpecific((p),"dScalar",(a))
#define dValidRealPointer(p,a) dValidPointerSpecific((p),"dReal",(a))

#define dNonNullElse(a,b) ((a)?(a):(b))
#define dStrlen(s,l) PetscStrlen((s),(l))
extern dErr dObjectGetComm(dObject obj,MPI_Comm *comm);

extern dErr dRealTableView(dInt m,dInt n,const dReal mat[],const char *name,dViewer viewer);

static inline dInt dMaxInt(dInt a,dInt b) { return (a > b) ? a : b; }
static inline dInt dMinInt(dInt a,dInt b) { return (a < b) ? a : b; }
static inline dReal dMax(dReal a,dReal b) { return (a > b) ? a : b; }
static inline dReal dMin(dReal a,dReal b) { return (a < b) ? a : b; }
static inline dReal dAbs(dScalar a) { return fabs(a); }
static inline dScalar dSqr(dScalar a) { return a * a; }
static inline dInt dSqrInt(dInt a) { return a * a; }
static inline dReal dSqrt(dReal a) { return sqrt(a); }
static inline dScalar dDotScalar3(const dScalar a[3],const dScalar b[3]) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static inline dScalar dColonSymScalar3(const dScalar a[6],const dScalar b[6])
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + 2*a[3]*b[3] + 2*a[4]*b[4] + 2*a[5]*b[5]; }
static inline void dTensorSymCompress3(const dScalar Du[9],dScalar Dv[6])
{Dv[0] = Du[0]; Dv[1] = Du[4]; Dv[2] = Du[8]; Dv[3] = (Du[1]+Du[3])/2; Dv[4] = (Du[2]+Du[6]); Dv[5] = (Du[5]+Du[7])/2;}
static inline void dTensorSymUncompress3(const dScalar Du[6],dScalar Dv[9])
{Dv[0] = Du[0]; Dv[4] = Du[1]; Dv[8] = Du[2]; Dv[1] = Dv[3] = Du[3]; Dv[2] = Dv[6] = Du[4]; Dv[5] = Dv[7] = Du[5];}

#define dGamma(a) tgamma(a) /* This is defined in math.h as of C99. */

#ifndef false
# define false PETSC_FALSE
#endif
#ifndef true
# define true PETSC_TRUE
#endif

/* stdbool.h has small (1 byte) bools, PETSc uses an enum which has few size guarantees, so we use it directly and keep
* it out of our public interface (which is why these typedefs are here and not in dohptype.h).
**/
typedef PetscBool dTruth;
typedef PetscBool  dBool;
#define dTRUE  PETSC_TRUE
#define dFALSE PETSC_FALSE

#define dMAX_PATH_LEN PETSC_MAX_PATH_LEN
#define dNAME_LEN     256
#define dSTR_LEN      256

#define dUNUSED PETSC_UNUSED

#define dCACHE_LINE    64l       /* my cache lines are 64 bytes long */
#define dRPCL dCACHE_LINE/sizeof(dReal)
#define dSPCL dCACHE_LINE/sizeof(dScalar)

#define dDEFAULT_ALIGN 16l       /* SSE instructions require 16 byte alignment */

#define dNextCacheAligned(p) dNextAlignedAddr(dCACHE_LINE,(p))
#define dNextAligned(p)      dNextAlignedAddr(dDEFAULT_ALIGN,(p))

/** Returns the next address satisfying the given alignment.
*
* This function cannot fail.
*
* @param alignment must be a power of 2
* @param ptr The pointer
*
* @return aligned address
*/
static inline void *dNextAlignedAddr(size_t alignment,void *ptr)
{
  const uintptr_t base = (uintptr_t)ptr;
  const uintptr_t mask = (uintptr_t)alignment-1;
  return (void*)((base + mask) & ~mask);
}

/* Needs to be a macro because of pointer arithmetic */
#define dMemClaim(mem,n,p) do {                         \
    (p) = dNextAligned(mem);                            \
    (mem) = dNextAligned((p) + (n));                    \
  } while (0)

#if defined(PETSC_USE_INFO)
# define dInfo(a,s,...) PetscInfo_Private(__func__,(a),(s),__VA_ARGS__)
#else
# define dInfo(a,s,...) 0
#endif

#if defined(PETSC_USE_DEBUG)
# define dMallocA2(n0,p0,n1,p1) (dMallocA((n0),(p0)) || dMallocA((n1),(p1)))
# define dMallocA3(n0,p0,n1,p1,n2,p2) (dMallocA((n0),(p0)) || dMallocA2((n1),(p1),(n2),(p2)))
# define dMallocA4(n0,p0,n1,p1,n2,p2,n3,p3) (dMallocA((n0),(p0)) || dMallocA3((n1),(p1),(n2),(p2),(n3),(p3)))
# define dMallocA5(n0,p0,n1,p1,n2,p2,n3,p3,n4,p4)                       \
  (dMallocA((n0),(p0)) || dMallocA4((n1),(p1),(n2),(p2),(n3),(p3),(n4),(p4)))
# define dMallocA6(n0,p0,n1,p1,n2,p2,n3,p3,n4,p4,n5,p5)                  \
  (dMallocA((n0),(p0)) || dMallocA5((n1),(p1),(n2),(p2),(n3),(p3),(n4),(p4),(n5),(p5)))
# define dMallocA7(n0,p0,n1,p1,n2,p2,n3,p3,n4,p4,n5,p5,n6,p6)            \
  (dMallocA((n0),(p0)) || dMallocA6((n1),(p1),(n2),(p2),(n3),(p3),(n4),(p4),(n5),(p5),(n6),(p6)))
# define dFree2(a,b) (dFree(a) || dFree(b))
# define dFree3(a,b,c) (dFree(a) || dFree(b) || dFree(c))
# define dFree4(a,b,c,d) (dFree(a) || dFree(b) || dFree(c) || dFree(d))
# define dFree5(a,b,c,d,e) (dFree(a) || dFree(b) || dFree(c) || dFree(d) || dFree(e))
# define dFree6(a,b,c,d,e,f) (dFree(a) || dFree(b) || dFree(c) || dFree(d) || dFree(e) || dFree(f))
# define dFree7(a,b,c,d,e,f,g) (dFree(a) || dFree(b) || dFree(c) || dFree(d) || dFree(e) || dFree(f) || dFree(g))
#else
# define dMallocA2(n0,p0,n1,p1)                                         \
  (dMalloc((n0)*sizeof(**(p0))+(n1)*sizeof(**(p1))+(dDEFAULT_ALIGN-1),(p0)) \
   || (*(p1) = dNextAligned(*(p0)+(n0)),0))
# define dMallocA3(n0,p0,n1,p1,n2,p2)                                   \
  (dMalloc((n0)*sizeof(**(p0))+(n1)*sizeof(**(p1))+(n2)*sizeof(**(p2))+2*(dDEFAULT_ALIGN-1),(p0)) \
   || (*(p1) = dNextAligned(*(p0)+(n0)),*(p2) = dNextAligned(*(p1)+(n1)),0))
# define dMallocA4(n0,p0,n1,p1,n2,p2,n3,p3)                             \
  (dMalloc((n0)*sizeof(**(p0))+(n1)*sizeof(**(p1))+(n2)*sizeof(**(p2))+(n3)*sizeof(**(p3))+3*(dDEFAULT_ALIGN-1),(p0)) \
   || (*(p1) = dNextAligned(*(p0)+(n0)),*(p2) = dNextAligned(*(p1)+(n1)),*(p3) = dNextAligned(*(p2)+(n2)),0))
# define dMallocA5(n0,p0,n1,p1,n2,p2,n3,p3,n4,p4)                       \
  (dMalloc((n0)*sizeof(**(p0))+(n1)*sizeof(**(p1))+(n2)*sizeof(**(p2))+(n3)*sizeof(**(p3))+(n4)*sizeof(**(p4))+4*(dDEFAULT_ALIGN-1),(p0)) \
   || (*(p1) = dNextAligned(*(p0)+(n0)),*(p2) = dNextAligned(*(p1)+(n1)),*(p3) = dNextAligned(*(p2)+(n2)),*(p4) = dNextAligned(*(p3)+(n3)),0))
# define dMallocA6(n0,p0,n1,p1,n2,p2,n3,p3,n4,p4,n5,p5)                 \
  (dMalloc((n0)*sizeof(**(p0))+(n1)*sizeof(**(p1))+(n2)*sizeof(**(p2))+(n3)*sizeof(**(p3))+(n4)*sizeof(**(p4))+(n5)*sizeof(**(p5))+5*(dDEFAULT_ALIGN-1),(p0)) \
   || (*(p1) = dNextAligned(*(p0)+(n0)),*(p2) = dNextAligned(*(p1)+(n1)),*(p3) = dNextAligned(*(p2)+(n2)),*(p4) = dNextAligned(*(p3)+(n3)),*(p5) = dNextAligned(*(p4)+(n4)),0))
# define dMallocA7(n0,p0,n1,p1,n2,p2,n3,p3,n4,p4,n5,p5,n6,p6)             \
  (dMalloc((n0)*sizeof(**(p0))+(n1)*sizeof(**(p1))+(n2)*sizeof(**(p2))+(n3)*sizeof(**(p3))+(n4)*sizeof(**(p4))+(n5)*sizeof(**(p5))+(n6)*sizeof(**(p6))+6*(dDEFAULT_ALIGN-1),(p0)) \
   || (*(p1) = dNextAligned(*(p0)+(n0)),*(p2) = dNextAligned(*(p1)+(n1)),*(p3) = dNextAligned(*(p2)+(n2)),*(p4) = dNextAligned(*(p3)+(n3)),*(p5) = dNextAligned(*(p4)+(n4)),*(p6) = dNextAligned(*(p5)+(n5)),0))
# define dFree2(a,b) dFree(a)
# define dFree3(a,b,c) dFree(a)
# define dFree4(a,b,c,d) dFree(a)
# define dFree5(a,b,c,d,e) dFree(a)
# define dFree6(a,b,c,d,e,f) dFree(a)
# define dFree7(a,b,c,d,e,f,g) dFree(a)
#endif

#if defined(PETSC_USE_DEBUG)
# define dFunctionBegin                                                 \
  do {                                                                  \
    if (petscstack && (petscstack->currentsize < PETSCSTACKSIZE)) {     \
      petscstack->function[petscstack->currentsize]  = __func__;        \
      petscstack->file[petscstack->currentsize]      = __FILE__;        \
      petscstack->directory[petscstack->currentsize] = __SDIR__;        \
      petscstack->line[petscstack->currentsize]      = __LINE__;        \
      petscstack->currentsize++;                                        \
    }} while (0)
# define dFunctionReturn(a)                     \
  do {                                          \
    PetscStackPop;                              \
    return(a);} while (0)

# define dFunctionReturnVoid()                  \
  do {                                          \
    PetscStackPop;                              \
    return;} while (0)
#else
# define dFunctionBegin do { } while (0)
# define dFunctionReturn(a) return (a)
# define dFunctienReturnVoid() return
#endif

#define dASSERT(cond) if (!(cond)) { dERROR(1,"Assertion failed: " #cond); }

#endif
