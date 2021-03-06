#if !defined _DOHP_H
#define _DOHP_H

#include "dohptype.h"
#include "dohputil.h"
#include <stdint.h>
#include <petscsys.h>

#define dPragmaQuote(a) _Pragma(#a)
#if defined dHAVE_PRAGMA_GCC
#  define dPragmaGCC(a) dPragmaQuote(GCC a)
#else
#  define dPragmaGCC(a)
#endif

#if defined dUSE_VALGRIND
#  include <valgrind/memcheck.h>
#  define dMakeMemUndefined(mem,bytes) do {             \
    memset(mem,0xff,bytes);                             \
    dPragmaGCC(diagnostic ignored "-Wconversion")       \
      VALGRIND_MAKE_MEM_UNDEFINED((mem),(bytes));       \
  } while (0)
#  define dCheckMemIsDefined(mem,bytes) do {            \
    dPragmaGCC(diagnostic ignored "-Wconversion")       \
      VALGRIND_CHECK_MEM_IS_DEFINED((mem),(bytes));     \
  } while (0)
#else
#  define dMakeMemUndefined(mem,bytes) do {     \
    if (dMemzero((mem),(bytes))) dERROR(PETSC_COMM_SELF,PETSC_ERR_MEMC,"Memory " #mem " is corrupt"); \
  } while (0)
#  define dCheckMemIsDefined(mem,bytes) do { } while (0)
#endif

#define dSTATUS_UNOWNED   (dEntStatus)0x1
#define dSTATUS_SHARED    (dEntStatus)0x2
#define dSTATUS_INTERFACE (dEntStatus)0x4
#define dSTATUS_GHOST     (dEntStatus)0x8

#define dDataTypeToITAPS(dtype,itype) (*(itype) = (dtype), 0)
#define dEntTopoToITAPS(dtopo,itopo) (*(itopo) = (dtopo), 0)
#define dEntTypeToITAPS(dtype,itype) (*(itype) = (dtype), 0)

#define dCHK(err) do {if (PetscUnlikely(err)) return PetscError(PETSC_COMM_SELF,__LINE__,__func__,__FILE__,__SDIR__,(err),PETSC_ERROR_REPEAT," ");} while (0)
#define dERROR(comm,n,...) return PetscError((comm),__LINE__,__func__,__FILE__,__SDIR__,(n),PETSC_ERROR_INITIAL,__VA_ARGS__)

#define dPrintf PetscPrintf
#define dMemcpy(a,b,c) PetscMemcpy(a,b,c)
#define dMemzero(a,b)  PetscMemzero(a,b)

#define dMalloc(a,b) PetscMalloc(a,b)
#define dNew(a,b) PetscNew(a,b)
#define dNewLog(a,b,c) PetscNewLog(a,b,c)
#define dFree(a) PetscFree(a)
#define dNewM(n,t,p) (dMalloc((n)*sizeof(t),(p)) || dMemzero(*(p),(n)*sizeof(t)))
#define dMallocM(n,t,p) (dMalloc((n)*sizeof(t),(p)))
#define dMallocA(n,p) (dMalloc((n)*sizeof(**(p)),(p)))
#define dCalloc(n,p) (dMalloc((n),(p)) || dMemzero(*(p),(n)))
#define dCallocA(n,p) (dCalloc((n)*sizeof(**(p)),(p)))

#define dNonNullElse(a,b) ((a)?(a):(b))
#define dStrlen(s,l) PetscStrlen((s),(l))
extern dErr dObjectGetComm(dObject obj,MPI_Comm *comm);

static inline dInt dMaxInt(dInt a,dInt b) { return (a > b) ? a : b; }
static inline dInt dMinInt(dInt a,dInt b) { return (a < b) ? a : b; }
static inline dReal dMax(dReal a,dReal b) { return (a > b) ? a : b; }
static inline dReal dMin(dReal a,dReal b) { return (a < b) ? a : b; }
static inline dReal dAbs(dScalar a) { return fabs(a); }
static inline dScalar dSqr(dScalar a) { return a * a; }
static inline dInt dSqrInt(dInt a) { return a * a; }
static inline dInt dSign(dReal a) { return a < 0 ? -1 : (a > 0); }
static inline dReal dSqrt(dReal a) { return sqrt(a); }
static inline dReal dPowReal(dReal a,dReal p) { return pow(a,p); }
static inline dReal dNormScalar3p(const dScalar a[3],const dReal p) { return dPowReal(dPowReal(dAbs(a[0]),p) + dPowReal(dAbs(a[1]),p) + dPowReal(dAbs(a[2]),p), 1./p); }
static inline dScalar dDotScalar3(const dScalar a[3],const dScalar b[3]) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static inline dScalar dDotScalarColumn3(const dScalar a[3],const dScalar b[9],dInt i) { return a[0]*b[0*3+i] + a[1]*b[1*3+i] + a[2]*b[2*3+i]; } // x.T * A[:,i]
static inline dScalar dColonSymScalar3(const dScalar a[6],const dScalar b[6])
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + 2*a[3]*b[3] + 2*a[4]*b[4] + 2*a[5]*b[5]; }
static inline void dTensorSymCompress3(const dScalar Du[9],dScalar Dv[6])
{Dv[0] = Du[0]; Dv[1] = Du[4]; Dv[2] = Du[8]; Dv[3] = (Du[1]+Du[3])/2; Dv[4] = (Du[2]+Du[6])/2; Dv[5] = (Du[5]+Du[7])/2;}
static inline void dTensorSymUncompress3(const dScalar Du[6],dScalar Dv[9])
{Dv[0] = Du[0]; Dv[4] = Du[1]; Dv[8] = Du[2]; Dv[1] = Dv[3] = Du[3]; Dv[2] = Dv[6] = Du[4]; Dv[5] = Dv[7] = Du[5];}
static inline void dTensorMultGESY3(dScalar Cf[9],const dScalar Af[9],const dScalar S[6])
{                               /* C = A*S, S has compressed symmetric storage */
  dScalar       (*C)[3] =       (dScalar(*)[3])Cf;
  const dScalar (*A)[3] = (const dScalar(*)[3])Af;
  for (dInt i=0; i<3; i++) {
    C[i][0] = A[i][0]*S[0] + A[i][1]*S[3] + A[i][2]*S[4];
    C[i][1] = A[i][0]*S[3] + A[i][1]*S[1] + A[i][2]*S[5];
    C[i][2] = A[i][0]*S[4] + A[i][1]*S[5] + A[i][2]*S[2];
  }
}
static inline void dTensorMultAddGESY3(dScalar Cf[9],const dScalar Af[9],const dScalar S[6])
{                               /* C += A*S, S has compressed symmetric storage */
  dScalar       (*C)[3] =       (dScalar(*)[3])Cf;
  const dScalar (*A)[3] = (const dScalar(*)[3])Af;
  for (dInt i=0; i<3; i++) {
    C[i][0] += A[i][0]*S[0] + A[i][1]*S[3] + A[i][2]*S[4];
    C[i][1] += A[i][0]*S[3] + A[i][1]*S[1] + A[i][2]*S[5];
    C[i][2] += A[i][0]*S[4] + A[i][1]*S[5] + A[i][2]*S[2];
  }
}

#define dGamma(a) tgamma(a) /* This is defined in math.h as of C99. */

/* stdbool.h has small (1 byte) bools, PETSc uses an enum which has few size guarantees, so we use it directly and keep
* it out of our public interface (which is why these typedefs are here and not in dohptype.h).
**/
typedef PetscBool dBool;
#define dTRUE  PETSC_TRUE
#define dFALSE PETSC_FALSE

#define dMAX_PATH_LEN PETSC_MAX_PATH_LEN
#define dNAME_LEN     256
#define dSTR_LEN      256

#define dDEFAULT_ALIGN PETSC_MEMALIGN

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
/* This is like PetscFunctionBegin, but does not check __FUNCT__ because we have C99 */
#  define dFunctionBegin                                                \
  do {                                                                  \
    PetscStack* petscstackp;                                            \
    petscstackp = (PetscStack*)PetscThreadLocalGetValue(petscstack);    \
    if (petscstackp && (petscstackp->currentsize < PETSCSTACKSIZE)) {   \
      petscstackp->function[petscstackp->currentsize]  = PETSC_FUNCTION_NAME; \
      petscstackp->file[petscstackp->currentsize]      = __FILE__;      \
      petscstackp->directory[petscstackp->currentsize] = __SDIR__;      \
      petscstackp->line[petscstackp->currentsize]      = __LINE__;      \
      petscstackp->petscroutine[petscstackp->currentsize] = PETSC_TRUE; \
      petscstackp->currentsize++;                                       \
    }                                                                   \
  } while (0)
#else
#  define dFunctionBegin do {} while(0)
#endif

#define dFunctionReturn(a) PetscFunctionReturn(a)

#define dASSERT(cond) if (!(cond)) { dERROR(PETSC_COMM_SELF,1,"Assertion failed: " #cond); }

#endif
