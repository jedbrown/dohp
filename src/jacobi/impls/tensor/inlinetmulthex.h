#if !defined _INLINETMULTHEX_H
#define _INLINETMULTHEX_H

#include "dohptype.h"

/**
* The core computational kernel.  Performs a tensor product operation with the matrices A[0..2].
*
* @param[in] D number of degrees of freedom
* @param[in] P array of length 3, dimensions of input block
* @param[in] Q array of length 3, dimensions of output block
* @param[in] A array of pointers, length 3, transformation matrix in each direction, shape(A[i])=(Q[i],P[i]), i=0,1,2
* @param[in] f input vector with size corresponding to \in
* @param[in,out] g output vector with size corresponding to \out
* @param[in] imode ADD_VALUES or INSERT_VALUES
*
* @return err
*/
static dErr TensorMult_Hex_nounroll(dInt D,const dInt P[3],const dInt Q[3],const dReal *A[3],const dScalar in[],dScalar out[restrict],InsertMode imode)
{
  dScalar amem[Q[0]*P[1]*P[2]*D],bmem[Q[0]*Q[1]*P[2]*D];
  const dReal (*restrict Ax)[P[0]] = (const dReal (*)[P[0]])A[0];
  const dReal (*restrict Ay)[P[1]] = (const dReal (*)[P[1]])A[1];
  const dReal (*restrict Az)[P[2]] = (const dReal (*)[P[2]])A[2];
  dErr err;

  dFunctionBegin;
  err = dMemzero(amem,sizeof(amem));dCHK(err);
  err = dMemzero(bmem,sizeof(bmem));dCHK(err);
  switch (imode) {
    case INSERT_VALUES:
      err = dMemzero(out,Q[0]*Q[1]*Q[2]*D*sizeof(out[0]));dCHK(err);
      break;
    case ADD_VALUES:
      break;
    default:
      dERROR(1,"Requested InsertMode %d not supported for this operation.",imode);
  }

  for (dInt l=0; l<Q[0]; l++) {
    const dInt JKD = P[1]*P[2]*D;
    const dScalar *restrict Axl = Ax[l];
    dScalar *restrict aa = &amem[l*JKD];
    for (dInt i=0; i<P[0]; i++) {
      const dScalar *restrict ff = &in[i*JKD];
      const dReal aax = Axl[i];
      for (dInt jkd=0; jkd<JKD; jkd++) {
        aa[jkd] += aax * ff[jkd];
      }
    }
  }
  for (dInt i=0; i<Q[0]; i++) {
    const dInt JKD = P[1]*P[2]*D,KD = P[2]*D;
    dScalar *restrict aa = &amem[i*JKD];
    for (dInt l=0; l<Q[1]; l++) {
      dScalar *restrict bb = &bmem[(i*Q[1]+l)*KD];
      const dReal *restrict Ayl = Ay[l];
      for (dInt j=0; j<P[1]; j++) {
        const dScalar *restrict aaa = &aa[j*KD];
        const dReal aay = Ayl[j];
        for (dInt kd=0; kd<KD; kd++) {
          bb[kd] += aay * aaa[kd];
        }
      }
    }
  }
  switch (D) {
    case 1:
      for (dInt ij=0; ij<Q[0]*Q[1]; ij++) {
        const dScalar *restrict bb = &bmem[ij*P[2]];
        dScalar *restrict cc = &out[ij*Q[2]];
        for (dInt l=0; l<Q[2]; l++) {
          const dScalar *restrict aaz = Az[l];
          dScalar s = 0;
          for (dInt k=0; k<P[2]; k++) {
            s += aaz[k] * bb[k];
          }
          cc[l] += s;
        }
      }
      break;
    default:
      for (dInt ij=0; ij<Q[0]*Q[1]; ij++) {
        const dScalar *restrict bb = &bmem[ij*P[2]*D];
        for (dInt l=0; l<Q[2]; l++) {
          dScalar *restrict cc = &out[(ij*Q[2]+l)*D];
          const dScalar *restrict aa = Az[l];
          dInt kd = 0;
          for (dInt k=0; k<P[2]; k++) {
            const dScalar aaz = aa[k];
            for (dInt d=0; d<D; d++,kd++) {
              cc[d] += aaz * bb[kd];
            }
          }
        }
      }
  }
  PetscLogFlops((Q[0]*P[0]*P[1]*P[2] + Q[0]*Q[1]*P[1]*P[2] + Q[0]*Q[1]*Q[2]*P[2])*D*2);
  dFunctionReturn(0);
}


#if defined __SSE3__
#include <pmmintrin.h>
#include "private/microbench.h"


#define D 1
#define P2 4
#define Q2 4
static dErr TensorMult_Hex_P4_Q4_D1(dInt D_is_1,const dInt P[3],const dInt Q[3],const dReal *A[3],const dScalar in[],dScalar out[restrict],InsertMode imode)
{
  dScalar amem[Q[0]*P[1]*P2*D],bmem[Q[0]*Q[1]*P2*D];
  const dReal (*restrict Ax)[P[0]] = (const dReal (*)[P[0]])A[0];
  const dReal (*restrict Ay)[P[1]] = (const dReal (*)[P[1]])A[1];
  const dReal (*restrict Az)[P2] = (const dReal (*)[P2])A[2];
  dErr err;
/* #define VERBOSE_TIMING 1 */
#ifdef VERBOSE_TIMING
  uint32_t time0,time1,time2,time3;
#endif

  dFunctionBegin;
  if (P[2] != P2 || Q[2] != Q2 || D_is_1 != D) dERROR(1,"input sizes do not agree with unrolled sizes");
  err = dMemzero(amem,sizeof(amem));dCHK(err);
  _mm_prefetch((const char*)in,_MM_HINT_T0);
  err = dMemzero(bmem,sizeof(bmem));dCHK(err);
  switch (imode) {
    case INSERT_VALUES:
      err = dMemzero(out,Q[0]*Q[1]*Q2*D*sizeof(out[0]));dCHK(err);
      break;
    case ADD_VALUES:
      break;
    default:
      dERROR(1,"Requested InsertMode %d not supported for this operation.",imode);
  }

#ifdef VERBOSE_TIMING
  time0 = read_time();
#endif

#if 0
  for (dInt l=0; l<Q[0]; l++) {
    const dInt JKD = P[1]*P2*D;
    const dScalar *restrict Axl = Ax[l];
    dScalar *restrict aa = &amem[l*JKD];
    for (dInt i=0; i<P[0]; i++) {
      const dScalar *restrict ff = &in[i*JKD];
      const dReal aax = Axl[i];
      for (dInt jkd=0; jkd<JKD; jkd++) {
        aa[jkd] += aax * ff[jkd];
      }
    }
  }
#elif 0
  for (dInt l=0; l<Q[0]; l+=2) {
    const dInt JKD = P[1]*P2*D;
    dScalar *restrict aa0 = &amem[l*JKD], *restrict aa1 = &amem[(l+1)*JKD];
    for (dInt i=0; i<P[0]; i+=2) {
      const dScalar *restrict ff0 = &in[i*JKD], *restrict ff1 = &in[(i+1)*JKD];
      const dReal aax00 = Ax[l][i],aax01 = Ax[l][i+1],aax10 = Ax[l+1][i],aax11 = Ax[l+1][i+1];
      for (dInt j=0; j<JKD; j+=P2*D) {
        aa0[j+0] += aax00 * ff0[j+0] + aax01 * ff1[j+0];
        aa0[j+1] += aax00 * ff0[j+1] + aax01 * ff1[j+1];
        aa0[j+2] += aax00 * ff0[j+2] + aax01 * ff1[j+2];
        aa0[j+3] += aax00 * ff0[j+3] + aax01 * ff1[j+3];
        aa1[j+0] += aax10 * ff0[j+0] + aax11 * ff1[j+0];
        aa1[j+1] += aax10 * ff0[j+1] + aax11 * ff1[j+1];
        aa1[j+2] += aax10 * ff0[j+2] + aax11 * ff1[j+2];
        aa1[j+3] += aax10 * ff0[j+3] + aax11 * ff1[j+3];
      }
    }
  }
#else
  for (dInt l=0; l<Q[0]; l+=2) {
    const dInt JKD = P[1]*P2*D;
    dScalar *restrict aa0 = &amem[l*JKD], *restrict aa1 = &amem[(l+1)*JKD];
    for (dInt i=0; i<P[0]; i+=2) {
      const dScalar *restrict ff0 = &in[i*JKD], *restrict ff1 = &in[(i+1)*JKD];
      //const dReal aax00 = Ax[l][i],aax01 = Ax[l][i+1],aax10 = Ax[l+1][i],aax11 = Ax[l+1][i+1];
      __m128d
        aax00 = _mm_loaddup_pd(&Ax[l][i]),
        aax01 = _mm_loaddup_pd(&Ax[l][i+1]),
        aax10 = _mm_loaddup_pd(&Ax[l+1][i]),
        aax11 = _mm_loaddup_pd(&Ax[l+1][i+1]);
      for (dInt j=0; j<JKD; j+=P2*D) {
        __m128d
          ff00 = _mm_load_pd(&ff0[j]),ff02 = _mm_load_pd(&ff0[j+2]),
          ff10 = _mm_load_pd(&ff1[j]),ff12 = _mm_load_pd(&ff1[j+2]),
          aa00 = _mm_load_pd(&aa0[j]),aa02 = _mm_load_pd(&aa0[j+2]),
          aa10 = _mm_load_pd(&aa1[j]),aa12 = _mm_load_pd(&aa1[j+2]);
        aa00 = _mm_add_pd(aa00,_mm_add_pd(_mm_mul_pd(aax00,ff00),_mm_mul_pd(aax01,ff10)));
        aa02 = _mm_add_pd(aa02,_mm_add_pd(_mm_mul_pd(aax00,ff02),_mm_mul_pd(aax01,ff12)));
        aa10 = _mm_add_pd(aa10,_mm_add_pd(_mm_mul_pd(aax10,ff00),_mm_mul_pd(aax11,ff10)));
        aa12 = _mm_add_pd(aa12,_mm_add_pd(_mm_mul_pd(aax10,ff02),_mm_mul_pd(aax11,ff12)));
        _mm_store_pd(&aa0[j+0],aa00);
        _mm_store_pd(&aa0[j+2],aa02);
        _mm_store_pd(&aa1[j+0],aa10);
        _mm_store_pd(&aa1[j+2],aa12);
      }
    }
  }
#endif

#ifdef VERBOSE_TIMING
  time1 = read_time();
#endif

  _mm_prefetch((const char*)out,_MM_HINT_T1); /* Make sure it is at least in L2 when we get to it after this loop */

#if 0
  for (dInt i=0; i<Q[0]; i++) {
    const dInt JKD = P[1]*P2*D,KD = P2*D;
    dScalar (*restrict aa)[KD] = (dScalar(*)[4])&amem[i*JKD];
    dScalar (*restrict bb)[KD] = (dScalar(*)[4])&bmem[(i*Q[1])*KD];
    for (dInt l=0; l<Q[1]; l++) {
      for (dInt j=0; j<P[1]; j++) {
        for (dInt k=0; k<KD; k++) {
          bb[l][k] += Ay[l][j] * aa[j][k];
        }
      }
    }
  }
#elif 0
  if (P[1] != 4) dERROR(1,"not supported");
  for (dInt i=0; i<Q[0]; i++) {
    const dInt JKD = P[1]*P2*D,KD = P2*D;
    dScalar (*restrict aa)[KD] = (dScalar(*)[4])&amem[i*JKD];
    dScalar (*restrict bb)[KD] = (dScalar(*)[4])&bmem[(i*Q[1])*KD];
    for (dInt l=0; l<Q[1]; l++) { /* \a k is unrolled down, \a j is unrolled across */
      bb[l][0] = Ay[l][0]*aa[0][0] + Ay[l][1]*aa[1][0] + Ay[l][2]*aa[2][0] + Ay[l][3]*aa[3][0];
      bb[l][1] = Ay[l][0]*aa[0][1] + Ay[l][1]*aa[1][1] + Ay[l][2]*aa[2][1] + Ay[l][3]*aa[3][1];
      bb[l][2] = Ay[l][0]*aa[0][2] + Ay[l][1]*aa[1][2] + Ay[l][2]*aa[2][2] + Ay[l][3]*aa[3][2];
      bb[l][3] = Ay[l][0]*aa[0][3] + Ay[l][1]*aa[1][3] + Ay[l][2]*aa[2][3] + Ay[l][3]*aa[3][3];
    }
  }
#elif 0
  for (dInt i=0; i<Q[0]; i++) {
    const dInt JKD = P[1]*P2*D,KD = P2*D;
    dScalar *restrict aa = &amem[i*JKD];
    for (dInt l=0; l<Q[1]; l++) {
      dScalar *restrict bb = &bmem[(i*Q[1]+l)*KD];
      const dReal *restrict Ayl = Ay[l];
      for (dInt j=0; j<P[1]; j++) {
        const dScalar *restrict aaa = &aa[j*KD];
        const dReal aay = Ayl[j];
        for (dInt kd=0; kd<KD; kd++) {
          bb[kd] += aay * aaa[kd];
        }
      }
    }
  }
#elif 0
  for (dInt i=0; i<Q[0]; i++) {
    const dInt JKD = P[1]*P2*D,KD = P2*D;
    dScalar *restrict aa = &amem[i*JKD];
    for (dInt l=0; l<Q[1]; l+=2) {
      dScalar *restrict bb0 = &bmem[(i*Q[1]+l)*KD],*restrict bb1 = &bmem[(i*Q[1]+l+1)*KD];
      for (dInt j=0; j<P[1]; j+=2) {
        const dScalar *restrict aaa0 = &aa[j*KD],*restrict aaa1 = &aa[(j+1)*KD];
        const dReal aay00 = Ay[l][j],aay01 = Ay[l][j+1],aay10 = Ay[l+1][j],aay11 = Ay[l+1][j+1];
        for (dInt kd=0; kd<KD; kd+=4) {
          bb0[0] += aay00 * aaa0[0] + aay01 * aaa1[0];
          bb0[1] += aay00 * aaa0[1] + aay01 * aaa1[1];
          bb0[2] += aay00 * aaa0[2] + aay01 * aaa1[2];
          bb0[3] += aay00 * aaa0[3] + aay01 * aaa1[3];
          bb1[0] += aay10 * aaa0[0] + aay11 * aaa1[0];
          bb1[1] += aay10 * aaa0[1] + aay11 * aaa1[1];
          bb1[2] += aay10 * aaa0[2] + aay11 * aaa1[2];
          bb1[3] += aay10 * aaa0[3] + aay11 * aaa1[3];
        }
      }
    }
  }
#elif 0
  for (dInt i=0; i<Q[0]; i++) {
    const dInt JKD = P[1]*P2*D,KD = P2*D;
    dScalar *restrict aa = &amem[i*JKD];
    for (dInt l=0; l<Q[1]; l+=2) {
      dScalar *restrict bb0 = &bmem[(i*Q[1]+l)*KD],*restrict bb1 = &bmem[(i*Q[1]+l+1)*KD];
      //if ((size_t)bb0 & 0xf || (size_t)bb1 & 0xf) dERROR(1,"bmem is unaligned");
      __m128d bb00 = _mm_setzero_pd(),bb02=bb00,bb10=bb00,bb12=bb00;
      for (dInt j=0; j<P[1]; j+=2) {
        const dScalar *restrict aaa0 = &aa[j*KD],*restrict aaa1 = &aa[(j+1)*KD];
        //if ((size_t)aaa0 & 0xf || (size_t)aaa1 & 0xf) dERROR(1,"aa is unaligned");
        __m128d
          aaa00 = _mm_load_pd(aaa0),aaa02 = _mm_load_pd(aaa0+2),aaa10 = _mm_load_pd(aaa1),aaa12 = _mm_load_pd(aaa1+2),
          aay00 = _mm_loaddup_pd(&Ay[l][j]),aay01 = _mm_loaddup_pd(&Ay[l][j+1]), aay10 = _mm_loaddup_pd(&Ay[l+1][j]),aay11 = _mm_loaddup_pd(&Ay[l+1][j+1]);
        bb00 = _mm_add_pd(bb00,_mm_add_pd(_mm_mul_pd(aay00,aaa00),_mm_mul_pd(aay01,aaa10)));
        bb02 = _mm_add_pd(bb02,_mm_add_pd(_mm_mul_pd(aay00,aaa02),_mm_mul_pd(aay01,aaa12)));
        bb10 = _mm_add_pd(bb10,_mm_add_pd(_mm_mul_pd(aay10,aaa00),_mm_mul_pd(aay11,aaa10)));
        bb12 = _mm_add_pd(bb12,_mm_add_pd(_mm_mul_pd(aay10,aaa02),_mm_mul_pd(aay11,aaa12)));
      }
      _mm_store_pd(bb0,bb00);
      _mm_store_pd(bb0+2,bb02);
      _mm_store_pd(bb1,bb10);
      _mm_store_pd(bb1+2,bb12);
    }
  }
#else
  if (P[1] != P2) dERROR(1,"not supported");
  for (dInt i=0; i<Q[0]; i++) {
    const dInt JKD = P[1]*P2*D,KD = P2*D;
    dScalar (*restrict aa)[KD] = (dScalar(*)[KD])&amem[i*JKD];
    dScalar (*restrict bb)[KD] = (dScalar(*)[KD])&bmem[(i*Q[1])*KD];
    __m128d
      aa00 = _mm_load_pd(&aa[0][0]),aa02 = _mm_load_pd(&aa[0][2]),aa10 = _mm_load_pd(&aa[1][0]),aa12 = _mm_load_pd(&aa[1][2]),
      aa20 = _mm_load_pd(&aa[2][0]),aa22 = _mm_load_pd(&aa[2][2]),aa30 = _mm_load_pd(&aa[3][0]),aa32 = _mm_load_pd(&aa[3][2]);
    for (dInt l=0; l<Q[1]; l++) {
      __m128d ayl0p = _mm_loaddup_pd(&Ay[l][0]),ayl1p = _mm_loaddup_pd(&Ay[l][1]),ayl2p = _mm_loaddup_pd(&Ay[l][2]),ayl3p = _mm_loaddup_pd(&Ay[l][3]);
      __m128d t0,t2;
      t0 = _mm_add_pd(_mm_add_pd(_mm_mul_pd(ayl0p,aa00),_mm_mul_pd(ayl1p,aa10)),  _mm_add_pd(_mm_mul_pd(ayl2p,aa20),_mm_mul_pd(ayl3p,aa30)));
      t2 = _mm_add_pd(_mm_add_pd(_mm_mul_pd(ayl0p,aa02),_mm_mul_pd(ayl1p,aa12)),  _mm_add_pd(_mm_mul_pd(ayl2p,aa22),_mm_mul_pd(ayl3p,aa32)));
      _mm_store_pd(&bb[l][0],t0);
      _mm_store_pd(&bb[l][2],t2);
    }
  }
#endif

#ifdef VERBOSE_TIMING
  time2 = read_time();
#endif

#if 0
  for (dInt ij=0; ij<Q[0]*Q[1]; ij++) {
    const dScalar *restrict bb = &bmem[ij*P2*D];
    for (dInt l=0; l<Q2; l++) {
      dScalar *restrict cc = &out[(ij*Q2+l)*D];
      const dScalar *restrict aa = Az[l];
      dInt kd = 0;
      for (dInt k=0; k<P2; k++) {
        const dScalar aaz = aa[k];
        for (dInt d=0; d<D; d++,kd++) {
          cc[d] += aaz * bb[kd];
        }
      }
    }
  }
#elif 0
  for (dInt ij=0; ij<Q[0]*Q[1]; ij++) {
    const dScalar *restrict bb = &bmem[ij*P2*D];
    dScalar *restrict cc = &out[(ij*Q2)*D];
    for (dInt l=0; l<Q2; l++) {
      cc[l] += Az[l][0]*bb[0] + Az[l][1]*bb[1] + Az[l][2]*bb[2] + Az[l][3]*bb[3];
      //for (dInt k=0; k<P2; k++) cc[l] += Az[l][k] * bb[k];
    }
  }
#elif 0
  for (dInt ij=0; ij<Q[0]*Q[1]; ij++) {
    const dScalar               /* This should be hoisted out, check assembly to confirm */
      aa00 = Az[0][0],aa01 = Az[0][1],aa02 = Az[0][2],aa03 = Az[0][3],
      aa10 = Az[1][0],aa11 = Az[1][1],aa12 = Az[1][2],aa13 = Az[1][3],
      aa20 = Az[2][0],aa21 = Az[2][1],aa22 = Az[2][2],aa23 = Az[2][3],
      aa30 = Az[3][0],aa31 = Az[3][1],aa32 = Az[3][2],aa33 = Az[3][3];
    const dScalar *restrict bb = &bmem[ij*P2*D],bb0=bb[0],bb1=bb[1],bb2=bb[2],bb3=bb[3];
    dScalar *restrict cc = &out[(ij*Q2)*D];
    cc[0] += aa00*bb0 + aa01*bb1 + aa02*bb2 + aa03*bb3;
    cc[1] += aa10*bb0 + aa11*bb1 + aa12*bb2 + aa13*bb3;
    cc[2] += aa20*bb0 + aa21*bb1 + aa22*bb2 + aa23*bb3;
    cc[3] += aa30*bb0 + aa31*bb1 + aa32*bb2 + aa33*bb3;
  }
#else
  {
    __m128d
      aa00 = _mm_load_pd(&Az[0][0]),aa02 = _mm_load_pd(&Az[0][2]),
      aa10 = _mm_load_pd(&Az[1][0]),aa12 = _mm_load_pd(&Az[1][2]),
      aa20 = _mm_load_pd(&Az[2][0]),aa22 = _mm_load_pd(&Az[2][2]),
      aa30 = _mm_load_pd(&Az[3][0]),aa32 = _mm_load_pd(&Az[3][2]);
    for (dInt ij=0; ij<Q[0]*Q[1]; ij++) {
      const dScalar *restrict bb = &bmem[ij*P2*D];
      dScalar *restrict cc = &out[ij*Q2*D];
      __m128d t0,t1,t2,t3,
        bb0 = _mm_load_pd(bb),bb2 = _mm_load_pd(bb+2),
        cc0 = _mm_load_pd(cc),cc2 = _mm_load_pd(cc+2);
      t0 = _mm_add_pd(_mm_mul_pd(aa00,bb0),_mm_mul_pd(aa02,bb2));
      t1 = _mm_add_pd(_mm_mul_pd(aa10,bb0),_mm_mul_pd(aa12,bb2));
      t2 = _mm_add_pd(_mm_mul_pd(aa20,bb0),_mm_mul_pd(aa22,bb2));
      t3 = _mm_add_pd(_mm_mul_pd(aa30,bb0),_mm_mul_pd(aa32,bb2));
      _mm_store_pd(cc+0,_mm_add_pd(cc0,_mm_hadd_pd(t0,t1)));
      _mm_store_pd(cc+2,_mm_add_pd(cc2,_mm_hadd_pd(t2,t3)));
    }
  }
#endif

#ifdef VERBOSE_TIMING
  time3 = read_time();
  printf("%s cycle counts %10ld %10ld %10ld\n",__func__,(long)(time1-time0),(long)(time2-time1),(long)(time3-time2));
#endif

  PetscLogFlops((Q[0]*P[0]*P[1]*P2 + Q[0]*Q[1]*P[1]*P2 + Q[0]*Q[1]*Q2*P2)*D*2);
  dFunctionReturn(0);
}

#endif
#endif
