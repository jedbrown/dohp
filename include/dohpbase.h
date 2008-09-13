#ifndef _DOHPBASE_H
#define _DOHPBASE_H

#include "dohptype.h"

PETSC_EXTERN_CXX_BEGIN

#define INLINE static inline

INLINE void dGeomVecPlus(const dReal a[],const dReal b[],dReal c[]) {dInt i; for (i = 0; i<3; i++) c[i] = a[i]+b[i];}
INLINE void dGeomVecMinus(const dReal a[],const dReal b[],dReal c[]) {dInt i; for (i = 0; i<3; i++) c[i] = a[i]-b[i];}
INLINE dReal dGeomDotProd(const dReal a[],const dReal b[]) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
/* c = a \cross b */
INLINE void dGeomCrossProd(const dReal a[],const dReal b[],dReal c[])
{c[0]=a[1]*b[2]-a[2]*b[1]; c[1]=a[2]*b[0]-a[0]*b[2]; c[2]=a[0]*b[1]+a[1]*b[0];}
/* coordinates for vertices are interlaced in a, result goes in b */
INLINE void dGeomVecMeanI(const dInt n,const dReal a[],dReal b[])
{dInt i,j; b[0]=b[1]=b[2]=0; for (i=0;i<n;i++) {for (j=0;j<3;j++) b[j]+=a[3*i+j]/n;}}
/* coordinates are interlaced in a, result in b */
INLINE void dGeomQuadFaceNormal(const dReal a[],dReal b[]) {
  dReal f[3],g[3],h[3];
  dGeomVecMinus(a+3,a,b); dGeomVecMinus(a+6,a+9,f); dGeomVecPlus(b,f,g); /* g = mean vector in x direction */
  dGeomVecMinus(a+9,a,b); dGeomVecMinus(a+6,a+3,f); dGeomVecPlus(b,f,h); /* h = mean vector in y direction */
  dGeomCrossProd(g,h,b);
}
INLINE dBool dGeomQuadParallel(const dReal a[],const dReal b[]) /* return true if Quad a \dot b > 0 */
{dReal f[3]; dGeomQuadFaceNormal(a,f); return (dGeomDotProd(b,f) > 0);}

INLINE void dGeomConvexComb_2_4(dReal x,dReal y,const dReal v[],const dInt p[],dReal f[])
{
  dInt i;
  for (i=0; i<3; i++) {
    f[i] = 0.25*((-1-x)*((-1-y)*v[p[0]*3+i] + (1+y)*v[p[3]*3+i]) + (1+x)*((-1-y)*v[p[1]*3+i] + (1+y)*v[p[2]*3+i]));
  }
}



PETSC_EXTERN_CXX_END
#endif  /* _DOHPBASE_H */
