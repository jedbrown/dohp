#ifndef __DOHPBASE_H
#define __DOHPBASE_H

PETSC_EXTERN_CXX_BEGIN

#define INLINE static inline

INLINE void GeomVecPlus(const PetscReal a[],const PetscReal b[],PetscReal c[]) {PetscInt i; for (i = 0; i<3; i++) c[i] = a[i]+b[i];}
INLINE void GeomVecMinus(const PetscReal a[],const PetscReal b[],PetscReal c[]) {PetscInt i; for (i = 0; i<3; i++) c[i] = a[i]-b[i];}
INLINE PetscReal GeomDotProd(const PetscReal a[],const PetscReal b[]) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
/* c = a \cross b */
INLINE void GeomCrossProd(const PetscReal a[],const PetscReal b[],PetscReal c[])
{c[0]=a[1]*b[2]-a[2]*b[1]; c[1]=a[2]*b[0]-a[0]*b[2]; c[2]=a[0]*b[1]+a[1]*b[0];}
/* coordinates for vertices are interlaced in a, result goes in b */
INLINE void GeomVecMeanI(const PetscInt n,const PetscReal a[],PetscReal b[])
{PetscInt i,j; b[0]=b[1]=b[2]=0; for (i=0;i<n;i++) {for (j=0;j<3;j++) b[j]+=a[3*i+j]/n;}}
/* coordinates are interlaced in a, result in b */
INLINE void GeomQuadFaceNormal(const PetscReal a[],PetscReal b[]) {
  PetscReal f[3],g[3],h[3];
  GeomVecMinus(a+3,a,b); GeomVecMinus(a+6,a+9,f); GeomVecPlus(b,f,g); /* g = mean vector in x direction */
  GeomVecMinus(a+9,a,b); GeomVecMinus(a+6,a+3,f); GeomVecPlus(b,f,h); /* h = mean vector in y direction */
  GeomCrossProd(g,h,b);
}
INLINE PetscTruth GeomQuadParallel(const PetscReal a[],const PetscReal b[]) /* return true if Quad a \dot b > 0 */
{PetscReal f[3]; GeomQuadFaceNormal(a,f); return (GeomDotProd(b,f) > 0);}

INLINE void GeomConvexComb_2_4(PetscReal x,PetscReal y,const PetscReal v[],const PetscInt p[],PetscReal f[])
{
  PetscInt i;
  for (i=0; i<3; i++) {
    f[i] = 0.25*((-1-x)*((-1-y)*v[p[0]*3+i] + (1+y)*v[p[3]*3+i]) + (1+x)*((-1-y)*v[p[1]*3+i] + (1+y)*v[p[2]*3+i]));
  }
}


PETSC_EXTERN_CXX_END
#endif  /* __DOHPBASE_H */
