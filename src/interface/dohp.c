#include "include/private/dohpimpl.h"
#include "src/inline/ilu.h"

PetscCookie DOHP_COOKIE;
PetscLogEvent DOHP_MatMult, DOHP_FunctionEval, DOHP_JacobianEval;

inline void DohpConvexComb_2_4(PetscReal x,PetscReal y,const PetscReal v[],const PetscInt p[],PetscReal f[])
{
  PetscInt i;
  for (i=0; i<3; i++) {
    f[i] = 0.25*((-1-x)*((-1-y)*v[p[0]*3+i] + (1+y)*v[p[3]*3+i]) + (1+x)*((-1-y)*v[p[1]*3+i] + (1+y)*v[p[2]*3+i]));
  }
}

/* Evaluate the Jacobian at quadrature points, this implementation simply uses the vertex
coordinates (placed in emap) to define the coordinate transformation.  Later implementations can use
higher order parametrically mapped elements or can optimize for affine mappings. */
#undef __FUNCT__
#define __FUNCT__ "DohpQuotientComputeElemJac_Hex"
PetscErrorCode DohpQuotientComputeElemJac_Hex(const DohpEMap_Hex *e,const DohpRule_Hex *rule,PetscReal *detJ,PetscReal *Jac,PetscReal *invJ)
{
  const DohpRule_Line *r=rule->l;
  PetscReal *J,f[6][3];
  PetscInt i,j,k,l,I;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (i=0; i<r[0].size; i++) {
    for (j=0; j<r[1].size; j++) {
      for (k=0; k<r[2].size; k++) {
        I = (i*r[1].size + j)*r[2].size + k; /* Index of the current quadrature point. */
        J = &Jac[I*9];                       /* To ease indexing while we build the Jacobian in the output array. */
        /* Compute coordinates at the projection of the quadrature point on each face. */
        DohpConvexComb_2_4(r[0].coord[i],r[2].coord[k],e->vtx,DohpHexQuad[0],f[0]);
        DohpConvexComb_2_4(r[1].coord[j],r[2].coord[k],e->vtx,DohpHexQuad[1],f[1]);
        DohpConvexComb_2_4(r[0].coord[i],r[2].coord[k],e->vtx,DohpHexQuad[2],f[2]);
        DohpConvexComb_2_4(r[1].coord[j],r[2].coord[k],e->vtx,DohpHexQuad[3],f[3]);
        DohpConvexComb_2_4(r[0].coord[i],r[1].coord[j],e->vtx,DohpHexQuad[4],f[4]);
        DohpConvexComb_2_4(r[0].coord[i],r[1].coord[j],e->vtx,DohpHexQuad[5],f[5]);
        for (l=0; l<3; l++) {
          J[0*3+l] = 0.5*(f[1][l] - f[3][l]);
          J[1*3+l] = 0.5*(f[2][l] - f[0][l]);
          J[2*3+l] = 0.5*(f[5][l] - f[4][l]);
        }
        detJ[I] = (J[0*3+0]*(J[1*3+1]*J[2*3+2] - J[1*3+2]*J[2*3+1]) +
                   J[0*3+1]*(J[1*3+2]*J[2*3+0] - J[1*3+0]*J[2*3+2]) +
                   J[0*3+2]*(J[1*3+0]*J[2*3+1] - J[1*3+1]*J[2*3+0]));
        ierr = PetscMemcpy(&invJ[I*9],J,9*sizeof(PetscReal));CHKERRQ(ierr);
        ierr = Kernel_A_gets_inverse_A_3(&invJ[I*9],0.0);CHKERRQ(ierr);
      }
    }
  }
  PetscFunctionReturn(0);
}
