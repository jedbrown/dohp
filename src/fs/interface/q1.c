#include "dohpfs.h"
#include "dohpgeom.h"
#include "private/fsimpl.h"     /* Only for dLOG_Q1HexComputeQuadrature!  Otherwise this is stand-alone. */

/** Compute basis functions and quadrature weights for Q1 element
*
* @param[in] x global coordinates of the 8 corners, canonical ordering
* @param[out] n number of quadrature points
* @param[out] inqx will point to coordinates of quadrature points
* @param[out] injw will point to quadrature weights times jacobian determinant (this should be used to weight integration)
* @param[out] inbasis will point to array of basis functions evaluated at quadrature points \c basis[i*8+j] is basis function \a j evaluated at quadrature point \a i
* @param[out] inderiv will point to array of derivatives evaluated at quadrature points \c deriv[(i*8+j)*3+k] is derivative in \a k direction of basis \a j evaluated at quadrature point \a i
**/
dErr dQ1HexComputeQuadrature(const dReal x[8][3],dInt *n,const dReal (**inqx)[3],const dReal **injw,const dReal **inbasis,const dReal **inderiv)
{
  static const dReal linecoords[2] = {-0.57735026918962562,0.57735026918962562}; /* 2-point Gauss quadrature */
  //static const dReal lineweight[2] = {1.0,1.0};
  //static const dReal linebasis[2][2] = {{0.78867513459481287,0.21132486540518708},{0.21132486540518708,0.78867513459481287}};
  //static const dReal linederiv[2][2] = {{-0.5,-0.5},{0.5,0.5}};
  static const dReal hexweight[8] = {1,1,1,1, 1,1,1,1};
  static dReal qx[2*2*2][3],jw[2*2*2],basis[8][2*2*2],deriv[8][2*2*2][3];
  dErr err;

  dFunctionBegin;
  err = PetscLogEventBegin(dLOG_Q1HexComputeQuadrature,0,0,0,0);dCHK(err);
  for (dInt i=0; i<2; i++) {
    const dReal q0 = linecoords[i],q0m = 0.125*(1-q0),q0p = 0.125*(1+q0),qmd = q0m,qpd = q0p;
    for (dInt j=0; j<2; j++) {
      const dReal q1 = linecoords[j],q1m = 1-q1,q1p = 1+q1;
      const dReal qmm = q0m*q1m;
      const dReal qpm = q0p*q1m;
      const dReal qmp = q0m*q1p;
      const dReal qpp = q0p*q1p;
      const dReal qdm = 0.125*q1m,qdp = 0.125*q1p;
      for (dInt k=0; k<2; k++) {
        const dInt p = (i*2+j)*2+k;
        const dReal q2 = linecoords[k],q2m = 1-q2,q2p = 1+q2; /* location of quadrature point in reference coordinates */
        const dReal qmmm = qmm*q2m,qmmp = qmm*q2p,qmpm = qmp*q2m,qmpp = qmp*q2p;
        const dReal qpmm = qpm*q2m,qpmp = qpm*q2p,qppm = qpp*q2m,qppp = qpp*q2p;
        const dReal qdmm = qdm*q2m,qdmp = qdm*q2p,qdpm = qdp*q2m,qdpp = qdp*q2p;
        const dReal qmdm = qmd*q2m,qmdp = qmd*q2p,qpdm = qpd*q2m,qpdp = qpd*q2p;
        const dReal hbasis[8] = {qmmm,qpmm,qppm,qmpm,qmmp,qpmp,qppp,qmpp};         /* Hex basis at this quadrature point */
        const dReal hderiv[3][8] = {{-qdmm,qdmm,qdpm,-qdpm,-qdmp,qdmp,qdpp,-qdpp}, /* Hex reference deriv at this quadrature point */
                                    {-qmdm,-qpdm,qpdm,qmdm,-qmdp,-qpdp,qpdp,qmdp},
                                    {-qmm,-qpm,-qpp,-qmp,qmm,qpm,qpp,qmp}};
        dReal J[3][3],Jinv[3][3],Jdet;
        for (dInt l=0; l<3; l++) {
          /* Set the global coordinates at each quadrature point */
          qx[p][l] = (+ x[0][l]*qmmm + x[1][l]*qpmm + x[2][l]*qppm + x[3][l]*qmpm + x[4][l]*qmmp + x[5][l]*qpmp + x[6][l]*qppp + x[7][l]*qmpp);
          /* Compute the sub-element Jacobian at every quadrature point */
          J[l][0] =  (- x[0][l]*qdmm + x[1][l]*qdmm + x[2][l]*qdpm - x[3][l]*qdpm - x[4][l]*qdmp + x[5][l]*qdmp + x[6][l]*qdpp - x[7][l]*qdpp);
          J[l][1] =  (- x[0][l]*qmdm - x[1][l]*qpdm + x[2][l]*qpdm + x[3][l]*qmdm - x[4][l]*qmdp - x[5][l]*qpdp + x[6][l]*qpdp + x[7][l]*qmdp);
          J[l][2] =  (- x[0][l]*qmm - x[1][l]*qpm - x[2][l]*qpp - x[3][l]*qmp + x[4][l]*qmm + x[5][l]*qpm + x[6][l]*qpp + x[7][l]*qmp);
        }
        err = dGeomInvert3(&J[0][0],&Jinv[0][0],&Jdet);dCHK(err);
        jw[p] = hexweight[p]*Jdet;
        for (dInt l=0; l<8; l++) { /* Loop over corners in canonical order */
          basis[p][l] = hbasis[l]; /* Basis from corner \a l evaluated at this quadrature point (\a p) */
          /* Derivatives of basis from corner \a l with respect to global coordinates, evaluated at this quadrature point (\a p) */
          deriv[p][l][0] = hderiv[0][l] * Jinv[0][0] + hderiv[1][l] * Jinv[1][0] + hderiv[2][l] * Jinv[2][0];
          deriv[p][l][1] = hderiv[0][l] * Jinv[0][1] + hderiv[1][l] * Jinv[1][1] + hderiv[2][l] * Jinv[2][1];
          deriv[p][l][2] = hderiv[0][l] * Jinv[0][2] + hderiv[1][l] * Jinv[1][2] + hderiv[2][l] * Jinv[2][2];
        }
      }
    }
  }
  err = PetscLogFlops(2*(4+2*(8+2*(2+16+3*60+42+1+3*5))));dCHK(err);
  *n = 2*2*2;
  *inqx = (const dReal(*)[3])&qx[0][0];
  *injw = jw;
  *inbasis = &basis[0][0];
  *inderiv = &deriv[0][0][0];
  err = PetscLogEventEnd(dLOG_Q1HexComputeQuadrature,0,0,0,0);dCHK(err);
  dFunctionReturn(0);
}
