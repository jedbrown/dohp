#include "include/private/dohpimpl.h"
#include "src/inline/ilu.h"

#if 0
/* Evaluate the Jacobian at quadrature points, this implementation simply uses the vertex
coordinates (placed in emap) to define the coordinate transformation.  Later implementations can use
higher order parametrically mapped elements or can optimize for affine mappings. */
#undef __FUNCT__
#define __FUNCT__ "dQuotientComputeElemJac_Hex"
dErr dQuotientComputeElemJac_Hex(const DohpEMap *e,const dRule *rule,dReal detJ[],dReal Jac[],dReal *invJ)
{
  const dRule_Line *r=rule->l;
  dReal *J,f[6][3];
  dInt i,j,k,l,I;
  dErr err;

  dFunctionBegin;

  for (i=0; i<r[0].size; i++) {
    for (j=0; j<r[1].size; j++) {
      for (k=0; k<r[2].size; k++) {
        I = (i*r[1].size + j)*r[2].size + k; /* Index of the current quadrature point. */
        J = &Jac[I*9];                       /* To ease indexing while we build the Jacobian in the output array. */
        /* Compute coordinates at the projection of the quadrature point on each face. */
        dGeomConvexComb_2_4(r[0].coord[i],r[2].coord[k],e->vtx,DohpHexQuad[0],f[0]);
        dGeomConvexComb_2_4(r[1].coord[j],r[2].coord[k],e->vtx,DohpHexQuad[1],f[1]);
        dGeomConvexComb_2_4(r[0].coord[i],r[2].coord[k],e->vtx,DohpHexQuad[2],f[2]);
        dGeomConvexComb_2_4(r[1].coord[j],r[2].coord[k],e->vtx,DohpHexQuad[3],f[3]);
        dGeomConvexComb_2_4(r[0].coord[i],r[1].coord[j],e->vtx,DohpHexQuad[4],f[4]);
        dGeomConvexComb_2_4(r[0].coord[i],r[1].coord[j],e->vtx,DohpHexQuad[5],f[5]);
        for (l=0; l<3; l++) {
          J[0*3+l] = 0.5*(f[1][l] - f[3][l]);
          J[1*3+l] = 0.5*(f[2][l] - f[0][l]);
          J[2*3+l] = 0.5*(f[5][l] - f[4][l]);
        }
        detJ[I] = (J[0*3+0]*(J[1*3+1]*J[2*3+2] - J[1*3+2]*J[2*3+1]) +
                   J[0*3+1]*(J[1*3+2]*J[2*3+0] - J[1*3+0]*J[2*3+2]) +
                   J[0*3+2]*(J[1*3+0]*J[2*3+1] - J[1*3+1]*J[2*3+0]));
        err = PetscMemcpy(&invJ[I*9],J,9*sizeof(dReal));dCHK(err);
        err = Kernel_A_gets_inverse_A_3(&invJ[I*9],0.0);dCHK(err);
      }
    }
  }
  dFunctionReturn(0);
}
#endif

#if 0                           /* **************************** Turn this stuff off for now */
#undef __FUNCT__
#define __FUNCT__ "DohpMFSSetUp"
dErr DohpMFSSetUp(DohpMFS mfs)
{
  dMesh m=mfs->mesh;
  iMesh_Instance mi=m->mi;
  iBase_TagHandle sizetag;
  MeshListInt s=MLZ;
  size_t namelen;
  dErr err;

  dFunctionBegin;
  /* Get the requested order of elements out of the mesh. */
  err = PetscStrlen(mfs->sizetagname,&namelen);dCHK(err);
  iMesh_getTagHandle(mi,mfs->sizetagname,&sizetag,&err,(int)namelen);dICHK(mi,err);
  iMesh_getIntArrData(mi,m->r.v,m->r.s,sizetag,&s.v,&s.a,&s.s,&err);dICHK(mi,err);

  /* Assign orders to the facet space by applying the minimum rule. */
  err = DohpMFSApplyMinimumRule(mfs,&s);dCHK(err);
  /* Set up the element operations to be compatible with the quadrature and requested approximation order. */
  err = DohpMFSSetUpElementBases(mfs,&s);dCHK(err);
  /* Set up the facet-element mapping functions to be compatible with the respective approximation orders. */
  err = DohpMFSSetUpElemFacetProjections(mfs);dCHK(err);

  err = MeshListFree(s);dCHK(err);
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpMFSApplyMinimumRule"
dErr DohpMFSApplyMinimumRule(DohpMFS mfs,const MeshListInt *s)
{
  dMesh m=mfs->mesh;
  /* iMesh_Instance mi=m->mi; */
  /* dErr err; */
#define SIZE(f,m,n)             /* Decrease the order of face number f to (m,n) */
  /* Element order is 3 int values starting at s->v[i*3] */

  /* Iterate over the regions in the support of the function space and set adjacent face basis size
  * to be no more than the region basis size.  Then iterate over the faces and set adjacent edge
  * basis size to be no more than the face.  Each boundary face will be conforming so it will only
  * get set once.  Each interior face will get hit at least twice, more if it is nonconforming.
  * Each edge will get hit multiple times.  The active vertices (those adjacent to active edges) are
  * always a node. */

  dFunctionBegin;
  for (dInt i=0; i<m->r.s; i++) {      /* Each element */
    SIZE(0,0,2);
    SIZE(1,1,2);
    SIZE(2,0,2);
    SIZE(3,1,2);
    SIZE(4,1,0);
    SIZE(5,0,1);
  }
  dFunctionReturn(0);
}
#endif  /* *************************** inactive functions */


#undef __FUNCT__
#define __FUNCT__ "DohpMFSSetUpElementBases"
dErr DohpMFSSetUpElementBases(DohpMFS mfs,const MeshListInt *s)
{

  dFunctionBegin;
  dERROR(1,"not implemented");
  dValidHeader(mfs,dMFS_COOKIE,1);
  dValidIntPointer(s,2);
  dFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DohpMFSSetUpElemFacetProjections"
dErr DohpMFSSetUpElemFacetProjections(DohpMFS mfs)
{

  dFunctionBegin;
  dERROR(1,"not implemented");
  dValidHeader(mfs,dMFS_COOKIE,1);
  dFunctionReturn(0);
}

/*
  To accomodate removal of Dirichlet degrees of freedom, we make the Dirichlet face degrees of freedom point to a
  special vector.  When applying the Jacobian, we put zeros in this vector because we are solving for the homogeneous
  part of the solution.  When we are evaluating residuals, we put the actual Dirichlet values in this vector which is
  consistent with splitting the solution into a homogeneous and inhomogeneous part.  See the design documentation.
*/
#undef __FUNCT__
#define __FUNCT__ "DohpMFSSetUpBoundaryTypes"
dErr DohpMFSSetUpBoundaryTypes(DohpMFS mfs)
{

  dFunctionBegin;
  dERROR(1,"not implemented");
  dValidHeader(mfs,dMFS_COOKIE,1);
  dFunctionReturn(0);
}
