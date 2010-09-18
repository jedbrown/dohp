#include <dohpfsimpl.h>
#include <dohpvec.h>

dErr dFSGetCoordinateFS(dFS fs,dFS *incfs)
{
  dErr err;

  dFunctionBegin;
  if (!fs->geometry.fs) {
    dFS cfs;
    dMesh mesh;
    dJacobi jacobi;
    dMeshTag dtag;
    char degreename[256],prefix[256];

    err = dFSGetMesh(fs,&mesh);dCHK(err);
    err = dFSGetJacobi(fs,&jacobi);dCHK(err);
    err = dFSCreate(((dObject)fs)->comm,&cfs);dCHK(err);
    err = dFSSetMesh(cfs,mesh,fs->set.active);dCHK(err);
    err = PetscSNPrintf(degreename,sizeof degreename,"%scfs_degree",((dObject)fs)->prefix);dCHK(err);
    err = dMeshCreateRuleTagIsotropic(mesh,fs->set.active,degreename,dPolynomialOrderCreate(0,1,1,1),&dtag);dCHK(err);
    err = dFSSetDegree(cfs,jacobi,dtag);dCHK(err);
    err = dFSSetRuleTag(cfs,jacobi,dtag);dCHK(err); /* FIXME: remove this attribute from dFS */
    err = dFSSetBlockSize(cfs,3);dCHK(err);
    err = PetscSNPrintf(prefix,sizeof prefix,"%scoord_",((dObject)fs)->prefix);dCHK(err);
    err = dFSSetOptionsPrefix(cfs,prefix);dCHK(err);
    err = dFSSetFromOptions(cfs);dCHK(err);
    fs->geometry.fs = cfs;
  }
  *incfs = fs->geometry.fs;
  dFunctionReturn(0);
}
