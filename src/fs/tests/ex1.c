static const char help[] = "Test the construction of dFS objects.\n";

#include "dohpfs.h"
#include "dohp.h"

static dErr createIsotropicIntTag(iMesh_Instance mi,dMeshESH set,dEntType type,dEntTopology topo,int n,const int v[],const char name[],dMeshTag *intag)
{
  MeshListEH ents=MLZ;
  dMeshTag tag;
  dInt *values;
  dErr err;

  dFunctionBegin;
  dValidPointer(intag,8);
  *intag = 0;
  iMesh_createTag(mi,name,n,iBase_INTEGER,&tag,&err,(int)strlen(name));dCHK(err);
  iMesh_getEntities(mi,set,type,topo,&ents.v,&ents.a,&ents.s,&err);dICHK(mi,err);
  err = dMalloc(n*ents.s*sizeof(*values),&values);dCHK(err);
  for (dInt i=0; i<ents.s; i++) {
    err = dMemcpy(&values[n*i],v,n*sizeof(v[0]));dCHK(err);
  }
  iMesh_setIntArrData(mi,ents.v,ents.s,tag,values,n*ents.s,&err);dICHK(mi,err);
  err = dFree(values);dCHK(err);
  *intag = tag;
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  const char pTagName[] = "dohp_partition";
  const dInt rsize[] = {5,5,5},dsize[] = {4,4,4};
  iMesh_Instance mi;
  dJacobi jac;
  dFS fs;
  dMesh mesh;
  dMeshESH domain;
  dMeshTag partition,rtag,dtag;
  MPI_Comm comm;
  PetscViewer viewer;
  dErr err;

  dFunctionBegin;
  err = PetscInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = dMeshCreate(comm,&mesh);dCHK(err);
  err = dMeshLoad(mesh,"dblock.h5m","");dCHK(err);
  err = dMeshView(mesh,viewer);dCHK(err);
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  iMesh_getRootSet(mi,&domain,&err);dICHK(mi,err);
  iMesh_getTagHandle(mi,pTagName,&partition,&err,strlen(pTagName));dICHK(mi,err);
  err = createIsotropicIntTag(mi,domain,iBase_REGION,iMesh_HEXAHEDRON,3,rsize,"region_rule",&rtag);
  err = createIsotropicIntTag(mi,domain,iBase_REGION,iMesh_HEXAHEDRON,3,dsize,"region_degree",&dtag);

  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);

  err = dFSCreate(comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,domain,partition);dCHK(err);
  err = dFSSetRuleTag(fs,jac,rtag);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  err = dFSSetFromOptions(fs);dCHK(err);

  err = dFSDestroy(fs);dCHK(err);
  err = dJacobiDestroy(jac);dCHK(err);
  err = dMeshDestroy(mesh);
  err = PetscFinalize();dCHK(err);
  dFunctionReturn(0);
}
