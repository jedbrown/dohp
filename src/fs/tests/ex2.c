static const char help[] = "Test the construction of dFS objects.\n";

#include <dohpfs.h>

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
  const dInt rsize[] = {5,5,5},dsize[] = {4,4,4};
  iMesh_Instance mi;
  dJacobi jac;
  dFS fs;
  dMesh mesh;
  dMeshESH domain;
  dMeshTag rtag,dtag,ownertag;
  MPI_Comm comm;
  PetscViewer viewer;
  PetscMPIInt rank;
  dErr err;

  err = dInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  err = MPI_Comm_rank(comm,&rank);dCHK(err);
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = dMeshCreate(comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"zdblock.h5m",NULL);dCHK(err);
  err = dMeshSetType(mesh,dMESHSERIAL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);

  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  iMesh_getRootSet(mi,&domain,&err);dICHK(mi,err);

  err = createIsotropicIntTag(mi,domain,iBase_REGION,iMesh_HEXAHEDRON,3,rsize,"region_rule",&rtag);
  err = createIsotropicIntTag(mi,domain,iBase_REGION,iMesh_HEXAHEDRON,3,dsize,"region_degree",&dtag);

  err = createIsotropicIntTag(mi,domain,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,1,&rank,"owner",&ownertag);

  err = dMeshTagBcast(mesh,ownertag);dCHK(err);
  //err = dMeshView(mesh,viewer);dCHK(err);

  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);

  err = dFSCreate(comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,domain);dCHK(err);
  err = dFSSetRuleTag(fs,jac,rtag);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  err = dFSSetFromOptions(fs);dCHK(err);

  {
    char outname[256];
    int size;
    MPI_Comm_size(comm,&size);
    snprintf(outname,sizeof(outname),"foo%d-%d.vtk",size,rank);
    iMesh_save(mi,0,outname,"",&err,sizeof(outname),0);dCHK(err);
  }

  err = dFSDestroy(&fs);dCHK(err);
  err = dJacobiDestroy(&jac);dCHK(err);
  err = dMeshDestroy(&mesh);
  err = dFinalize();dCHK(err);
  return 0;
}
