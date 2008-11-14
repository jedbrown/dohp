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
  const dInt rsize[] = {5,5,5},dsize[] = {4,4,4},dall[] = {3,3,3};
  iMesh_Instance mi;
  dJacobi jac;
  dFS fs;
  dMesh mesh;
  dMeshESH domain;
  dMeshTag rtag,dtag,ownertag,dalltag;
  MPI_Comm comm;
  PetscViewer viewer;
  PetscMPIInt rank;
  dErr err;

  err = PetscInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  err = MPI_Comm_rank(comm,&rank);dCHK(err);
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = dMeshCreate(comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);

  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  iMesh_getRootSet(mi,&domain,&err);dICHK(mi,err);

  err = createIsotropicIntTag(mi,domain,iBase_REGION,iMesh_HEXAHEDRON,3,rsize,"region_rule",&rtag);dCHK(err);
  err = createIsotropicIntTag(mi,domain,iBase_REGION,iMesh_HEXAHEDRON,3,dsize,"region_degree",&dtag);dCHK(err);
  err = createIsotropicIntTag(mi,domain,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,3,dall,"all_degree",&dalltag);dCHK(err);
  err = createIsotropicIntTag(mi,domain,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,1,&rank,"owner",&ownertag);dCHK(err);

  {
    dInt i,n,*deg;
    dMeshEH *ents;
    dEntTopology *topo;
    CHKMEMQ;
    err = dMeshGetNumEnts(mesh,domain,dTYPE_ALL,dTOPO_ALL,&n);
    err = dMallocA3(n,&ents,n,&topo,3*n,&deg);dCHK(err);
    {
      dInt used;
      err = dMeshGetEnts(mesh,domain,dTYPE_ALL,dTOPO_ALL,ents,n,&used);dCHK(err);
      if (used != n) dERROR(1,"wrong number of entities");
    }
    err = dMeshGetTopo(mesh,n,ents,topo);dCHK(err);
    for (i=0; i<n; i++) {
      switch (topo[i]) {
        case dTOPO_POINT: deg[3*i] = deg[3*i+1] = deg[3*i+2] = 1; break;
        case dTOPO_LINE:  deg[3*i] = 5; deg[3*i+1] = deg[3*i+2] = 1; break;
        case dTOPO_QUAD:  deg[3*i] = deg[3*i+1] = 6; deg[3*i+2] = 1; break;
        case dTOPO_HEX:   deg[3*i] = deg[3*i+1] = deg[3*i+2] = 7; break;
        default: dERROR(1,"Topology not supported");
      }
      CHKMEMQ;
    }
    err = dMeshTagSetData(mesh,dalltag,ents,n,deg,3*n,dDATA_INT);dCHK(err);
    err = dFree3(ents,topo,deg);dCHK(err);
  }

  err = dMeshTagBcast(mesh,ownertag);dCHK(err);

  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);

  err = dFSCreate(comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,domain);dCHK(err);
  err = dFSSetRuleTag(fs,jac,rtag);dCHK(err);
  err = dFSSetDegree(fs,jac,dalltag);dCHK(err);
  err = dFSSetFromOptions(fs);dCHK(err);

  {
    char outname[256];
    int size;
    MPI_Comm_size(comm,&size);
    snprintf(outname,sizeof(outname),"foo%d-%d.vtk",size,rank);
    iMesh_save(mi,0,outname,"",&err,sizeof(outname),0);dCHK(err);
  }

  err = dFSDestroy(fs);dCHK(err);
  err = dJacobiDestroy(jac);dCHK(err);
  err = dMeshDestroy(mesh);
  err = PetscFinalize();dCHK(err);
  return 0;
}
