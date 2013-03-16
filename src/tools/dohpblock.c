static const char help[] = "Create a hexahedral mesh of a block domain with full connectivity.\n";

#include <dohpmesh.h>
#include <dohp.h>
#include <dohpsys.h>

int main(int argc, char *argv[])
{
  dErr err;
  dMesh mesh;
  PetscBool do_geom;
  char outfile[256] = "dblock.h5m";
  char outopts[]="";
  char outgeom[256] = "dblock.brep";
  iBase_EntitySetHandle root;
  iMesh_Instance mi;

  dFunctionBegin;
  err = dInitialize(&argc,&argv,(char *)0,help);dCHK(err);

  err = dMeshCreate(PETSC_COMM_WORLD,&mesh);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshGetRoot(mesh,&root);dCHK(err);

  err = PetscOptionsBegin(PetscObjectComm((PetscObject)mesh),NULL,"dohpblock: create cartesian meshes",NULL);dCHK(err);
  {
    err = PetscOptionsString("-o","outfile","none",outfile,outfile,sizeof(outfile),NULL);dCHK(err);
    err = PetscOptionsBool("-do_geom","create geometric models","none",do_geom=dTRUE,&do_geom,NULL);dCHK(err);
    if (do_geom) {
      err = PetscOptionsString("-ogeom","outfile for geometry","none",outgeom,outgeom,sizeof(outgeom),NULL);dCHK(err);
    }
  }
  err = PetscOptionsEnd();

  err = dMeshGenerateBlock(mesh,root,do_geom);dCHK(err);
#ifdef dHAVE_ITAPS_REL
  if (do_geom) {
    const char geom_save_options[] = ";TYPE=OCC;";
    iGeom_Instance geom;
    err = dMeshGetGeometryRelation(mesh,&geom,NULL);dCHK(err);
    iGeom_save(geom,outgeom,geom_save_options,&err,sizeof outgeom,sizeof geom_save_options);dIGCHK(geom,err);
  }
#endif

  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  iMesh_save(mi,root,outfile,outopts,&err,(int)sizeof(outfile),(int)strlen(outopts));dICHK(mi,err);
  err = dMeshDestroy(&mesh);dCHK(err);
  err = dFinalize();dCHK(err);
  dFunctionReturn(0);
}
