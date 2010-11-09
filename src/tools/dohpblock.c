static const char help[] = "Create a hexahedral mesh of a block domain with full connectivity.\n";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iGeom.h>
#include <iRel.h>
#include <dohpmesh.h>
#include <dohpgeom.h>
#include <dohpsys.h>

#define CHK(err) if (err) { printf("iMesh error at %s:%d\n", __FILE__, __LINE__); exit(1); }
#define ERR(str) { printf("error at %s:%d ",__FILE__,__LINE__); exit(2); }
#define NAME_LEN 256
#define DATA_LEN 256

typedef struct {
  double x0,x1,y0,y1,z0,z1;
} Box;
static void AddToFace(int **face,int *facecount,int f,int i) { face[f][facecount[f]++] = i; }
static dErr CommitToFaceSets(iMesh_Instance mesh,dMeshEH *ents,int **face,int *facecount,dMeshESH *facesets,dMeshEH *entbuf)
{
  int err;
  for (int f=0; f<6; f++) {
    for (int i=0; i<facecount[f]; i++) entbuf[i] = ents[face[f][i]];
    iMesh_addEntArrToSet(mesh,entbuf,facecount[f],facesets[f],&err);dICHK(mesh,err);
    facecount[f] = 0;           /* Clear for future use */
  }
  return 0;
}
static dErr PrintBoundingBox(iGeom_Instance geom,iBase_EntityHandle gent,const char *name)
{
  dErr err;
  double x0,y0,z0,x1,y1,z1;
  iGeom_getEntBoundBox(geom,gent,&x0,&y0,&z0,&x1,&y1,&z1,&err);dIGCHK(geom,err);
  printf("Geom `%s' bounding box (%f,%f)x(%f,%f)x(%f,%f)\n",name,x0,x1,y0,y1,z0,z1);
  return 0;
}

static dErr doMaterial(iMesh_Instance mesh)
{
  static const char matSetName[] = "MAT_SET",matNumName[] = "MAT_NUM";
  dMeshTag matSetTag,matNumTag;
  dMeshESH mat[2];
  dMeshEH *ents;
  MeshListEH r=MLZ,v=MLZ;
  MeshListInt rvo=MLZ;
  MeshListReal x=MLZ;
  dReal fx,center[3],*matnum;
  dInt nents;
  dErr err;

  dFunctionBegin;
  iMesh_createTag(mesh,matSetName,1,iBase_INTEGER,&matSetTag,&err,sizeof(matSetName));dICHK(mesh,err);
  iMesh_createTag(mesh,matNumName,1,iBase_DOUBLE,&matNumTag,&err,sizeof(matNumName));dICHK(mesh,err);
  iMesh_getEntities(mesh,0,iBase_REGION,iMesh_ALL_TOPOLOGIES,MLREF(r),&err);dICHK(mesh,err);
  iMesh_getEntArrAdj(mesh,r.v,r.s,iBase_VERTEX,MLREF(v),MLREF(rvo),&err);dICHK(mesh,err);
  iMesh_getVtxArrCoords(mesh,v.v,v.s,iBase_INTERLEAVED,MLREF(x),&err);dICHK(mesh,err);
  err = dMalloc(r.s*sizeof(ents[0]),&ents);dCHK(err);
  err = dMalloc(r.s*sizeof(matnum[0]),&matnum);dCHK(err);
  for (dInt i=0; i<2; i++) {
    iMesh_createEntSet(mesh,0,&mat[i],&err);dICHK(mesh,err);
    iMesh_setEntSetData(mesh,mat[i],matSetTag,(char*)&i,sizeof(i),&err);dICHK(mesh,err);
    nents = 0;
    for (dInt j=0; j<r.s; j++) {
      dGeomVecMeanI(8,x.v+3*rvo.v[j],center);
      fx = sqrt(dGeomDotProd(center,center)); /* material 0 if inside the unit ball, else material 1 */
      if (i == (fx < 1.0) ? 0 : 1) {
        ents[nents] = r.v[j];
        matnum[nents] = 1.0 * i;
        nents++;
      }
    }
    iMesh_addEntArrToSet(mesh,ents,nents,mat[i],&err);dICHK(mesh,err);
    iMesh_setArrData(mesh,ents,nents,matNumTag,(char*)matnum,nents*(int)sizeof(matnum[0]),&err);dICHK(mesh,err);
  }
  dFree(ents);
  MeshListFree(r); MeshListFree(v); MeshListFree(rvo); MeshListFree(x);
  dFunctionReturn(0);
}

static dErr doGlobalNumber(iMesh_Instance mesh)
{
  MeshListEH ents=MLZ;
  int owned,offset,*number;
  dMeshTag idTag;
  dErr err;

  dFunctionBegin;
  iMesh_getEntities(mesh,0,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,MLREF(ents),&err);dICHK(mesh,err);
  err = dMalloc(ents.s*sizeof(number[0]),&number);dCHK(err);
  owned = ents.s; offset = 0;
  for (int i=0; i<owned; i++) {
    number[i] = offset + i;
  }
  iMesh_createTag(mesh,"dohp_global_number",1,iBase_INTEGER,&idTag,&err,sizeof("dohp_global_number"));dICHK(mesh,err);
  iMesh_setIntArrData(mesh,ents.v,owned,idTag,number,owned,&err);dICHK(mesh,err);
  err = dFree(number);dCHK(err);
  MeshListFree(ents);
  dFunctionReturn(0);
}

static dErr doGlobalID(iMesh_Instance mesh)
{
  MeshListEH ents=MLZ;
  MeshListInt type=MLZ;
  int count[4] = {0,0,0,0};
  int owned,offset,*number;
  dMeshTag idTag;
  dErr err;

  dFunctionBegin;
  iMesh_getEntities(mesh,0,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,MLREF(ents),&err);dICHK(mesh,err);
  iMesh_getEntArrType(mesh,ents.v,ents.s,MLREF(type),&err);dICHK(mesh,err);
  err = dMalloc(ents.s*sizeof(number[0]),&number);dCHK(err);
  owned = ents.s; offset = 0;
  for (int i=0; i<owned; i++) {
    number[i] = count[type.v[i]]++;
  }
  iMesh_getTagHandle(mesh,"GLOBAL_ID",&idTag,&err,sizeof("GLOBAL_ID"));dICHK(mesh,err);
  iMesh_setIntArrData(mesh,ents.v,owned,idTag,number,owned,&err);dICHK(mesh,err);
  err = dFree(number);dCHK(err);
  MeshListFree(ents); MeshListFree(type);
  dFunctionReturn(0);
}


static dErr createUniformTags(iMesh_Instance mesh)
{
  dMeshTag itag,rtag;
  MeshListEH ents=MLZ;
  int *idata;
  double *rdata;
  dErr err;

  dFunctionBegin;
  iMesh_getEntities(mesh,0,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,MLREF(ents),&err);dICHK(mesh,err);
  err = dMalloc(ents.s*sizeof(idata[0]),&idata);dCHK(err);
  err = dMalloc(ents.s*sizeof(rdata[0]),&rdata);dCHK(err);
  for (dInt i=0; i<ents.s; i++) {
    idata[i] = -i;
    rdata[i] = -1.0*i;
  }
  iMesh_createTag(mesh,"UNIFORM_INT",1,iBase_INTEGER,&itag,&err,sizeof("UNIFORM_INT"));dICHK(mesh,err);
  iMesh_createTag(mesh,"UNIFORM_REAL",1,iBase_DOUBLE,&rtag,&err,sizeof("UNIFORM_REAL"));dICHK(mesh,err);
  iMesh_setIntArrData(mesh,ents.v,ents.s,itag,idata,ents.s,&err);dICHK(mesh,err);
  iMesh_setDblArrData(mesh,ents.v,ents.s,rtag,rdata,ents.s,&err);dICHK(mesh,err);
  MeshListFree(ents); dFree(idata); dFree(rdata);
  dFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char *argv[])
{
  const char outopts[]="",pTagName[]="OWNING_PART", pSetName[]="PARALLEL_PARTITION";
  PetscBool  do_geom=1,assoc_with_brick=0,do_color_bdy=0,do_material = 1,do_uniform = 1,do_global_number = 0,do_global_id = 1;
  PetscBool  do_partition = 1,do_pressure = 0,do_faces = 1,do_edges = 1;
  char outfile[256] = "dblock.h5m";
  char outgeom[256] = "dblock.brep";
  dInt verbose = 1;
  iMesh_Instance mesh;
  iBase_EntityHandle *entbuf;
  iBase_EntitySetHandle root,facesets[6];
  iBase_TagHandle pTag;
  MeshListEH v=MLZ,e=MLZ,f=MLZ,r=MLZ,c=MLZ;
  MeshListReal x=MLZ;
  MeshListInt s=MLZ,part=MLZ;
  dInt *face[6],facecount[6]={0};
  int err,i,j,k,m,n,p,M,N,P,I,J,K,order=iBase_INTERLEAVED;
  Box box;

  dFunctionBegin;
  err = dInitialize(&argc,&argv,(char *)0,help);dCHK(err);

  err = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"dohpblock: create cartesian meshes",NULL);dCHK(err);
  {
    char boxstr[256] = "-1:1,-1:1,-1:1",mnp[256] = "9,9,9",MNP[256] = "2,2,2";
    err = PetscOptionsInt("-verbose","verbosity of output","none",verbose,&verbose,NULL);dCHK(err);
    err = PetscOptionsBool("-do_geom","create geometric models","none",do_geom,&do_geom,NULL);dCHK(err);
    if (do_geom) {
      err = PetscOptionsBool("-assoc_with_brick","associate boundaries with brick","none",assoc_with_brick,&assoc_with_brick,NULL);dCHK(err);
    }
    err = PetscOptionsBool("-do_color_bdy","color boundary sets","none",do_color_bdy,&do_color_bdy,NULL);dCHK(err);
    err = PetscOptionsBool("-do_material","create material sets","none",do_material,&do_material,NULL);dCHK(err);
    err = PetscOptionsBool("-do_uniform","create uniform sets","none",do_uniform,&do_uniform,NULL);dCHK(err);
    err = PetscOptionsBool("-do_global_number","create global_number tags","none",do_global_number,&do_global_number,NULL);dCHK(err);
    err = PetscOptionsBool("-do_global_id","create GLOBAL_ID tags","none",do_global_id,&do_global_id,NULL);dCHK(err);
    err = PetscOptionsBool("-do_partition","create partition sets","none",do_partition,&do_partition,NULL);dCHK(err);
    err = PetscOptionsBool("-do_pressure","create pressure sets","none",do_pressure,&do_pressure,NULL);dCHK(err);
    err = PetscOptionsBool("-do_faces","create face entities","none",do_faces,&do_faces,NULL);dCHK(err);
    err = PetscOptionsBool("-do_edges","create face entities","none",do_edges,&do_edges,NULL);dCHK(err);
    err = PetscOptionsString("-box","box x0:x1,y0:y1,z0:z1","none",boxstr,boxstr,sizeof(boxstr),NULL);dCHK(err);
    err = PetscOptionsString("-mnp","number of points m,n,p","none",mnp,mnp,sizeof(mnp),NULL);dCHK(err);
    err = PetscOptionsString("-procs_mnp","number of procs M,N,P","none",MNP,MNP,sizeof(MNP),NULL);dCHK(err);
    err = PetscOptionsString("-o","outfile","none",outfile,outfile,sizeof(outfile),NULL);dCHK(err);
    if (do_geom) {
      err = PetscOptionsString("-ogeom","outfile for geometry","none",outgeom,outgeom,sizeof(outgeom),NULL);dCHK(err);
    }
    i = sscanf(boxstr,"%lf:%lf,%lf:%lf,%lf:%lf",&box.x0,&box.x1,&box.y0,&box.y1,&box.z0,&box.z1);
    if (i != 6) dERROR(PETSC_COMM_SELF,1,"Failed to parse bounding box.");
    i = sscanf(mnp,"%d,%d,%d",&m,&n,&p);
    if (i != 3) dERROR(PETSC_COMM_SELF,1,"Failed to parse size.");
    i = sscanf(MNP,"%d,%d,%d",&M,&N,&P);
    if (i != 3) dERROR(PETSC_COMM_SELF,1,"Failed to parse partition size.");
  }
  err = PetscOptionsEnd();

  iMesh_newMesh("",&mesh,&err,0);dICHK(mesh,err);
  iMesh_getRootSet(mesh,&root,&err);dICHK(mesh,err);

  /* Allocate buffers */
  err = dMallocA(m*n*p*3,&entbuf);dCHK(err); /* More than enough to hold all entities of any given type */
  for (i=0; i<6; i++) {
    int n2max = dSqrInt(dMaxInt(m,dMaxInt(n,p)));
    err = dMallocA(2*n2max,&face[i]);dCHK(err);
    iMesh_createEntSet(mesh,0,&facesets[i],&err);dICHK(mesh,err);
  }

  /* Create vertices */
  x.a = x.s = m*n*p*3; x.v = malloc(x.a*sizeof(double));
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      for (k=0; k<p; k++) {
        I = (i*n+j)*p+k;
        if (i==0) AddToFace(face,facecount,3,I);
        else if (i==m-1) AddToFace(face,facecount,1,I);
        else if (j==0) AddToFace(face,facecount,0,I);
        else if (j==n-1) AddToFace(face,facecount,2,I);
        else if (k==0) AddToFace(face,facecount,4,I);
        else if (k==p-1) AddToFace(face,facecount,5,I);
        x.v[3*I+0] = box.x0 + (box.x1-box.x0)*(1.*i/(m-1));
        x.v[3*I+1] = box.y0 + (box.y1-box.y0)*(1.*j/(n-1));
        x.v[3*I+2] = box.z0 + (box.z1-box.z0)*(1.*k/(p-1));
      }
    }
  }
  iMesh_createVtxArr(mesh,m*n*p,order,x.v,x.s,&v.v,&v.a,&v.s,&err);dICHK(mesh,err);
  err = CommitToFaceSets(mesh,v.v,face,facecount,facesets,entbuf);
  MeshListFree(x);

  /* Create regions */
  c.a = c.s = (m-1)*(n-1)*(p-1)*8; c.v = malloc(c.a*sizeof(iBase_EntityHandle)); /* connectivity */
  I=0;
  for (i=0; i<m-1; i++) {
    for (j=0; j<n-1; j++) {
      for (k=0; k<p-1; k++) {
        c.v[I++] = v.v[((i+0)*n+(j+0))*p+(k+0)];
        c.v[I++] = v.v[((i+1)*n+(j+0))*p+(k+0)];
        c.v[I++] = v.v[((i+1)*n+(j+1))*p+(k+0)];
        c.v[I++] = v.v[((i+0)*n+(j+1))*p+(k+0)];
        c.v[I++] = v.v[((i+0)*n+(j+0))*p+(k+1)];
        c.v[I++] = v.v[((i+1)*n+(j+0))*p+(k+1)];
        c.v[I++] = v.v[((i+1)*n+(j+1))*p+(k+1)];
        c.v[I++] = v.v[((i+0)*n+(j+1))*p+(k+1)];
      }
    }
  }
  if (I != c.s) dERROR(PETSC_COMM_SELF,1,"Wrong number of regions.");
  iMesh_createEntArr(mesh,iMesh_HEXAHEDRON,c.v,c.s,&r.v,&r.a,&r.s,&s.v,&s.a,&s.s,&err);dICHK(mesh,err);
  if (r.s != (m-1)*(n-1)*(p-1)) dERROR(PETSC_COMM_SELF,1,"Wrong number of regions created.");
  printf("region size %d, status size %d\n",r.s,s.s);

  if (do_global_number) {err = doGlobalNumber(mesh);dCHK(err);}
  if (do_global_id) {err = doGlobalID(mesh);dCHK(err);}

  if (do_partition) {           /* Partition tags */
    /* Create partition. */
    part.a = part.s = r.s; part.v = malloc(part.a*sizeof(int));
    for (i=0; i<m-1; i++) {
      for (j=0; j<n-1; j++) {
        for (k=0; k<p-1; k++) {
          I = i*M/(m-1); J = j*N/(n-1); K = k*P/(p-1);
          part.v[(i*(n-1)+j)*(p-1)+k] = (I*N+J)*P+K;
        }
      }
    }
    /* MATERIAL_SET is a special name associated with all iMesh instances
    * If we are using a different name, we can assume it is not special. */
    if (strcmp(pTagName,"MATERIAL_SET")) {
      iMesh_createTag(mesh,pTagName,1,iBase_INTEGER,&pTag,&err,sizeof(pTagName));dICHK(mesh,err);
    } else {
      iMesh_getTagHandle(mesh,"MATERIAL_SET",&pTag,&err,sizeof("MATERIAL_SET"));dICHK(mesh,err);
    }
    iMesh_setIntArrData(mesh,r.v,r.s,pTag,part.v,part.s,&err);dICHK(mesh,err);
    MeshListFree(part);
  }

  if (do_partition)             /* Partition sets */
  {
    int ii,jj,kk;
    iBase_EntitySetHandle partset;
    iBase_EntityHandle *entp;
    /* reuse some stuff, set up the a partition set */
    iMesh_createTag(mesh,pSetName,1,iBase_INTEGER,&pTag,&err,sizeof(pSetName));dICHK(mesh,err);
    for (i=0; i<M; i++) {
      for (j=0; j<N; j++) {
        for (k=0; k<P; k++) {
          iMesh_createEntSet(mesh,0,&partset,&err);dICHK(mesh,err);
          entp = entbuf;
          for (ii=i*(m-1)/M; ii<(i+1)*(m-1)/M; ii++) {
            for (jj=j*(n-1)/N; jj<(j+1)*(n-1)/N; jj++) {
              for (kk=k*(p-1)/P; kk<(k+1)*(p-1)/P; kk++) {
                *entp++ = r.v[(ii*(n-1)+jj)*(p-1)+kk];
              }
            }
          }
          printf("part[%d (%d,%d,%d)] has %d regions\n",(i*N+j)*P+k,i,j,k,(int)(entp-entbuf));
          iMesh_addEntArrToSet(mesh,entbuf,(int)(entp-entbuf),partset,&err);dICHK(mesh,err);
          iMesh_setEntSetIntData(mesh,partset,pTag,(i*N+j)*P+k,&err);dICHK(mesh,err);
        }
      }
    }
  }
  MeshListFree(r); MeshListFree(s); MeshListFree(c);

  if (do_faces)
  {
    /* Create faces */
    c.a = c.s = 4*((m-1)*(n-1)*p + (m-1)*n*(p-1) + m*(n-1)*(p-1)); c.v = malloc(c.a*sizeof(iBase_EntityHandle));
    I = 0;
    for (i=0; i<m-1; i++) {     /* Faces with normal pointing in positive z direction */
      for (j=0; j<n-1; j++) {
        for (k=0; k<p; k++) {
          if (k==0) AddToFace(face,facecount,4,I/4);
          if (k==p-1) AddToFace(face,facecount,5,I/4);
          c.v[I++] = v.v[((i+0)*n+(j+0))*p+k];
          c.v[I++] = v.v[((i+1)*n+(j+0))*p+k];
          c.v[I++] = v.v[((i+1)*n+(j+1))*p+k];
          c.v[I++] = v.v[((i+0)*n+(j+1))*p+k];
        }
      }
    }
    for (i=0; i<m-1; i++) {     /* Faces with normal pointing in negative y direction */
      for (j=0; j<n; j++) {
        for (k=0; k<p-1; k++) {
          if (j==0) AddToFace(face,facecount,0,I/4);
          if (j==n-1) AddToFace(face,facecount,2,I/4);
          c.v[I++] = v.v[((i+0)*n+j)*p+(k+0)];
          c.v[I++] = v.v[((i+1)*n+j)*p+(k+0)];
          c.v[I++] = v.v[((i+1)*n+j)*p+(k+1)];
          c.v[I++] = v.v[((i+0)*n+j)*p+(k+1)];
        }
      }
    }
    for (i=0; i<m; i++) {       /* Faces with normal pointing in positive x direction */
      for (j=0; j<n-1; j++) {
        for (k=0; k<p-1; k++) {
          if (i==0) AddToFace(face,facecount,3,I/4);
          if (i==m-1) AddToFace(face,facecount,1,I/4);
          c.v[I++] = v.v[(i*n+(j+0))*p+(k+0)];
          c.v[I++] = v.v[(i*n+(j+1))*p+(k+0)];
          c.v[I++] = v.v[(i*n+(j+1))*p+(k+1)];
          c.v[I++] = v.v[(i*n+(j+0))*p+(k+1)];
        }
      }
    }
    if (I != c.s) dERROR(PETSC_COMM_SELF,1, "Wrong number of faces.");
    iMesh_createEntArr(mesh,iMesh_QUADRILATERAL,c.v,c.s,&f.v,&f.a,&f.s,&s.v,&s.a,&s.s,&err);dICHK(mesh,err);
    err = CommitToFaceSets(mesh,f.v,face,facecount,facesets,entbuf);dCHK(err);
    printf("face size %d, status size %d\n",f.s,s.s);
    MeshListFree(f); MeshListFree(s); MeshListFree(c);
  }

  if (do_edges)
  {
    /* Create edges */
    c.a = c.s = 2*(m*n*(p-1) + m*(n-1)*p + (m-1)*n*p); c.v = malloc(c.a*sizeof(iBase_EntityHandle));
    I = 0;
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
        for (k=0; k<p-1; k++) {
          if (i==0) AddToFace(face,facecount,0,I/2);
          else if (i==m-1) AddToFace(face,facecount,2,I/2);
          else if (j==0) AddToFace(face,facecount,3,I/2);
          else if (j==n-1) AddToFace(face,facecount,1,I/2);
          c.v[I++] = v.v[(i*n+j)*p+(k+0)];
          c.v[I++] = v.v[(i*n+j)*p+(k+1)];
        }
      }
    }
    for (i=0; i<m; i++) {
      for (j=0; j<n-1; j++) {
        for (k=0; k<p; k++) {
          if (i==0) AddToFace(face,facecount,0,I/2);
          else if (i==m-1) AddToFace(face,facecount,2,I/2);
          else if (k==0) AddToFace(face,facecount,4,I/2);
          else if (k==p-1) AddToFace(face,facecount,5,I/2);
          c.v[I++] = v.v[(i*n+(j+0))*p+k];
          c.v[I++] = v.v[(i*n+(j+1))*p+k];
        }
      }
    }
    for (i=0; i<m-1; i++) {
      for (j=0; j<n; j++) {
        for (k=0; k<p; k++) {
          if (j==0) AddToFace(face,facecount,3,I/2);
          else if (j==n-1) AddToFace(face,facecount,1,I/2);
          else if (k==0) AddToFace(face,facecount,4,I/2);
          else if (k==p-1) AddToFace(face,facecount,5,I/2);
          c.v[I++] = v.v[((i+0)*n+j)*p+k];
          c.v[I++] = v.v[((i+1)*n+j)*p+k];
        }
      }
    }
    if (I != c.s) dERROR(PETSC_COMM_SELF,1, "Wrong number of edges.");
    iMesh_createEntArr(mesh,iMesh_LINE_SEGMENT,c.v,c.s, &e.v,&e.a,&e.s, &s.v,&s.a,&s.s,&err);dICHK(mesh,err);
    err = CommitToFaceSets(mesh,e.v,face,facecount,facesets,entbuf);dCHK(err);
    printf("edge size %d, status size %d\n",e.s,s.s);
    MeshListFree(e); MeshListFree(s); MeshListFree(c);
  }

  /* We are done with the master vertex record. */
  MeshListFree(v);

  /* Create boundary sets, these are not related to geometry here */
  {
    dMeshESH wallset,topset,bottomset,senseSet;
    iBase_TagHandle bdyTag,senseTag;
    iMesh_getTagHandle(mesh,"NEUMANN_SET",&bdyTag,&err,sizeof("NEUMANN_SET"));dICHK(mesh,err);
    iMesh_createTag(mesh,"SENSE",1,iBase_INTEGER,&senseTag,&err,sizeof "SENSE");dICHK(mesh,err);
    iMesh_createEntSet(mesh,0,&wallset,&err);dICHK(mesh,err);
    iMesh_createEntSet(mesh,0,&topset,&err);dICHK(mesh,err);
    iMesh_createEntSet(mesh,0,&bottomset,&err);dICHK(mesh,err);
    iMesh_setEntSetIntData(mesh,wallset,bdyTag,100,&err);dICHK(mesh,err);
    iMesh_setEntSetIntData(mesh,topset,bdyTag,200,&err);dICHK(mesh,err);
    iMesh_setEntSetIntData(mesh,bottomset,bdyTag,300,&err);dICHK(mesh,err);
    for (i=0; i<4; i++) {iMesh_addEntSet(mesh,facesets[i],wallset,&err);dICHK(mesh,err);}
    iMesh_addEntSet(mesh,facesets[5],topset,&err);dICHK(mesh,err);
    iMesh_addEntSet(mesh,facesets[4],bottomset,&err);dICHK(mesh,err);

    /* Deal with SENSE on the walls */
    iMesh_createEntSet(mesh,0,&senseSet,&err);dICHK(mesh,err);
    iMesh_addEntSet(mesh,facesets[2],senseSet,&err);dICHK(mesh,err);
    iMesh_addEntSet(mesh,facesets[3],senseSet,&err);dICHK(mesh,err);
    iMesh_setEntSetIntData(mesh,senseSet,senseTag,-1,&err);dICHK(mesh,err);
    iMesh_addEntSet(mesh,senseSet,wallset,&err);dICHK(mesh,err);

    /* Deal with SENSE on the bottom */
    iMesh_createEntSet(mesh,0,&senseSet,&err);dICHK(mesh,err);
    iMesh_addEntSet(mesh,facesets[4],senseSet,&err);dICHK(mesh,err);
    iMesh_setEntSetIntData(mesh,senseSet,senseTag,-1,&err);dICHK(mesh,err);
    iMesh_addEntSet(mesh,senseSet,bottomset,&err);dICHK(mesh,err);

    for (i=0; i<6; i++) {err = dFree(face[i]);}
    err = dFree(entbuf);dCHK(err);
  }

  if (do_material) {err = doMaterial(mesh);dCHK(err);}

  /* Add a real valued tag over the vertices. */
  if (do_pressure) {
    static const char *myTagName = "pressure";
    iBase_TagHandle myTag;
    double *myData;

    iMesh_getEntities(mesh,0,iBase_VERTEX,iMesh_POINT,&v.v,&v.a,&v.s,&err);dICHK(mesh,err);
    iMesh_createTag(mesh,myTagName,1,iBase_DOUBLE,&myTag,&err,(int)strlen(myTagName));dICHK(mesh,err);
    err = PetscMalloc(v.s*sizeof(double),&myData);dCHK(err);
    for (i=0; i<v.s; i++) { myData[i] = 1.0 * i; }
    iMesh_setDblArrData(mesh,v.v,v.s,myTag,myData,v.s,&err);dICHK(mesh,err);
    err = PetscFree(myData);dCHK(err);
    MeshListFree(v);
  }

  if (do_uniform) {err = createUniformTags(mesh);dCHK(err);}

  if (do_geom)
  {
    char geom_options[] = ";ENGINE=OCC;",
      geom_save_options[] = ";TYPE=OCC;";
    iGeom_Instance geom;
    iRel_Instance assoc;
    iRel_RelationHandle rel;
    iBase_EntityHandle brick;
    iGeom_newGeom(geom_options,&geom,&err,sizeof geom_options);dIGCHK(geom,err);
    iRel_newAssoc(0,&assoc,&err,0);dIRCHK(assoc,err);
    iRel_createAssociation(assoc,geom,0,iRel_IGEOM_IFACE,mesh,1,iRel_IMESH_IFACE,&rel,&err);dIGCHK(assoc,err);
    iGeom_createBrick(geom,box.x1-box.x0,box.y1-box.y0,box.z1-box.z0,&brick,&err);dIGCHK(geom,err);
    iGeom_moveEnt(geom,brick,0.5*(box.x0+box.x1),0.5*(box.y0+box.y1),0.5*(box.z0+box.z1),&err);dIGCHK(geom,err);
    err = PrintBoundingBox(geom,brick,"brick");dCHK(err);
    {
      iBase_EntityHandle gface[6],*gface_p=gface;
      int gface_a=6,gface_s;
      iGeom_getEntAdj(geom,brick,2,&gface_p,&gface_a,&gface_s,&err);dIGCHK(geom,err);
      for (i=0; i<6; i++) {
        char name[20];
        sprintf(name,"face_%d",i);
        err = PrintBoundingBox(geom,gface[i],name);dCHK(err);
      }
      if (assoc_with_brick) {
        for (i=0; i<6; i++) {
          iRel_setEntSetAssociation(assoc,rel,brick,facesets[i],&err);dIRCHK(assoc,err);
        }
      } else {
        /* Set associations.  With the current Lasso implementation, these will not be saved */
        iRel_setEntSetAssociation(assoc,rel,gface[0],facesets[3],&err);dIRCHK(assoc,err);
        iRel_setEntSetAssociation(assoc,rel,gface[1],facesets[1],&err);dIRCHK(assoc,err);
        iRel_setEntSetAssociation(assoc,rel,gface[2],facesets[0],&err);dIRCHK(assoc,err);
        iRel_setEntSetAssociation(assoc,rel,gface[3],facesets[2],&err);dIRCHK(assoc,err);
        iRel_setEntSetAssociation(assoc,rel,gface[4],facesets[4],&err);dIRCHK(assoc,err);
        iRel_setEntSetAssociation(assoc,rel,gface[5],facesets[5],&err);dIRCHK(assoc,err);
      }
    }
    {
      dMeshTag meshGlobalIDTag,meshGeomDimTag,geomGlobalIDTag;
      /* Manually set association tags, these are set so that the associations above can be inferred. */
      iMesh_getTagHandle(mesh,"GLOBAL_ID",&meshGlobalIDTag,&err,sizeof "GLOBAL_ID");dICHK(mesh,err);
      iMesh_getTagHandle(mesh,"GEOM_DIMENSION",&meshGeomDimTag,&err,sizeof "GEOM_DIMENSION");dICHK(mesh,err);
      iGeom_getTagHandle(geom,"GLOBAL_ID",&geomGlobalIDTag,&err,sizeof "GLOBAL_ID");dIGCHK(geom,err);
      for (i=0; i<6; i++) {
        iBase_EntityHandle gface;
        int gid,gdim;
        iRel_getSetEntAssociation(assoc,rel,facesets[i],1,&gface,&err);dIRCHK(assoc,err);
        iGeom_getEntType(geom,gface,&gdim,&err);dIGCHK(geom,err);
        if (gdim != 2) dERROR(PETSC_COMM_SELF,1,"Geometric dimension is %d, expected 2",gdim);
        iGeom_getIntData(geom,gface,geomGlobalIDTag,&gid,&err);dIGCHK(geom,err);
        iMesh_setEntSetIntData(mesh,facesets[i],meshGeomDimTag,2,&err);dICHK(mesh,err);
        /* If the following line is disabled, Lasso will pick up the wrong relations, but at least they will still be with
        * surfaces.  Wouldn't it be better to not find relations? */
        iMesh_setEntSetIntData(mesh,facesets[i],meshGlobalIDTag,gid,&err);dICHK(mesh,err);
      }
    }
    iGeom_save(geom,outgeom,geom_save_options,&err,sizeof outgeom,sizeof geom_save_options);dIGCHK(geom,err);
    iGeom_dtor(geom,&err);dIGCHK(geom,err);
    iRel_dtor(assoc,&err);dIRCHK(assoc,err);
  }

  iMesh_save(mesh,0,outfile,outopts,&err,(int)sizeof(outfile),(int)strlen(outopts));dICHK(mesh,err);
  err = dFinalize();dCHK(err);
  dFunctionReturn(0);
}
