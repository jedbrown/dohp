static const char help[] = "Create a hexahedral mesh of a block domain with full connectivity.\n";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "private/dohpimpl.h"
#include "dohpgeom.h"

#define CHK(err) if (err) { printf("iMesh error at %s:%d\n", __FILE__, __LINE__); exit(1); }
#define ERR(str) { printf("error at %s:%d ",__FILE__,__LINE__); exit(2); }
#define NAME_LEN 256
#define DATA_LEN 256

typedef struct {
  double x0,x1,y0,y1,z0,z1;
} Box;
static int on_wall(Box *b,double xx[])
{
  return (xx[0] == b->x0 || xx[0] == b->x1 || xx[1] == b->y0 || xx[1] == b->y1 || xx[2] == b->z0);
}
static int on_lid(Box *b,double xx[])  {
  return (xx[2] == b->z1);
}
typedef int OnBdyFunc(Box *box,double[]);


static dErr doMaterial(iMesh_Instance mesh)
{
  static const char matSetName[] = "MAT_SET",matNumName[] = "MAT_NUM";
  int order = iBase_INTERLEAVED;
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
  iMesh_createTag(mesh,matSetName,1,iBase_INTEGER,&matSetTag,&err,strlen(matSetName));dICHK(mesh,err);
  iMesh_createTag(mesh,matNumName,1,iBase_DOUBLE,&matNumTag,&err,strlen(matNumName));dICHK(mesh,err);
  iMesh_getEntities(mesh,0,iBase_REGION,iMesh_ALL_TOPOLOGIES,MLREF(r),&err);dICHK(mesh,err);
  iMesh_getEntArrAdj(mesh,r.v,r.s,iBase_VERTEX,MLREF(v),MLREF(rvo),&err);dICHK(mesh,err);
  iMesh_getVtxArrCoords(mesh,v.v,v.s,&order,MLREF(x),&err);dICHK(mesh,err);
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
    iMesh_addEntArrToSet(mesh,ents,nents,&mat[i],&err);dICHK(mesh,err);
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
  iMesh_getTagHandle(mesh,"GLOBAL_ID",&idTag,&err,strlen("GLOBAL_ID"));dICHK(mesh,err);
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
  PetscTruth do_bdy = 0,do_material = 1,do_uniform = 1,do_global_number = 0,do_global_id = 1;
  PetscTruth do_partition = 1,do_pressure = 0,do_faces = 1,do_edges = 1,do_orient = 0;
  char outfile[256] = "dblock.h5m";
  dInt verbose = 1;
  iMesh_Instance mesh;
  iBase_EntitySetHandle root;
  iBase_TagHandle pTag,feOrientTag,feAdjTag,rfOrientTag,rfAdjTag;
  MeshListEH v=MLZ,e=MLZ,f=MLZ,r=MLZ,c=MLZ,ev=MLZ,fv=MLZ,rv=MLZ;
  MeshListReal x=MLZ;
  MeshListInt s=MLZ,part=MLZ,feo=MLZ,in=MLZ,fvo=MLZ,evo=MLZ,rfo=MLZ,rvo=MLZ;
  dInt *feOrient,*rfOrient;
  dInt feOrientSize,rfOrientSize;
  int err,i,j,k,m,n,p,M,N,P,I,J,K,order=iBase_INTERLEAVED;
  Box box;

  dFunctionBegin;
  err = PetscInitialize(&argc,&argv,(char *)0,help);dCHK(err);

  err = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"dohpblock: create cartesian meshes",NULL);dCHK(err);
  {
    char boxstr[256] = "-1:1,-1:1,-1:1",mnp[256] = "9,9,9",MNP[256] = "2,2,2";
    err = PetscOptionsInt("-verbose","verbosity of output","none",verbose,&verbose,NULL);
    err = PetscOptionsTruth("-do_bdy","create boundary sets","none",do_bdy,&do_bdy,NULL);dCHK(err);
    err = PetscOptionsTruth("-do_material","create material sets","none",do_material,&do_material,NULL);dCHK(err);
    err = PetscOptionsTruth("-do_uniform","create uniform sets","none",do_uniform,&do_uniform,NULL);dCHK(err);
    err = PetscOptionsTruth("-do_global_number","create global_number tags","none",do_global_number,&do_global_number,NULL);dCHK(err);
    err = PetscOptionsTruth("-do_global_id","create GLOBAL_ID tags","none",do_global_id,&do_global_id,NULL);dCHK(err);
    err = PetscOptionsTruth("-do_partition","create partition sets","none",do_partition,&do_partition,NULL);dCHK(err);
    err = PetscOptionsTruth("-do_pressure","create pressure sets","none",do_pressure,&do_pressure,NULL);dCHK(err);
    err = PetscOptionsTruth("-do_faces","create face entities","none",do_faces,&do_faces,NULL);dCHK(err);
    err = PetscOptionsTruth("-do_edges","create face entities","none",do_edges,&do_edges,NULL);dCHK(err);
    if (do_faces && do_edges) {
      err = PetscOptionsTruth("-do_orient","create explicit orientation for faces and edges","none",do_orient,&do_orient,NULL);dCHK(err);
    }
    err = PetscOptionsString("-box","box x0:x1,y0:y1,z0:z1","none",boxstr,boxstr,sizeof(boxstr),NULL);dCHK(err);
    err = PetscOptionsString("-mnp","number of points m,n,p","none",mnp,mnp,sizeof(mnp),NULL);dCHK(err);
    err = PetscOptionsString("-procs_mnp","number of procs M,N,P","none",MNP,MNP,sizeof(MNP),NULL);dCHK(err);
    err = PetscOptionsString("-o","outfile","none",outfile,outfile,sizeof(outfile),NULL);dCHK(err);
    i = sscanf(boxstr,"%lf:%lf,%lf:%lf,%lf:%lf",&box.x0,&box.x1,&box.y0,&box.y1,&box.z0,&box.z1);
    if (i != 6) dERROR(1,"Failed to parse bounding box.");
    i = sscanf(mnp,"%d,%d,%d",&m,&n,&p);
    if (i != 3) dERROR(1,"Failed to parse size.");
    i = sscanf(MNP,"%d,%d,%d",&M,&N,&P);
    if (i != 3) dERROR(1,"Failed to parse partition size.");
  }
  err = PetscOptionsEnd();

  iMesh_newMesh("", &mesh, &err, 0);dICHK(mesh,err);
  iMesh_getRootSet(mesh, &root, &err);dICHK(mesh,err);

  /* Create vertices */
  x.a = x.s = m*n*p*3; x.v = malloc(x.a*sizeof(double));
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      for (k=0; k<p; k++) {
        I = (i*n+j)*p+k;
        x.v[3*I+0] = box.x0 + (box.x1-box.x0)*(1.*i/(m-1));
        x.v[3*I+1] = box.y0 + (box.y1-box.y0)*(1.*j/(n-1));
        x.v[3*I+2] = box.z0 + (box.z1-box.z0)*(1.*k/(p-1));
      }
    }
  }
  iMesh_createVtxArr(mesh,m*n*p,order,x.v,x.s,&v.v,&v.a,&v.s,&err);dICHK(mesh,err);
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
  if (I != c.s) dERROR(1,"Wrong number of regions.");
  iMesh_createEntArr(mesh,iMesh_HEXAHEDRON,c.v,c.s,&r.v,&r.a,&r.s,&s.v,&s.a,&s.s,&err);dICHK(mesh,err);
  if (r.s != (m-1)*(n-1)*(p-1)) dERROR(1,"Wrong number of regions created.");
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
      iMesh_getTagHandle(mesh,"MATERIAL_SET",&pTag,&err,strlen("MATERIAL_SET"));dICHK(mesh,err);
    }
    iMesh_setIntArrData(mesh,r.v,r.s,pTag,part.v,part.s,&err);dICHK(mesh,err);
    MeshListFree(part);
  }

  if (do_partition)             /* Partition sets */
  {
    int ii,jj,kk;
    iBase_EntitySetHandle partset;
    iBase_EntityHandle *entbuf,*entp;
    /* reuse some stuff, set up the a partition set */
    err = dMallocA(r.s,&entbuf);dCHK(err);
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
          iMesh_addEntArrToSet(mesh,entbuf,(int)(entp-entbuf),&partset,&err);dICHK(mesh,err);
          iMesh_setEntSetIntData(mesh,partset,pTag,(i*N+j)*P+k,&err);dICHK(mesh,err);
        }
      }
    }
    err = dFree(entbuf);dCHK(err);
  }
  MeshListFree(r); MeshListFree(s); MeshListFree(c);

  if (do_faces)
  {
    /* Create faces */
    c.a = c.s = 4*((m-1)*(n-1)*p + (m-1)*n*(p-1) + m*(n-1)*(p-1)); c.v = malloc(c.a*sizeof(iBase_EntityHandle));
    I = 0;
    for (i=0; i<m-1; i++) {
      for (j=0; j<n-1; j++) {
        for (k=0; k<p; k++) {
          c.v[I++] = v.v[((i+0)*n+(j+0))*p+k];
          c.v[I++] = v.v[((i+1)*n+(j+0))*p+k];
          c.v[I++] = v.v[((i+1)*n+(j+1))*p+k];
          c.v[I++] = v.v[((i+0)*n+(j+1))*p+k];
        }
      }
    }
    for (i=0; i<m-1; i++) {
      for (j=0; j<n; j++) {
        for (k=0; k<p-1; k++) {
          c.v[I++] = v.v[((i+0)*n+j)*p+(k+0)];
          c.v[I++] = v.v[((i+1)*n+j)*p+(k+0)];
          c.v[I++] = v.v[((i+1)*n+j)*p+(k+1)];
          c.v[I++] = v.v[((i+0)*n+j)*p+(k+1)];
        }
      }
    }
    for (i=0; i<m; i++) {
      for (j=0; j<n-1; j++) {
        for (k=0; k<p-1; k++) {
          c.v[I++] = v.v[(i*n+(j+0))*p+(k+0)];
          c.v[I++] = v.v[(i*n+(j+1))*p+(k+0)];
          c.v[I++] = v.v[(i*n+(j+1))*p+(k+1)];
          c.v[I++] = v.v[(i*n+(j+0))*p+(k+1)];
        }
      }
    }
    if (I != c.s) dERROR(1, "Wrong number of faces.");
    iMesh_createEntArr(mesh,iMesh_QUADRILATERAL,c.v,c.s,&f.v,&f.a,&f.s,&s.v,&s.a,&s.s,&err);dICHK(mesh,err);
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
          c.v[I++] = v.v[(i*n+j)*p+(k+0)];
          c.v[I++] = v.v[(i*n+j)*p+(k+1)];
        }
      }
    }
    for (i=0; i<m; i++) {
      for (j=0; j<n-1; j++) {
        for (k=0; k<p; k++) {
          c.v[I++] = v.v[(i*n+(j+0))*p+k];
          c.v[I++] = v.v[(i*n+(j+1))*p+k];
        }
      }
    }
    for (i=0; i<m-1; i++) {
      for (j=0; j<n; j++) {
        for (k=0; k<p; k++) {
          c.v[I++] = v.v[((i+0)*n+j)*p+k];
          c.v[I++] = v.v[((i+1)*n+j)*p+k];
        }
      }
    }
    if (I != c.s) dERROR(1, "Wrong number of edges.");
    iMesh_createEntArr(mesh,iMesh_LINE_SEGMENT,c.v,c.s, &e.v,&e.a,&e.s, &s.v,&s.a,&s.s,&err);dICHK(mesh,err);
    printf("edge size %d, status size %d\n",e.s,s.s);
    MeshListFree(e); MeshListFree(s); MeshListFree(c);
  }

  /* We are done with the master vertex record and ready to set up orientations. */
  MeshListFree(v);

  if (do_orient)
  {
    /* Orient Edges with respect to Faces */
    iMesh_getEntities(mesh,0,iBase_FACE,iMesh_QUADRILATERAL,&f.v,&f.a,&f.s,&err);dICHK(mesh,err);
    iMesh_getAdjEntities(mesh,0,iBase_FACE,iMesh_QUADRILATERAL,iBase_EDGE,
                         &e.v,&e.a,&e.s, &feo.v,&feo.a,&feo.s, &in.v,&in.a,&in.s, &err);dICHK(mesh,err);
    MeshListFree(in);
    iMesh_getEntArrAdj(mesh,f.v,f.s,iBase_VERTEX,&fv.v,&fv.a,&fv.s,&fvo.v,&fvo.a,&fvo.s,&err);dICHK(mesh,err);
    iMesh_getEntArrAdj(mesh,e.v,e.s,iBase_VERTEX,&ev.v,&ev.a,&ev.s,&evo.v,&evo.a,&evo.s,&err);dICHK(mesh,err);
    printf("f %d, e %d, fv %d, feo %d, fvo %d, ev %d, evo %d\n", f.s,e.s,fv.s,feo.s,fvo.s,ev.s,evo.s);
    feOrientSize = e.s;
    err = PetscMalloc(feOrientSize*sizeof(feOrient[0]),&feOrient);dCHK(err);
    for (i=0; i<f.s; i++) {      /* Loop over faces */
      if (verbose > 2) {
        printf("face[%d] vertex: %ld %ld %ld %ld\n",i,
               (long)fv.v[fvo.v[i]],  (long)fv.v[fvo.v[i]+1],
               (long)fv.v[fvo.v[i]+2],(long)fv.v[fvo.v[i]+3]);
      }
      for (j=0; j<4; j++) {       /* Loop over edges adjacent to this face */
        if (verbose > 2) {
          printf("edge[%d][%d] %ld %ld\n",i,j,(long)ev.v[evo.v[feo.v[i]+j]],(long)ev.v[evo.v[feo.v[i]+j]+1]);
        }
        err = dGeomOrientFindPerm_QuadLine(&fv.v[fvo.v[i]],&ev.v[evo.v[feo.v[i]+j]],j,&feOrient[feo.v[i]+j]);dCHK(err);
      }
    }
    iMesh_createTag(mesh,dTAG_ORIENT_FACE_EDGE,4*sizeof(feOrient[0]),iBase_BYTES,
                    &feOrientTag,&err,strlen(dTAG_ORIENT_FACE_EDGE));dICHK(mesh,err);
    iMesh_setArrData(mesh,f.v,f.s,feOrientTag,(const char*)feOrient,feOrientSize*(dInt)sizeof(feOrient[0]),&err);dICHK(mesh,err);
    iMesh_createTag(mesh,dTAG_ADJ_FACE_EDGE,4,iBase_ENTITY_HANDLE,&feAdjTag,&err,strlen(dTAG_ADJ_FACE_EDGE));dICHK(mesh,err);
    iMesh_setEHArrData(mesh,f.v,f.s,feAdjTag,e.v,e.s,&err);dICHK(mesh,err);

    MeshListFree(f); MeshListFree(e); MeshListFree(feo); MeshListFree(fv); MeshListFree(fvo);
    MeshListFree(ev); MeshListFree(evo);
    err = PetscFree(feOrient);dCHK(err);

    /* Orient Faces with respect to Regions. */
    iMesh_getEntities(mesh,0,iBase_REGION,iMesh_HEXAHEDRON,&r.v,&r.a,&r.s,&err);dICHK(mesh,err);
    iMesh_getAdjEntities(mesh,0,iBase_REGION,iMesh_HEXAHEDRON,iBase_FACE,
                         &f.v,&f.a,&f.s, &rfo.v,&rfo.a,&rfo.s, &in.v,&in.a,&in.s, &err);dICHK(mesh,err);
    MeshListFree(in);
    iMesh_getEntArrAdj(mesh,r.v,r.s,iBase_VERTEX,&rv.v,&rv.a,&rv.s,&rvo.v,&rvo.a,&rvo.s,&err);dICHK(mesh,err);
    iMesh_getEntArrAdj(mesh,f.v,f.s,iBase_VERTEX,&fv.v,&fv.a,&fv.s,&fvo.v,&fvo.a,&fvo.s,&err);dICHK(mesh,err);
    printf("r %d, f %d, rv %d, rfo %d, rvo %d, fv %d, fvo %d\n", r.s,f.s,rv.s,rfo.s,rvo.s,fv.s,fvo.s);
    rfOrientSize = f.s;
    err = PetscMalloc(rfOrientSize*sizeof(rfOrient[0]),&rfOrient);dCHK(err);
    for (i=0; i<r.s && i<1e8; i++) {
      if (verbose > 2) {
        printf("region[%d]",i);     /* Vertices of this region */
        for (j=0; j<8; j++) { if (j%4==0) printf("  "); printf(" %3ld",(long)rv.v[rvo.v[i]+j]); }
        printf("\n");
      }
      for (j=0; j<6; j++) {       /* Faces of this region */
        if (verbose > 2) {
          printf("face[%d][%d] %3ld %3ld %3ld %3ld\n",i,j,(long)fv.v[fvo.v[rfo.v[i]+j]],
                 (long)fv.v[fvo.v[rfo.v[i]+j]+1],(long)fv.v[fvo.v[rfo.v[i]+j]+2],(long)fv.v[fvo.v[rfo.v[i]+j]+3]);
        }
        err = dGeomOrientFindPerm_HexQuad(&rv.v[rvo.v[i]],&fv.v[fvo.v[rfo.v[i]+j]],j,&rfOrient[rfo.v[i]+j]);dCHK(err);
      }
    }
    iMesh_createTag(mesh,dTAG_ORIENT_REGION_FACE,6*sizeof(rfOrient[0]),iBase_BYTES,
                    &rfOrientTag,&err,strlen(dTAG_ORIENT_REGION_FACE));dICHK(mesh,err);
    iMesh_setArrData(mesh,r.v,r.s,rfOrientTag,(const char*)rfOrient,rfOrientSize*(dInt)sizeof(rfOrient[0]),&err);dICHK(mesh,err);
    iMesh_createTag(mesh,dTAG_ADJ_REGION_FACE,6,iBase_ENTITY_HANDLE,&rfAdjTag,&err,strlen(dTAG_ADJ_REGION_FACE));dICHK(mesh,err);
    iMesh_setEHArrData(mesh,r.v,r.s,rfAdjTag,f.v,f.s,&err);dICHK(mesh,err);
    MeshListFree(r); MeshListFree(f); MeshListFree(rfo); MeshListFree(rv); MeshListFree(rvo);
    MeshListFree(fv); MeshListFree(fvo);
    err = PetscFree(rfOrient);dCHK(err);
  }

  /* Create the boundary parent set. */
  if (do_bdy)
  {
    const int nbdy = 3;
    const char bdyName[3][dNAME_LEN] = {dBDY_ROOT,"WALL","LID"};
    char normalName[dNAME_LEN];
    OnBdyFunc *onBdy[] = {0,on_wall,on_lid};/* We should never call onBdy[0], we use union instead. */
    iBase_EntitySetHandle bdy[3];
    iBase_TagHandle nameTag,bdyNormal[3]=MLZ;
    MeshListEH b=MLZ;
    MeshListInt off=MLZ,top=MLZ;
    dInt type,onbdy,flg;
    double *xx,y[3];
    char normal;

    /* Create entity sets to hold boundary facets, all inherit from the root boundary set. */
    for (i=0; i<nbdy; i++) {
      iMesh_createEntSet(mesh,0,&bdy[i],&err);dICHK(mesh,err);
      if (i > 0) {
        iMesh_addPrntChld(mesh,&bdy[0],&bdy[i],&err);dICHK(mesh,err);
      }
    }
    /* Set name tags over all the boundary entity sets */
    iMesh_createTag(mesh,dENT_SET_NAME,dNAME_LEN,iBase_BYTES,&nameTag,&err,strlen(dENT_SET_NAME));dICHK(mesh,err);
    for (i=0; i<nbdy; i++) {
      iMesh_setEntSetData(mesh,bdy[i],nameTag,bdyName[i],(int)strlen(bdyName[i]),&err);dICHK(mesh,err);
      if (i > 0) { /* For each set other than root, create a tag indicating whether face normal points in or out of the region. */
        err = PetscSNPrintf(normalName,sizeof(normalName),"%s+%s",dTAG_BDY_NORMAL,bdyName[i]);dCHK(err);
        iMesh_createTag(mesh,normalName,1,iBase_BYTES,&bdyNormal[i],&err,(int)strlen(normalName));dICHK(mesh,err);
      }
    }
    /* Put the boundary vertices in the appropriate sets, defined by their coordinates */
    iMesh_getEntities(mesh,0,iBase_VERTEX,iMesh_POINT,&v.v,&v.a,&v.s,&err);dICHK(mesh,err);
    iMesh_getVtxArrCoords(mesh,v.v,v.s,&order,&x.v,&x.a,&x.s,&err);dICHK(mesh,err);
    if (3*v.s != x.s) dERROR(1,"Wrong number of coordinates returned.");
    MeshListMalloc(v.s,b);
    for (i=1; i<nbdy; i++) {  /* All but the boundary root set */
      b.s = 0;
      for (j=0; j<v.s; j++) { /* All the vertices */
        if (onBdy[i](&box,&x.v[3*j])) b.v[b.s++] = v.v[j];
      }
      iMesh_addEntArrToSet(mesh,b.v,b.s,&bdy[i],&err);dICHK(mesh,err);
    }
    MeshListFree(b); MeshListFree(v); MeshListFree(x);
    /* Get entities of each dimension and adjacent vertices, check the vertices for membership, add entities to appropriate boundary sets */
    for (type=iBase_EDGE; type<=iBase_REGION; type++) {
      iMesh_getEntities(mesh,0,type,iMesh_ALL_TOPOLOGIES,&e.v,&e.a,&e.s,&err);dICHK(mesh,err);
      iMesh_getEntArrAdj(mesh,e.v,e.s,iBase_VERTEX,&v.v,&v.a,&v.s,&off.v,&off.a,&off.s,&err);dICHK(mesh,err);
      iMesh_getVtxArrCoords(mesh,v.v,v.s,&order,&x.v,&x.a,&x.s,&err);dICHK(mesh,err);
      iMesh_getEntArrTopo(mesh,e.v,e.s,&top.v,&top.a,&top.s,&err);dICHK(mesh,err);
      for (i=1; i<nbdy; i++) { /* all boundaries but root */
        for (j=0; j<e.s; j++) {   /* for each entity */
          onbdy = 1;
          /* An entity is in the set if all its vertices are in the set */
          for (k=off.v[j]; k<off.v[j+1]; k++) { /* for each vertex adjacent to the entity */
            iMesh_isEntContained(mesh,bdy[i],v.v[k],&flg,&err);dICHK(mesh,err);
            onbdy = onbdy && flg;
          }
          if (onbdy)  {
            iMesh_addEntToSet(mesh,e.v[j],&bdy[i],&err);dICHK(mesh,err);
            if (top.v[j] == iMesh_QUADRILATERAL) { /* Set the normal indicator for this boundary.  Normal=0 if face normal points outward. */
              xx = &x.v[off.v[j]*3]; /* The base for points of this face */
              dGeomVecMeanI(4,xx,y);  /* coordinate of center of face */
              normal = (char)!dGeomQuadParallel(xx,y); /* If the normal is parallel to the center vector, we say the normal faces out. */
              iMesh_setData(mesh,e.v[j],bdyNormal[i],&normal,1,&err);dICHK(mesh,err);
            }
          }
        }
      }
      MeshListFree(e); MeshListFree(v); MeshListFree(off); MeshListFree(x); MeshListFree(top);
    }
    /* add all the boundary sets to the root boundary set */
    for (i=1; i<nbdy; i++) {
      iMesh_addEntSet(mesh,bdy[i],&bdy[0],&err);dICHK(mesh,err);
    }
  }

                                /* Color the boundary entities */
  if (do_bdy)
  {
    iBase_TagHandle nameTag,bdyNum;
    iBase_EntitySetHandle bdyRoot=0;
    MeshListESH allESH=MLZ,bdy=MLZ;
    MeshListData name=MLZ;
    MeshListInt type=MLZ;
    dBool match;

    iMesh_getTagHandle(mesh,dENT_SET_NAME,&nameTag,&err,strlen(dENT_SET_NAME));dICHK(mesh,err);
    iMesh_getEntSets(mesh,0,0,&allESH.v,&allESH.a,&allESH.s,&err);dICHK(mesh,err);
    for (i=0; i<allESH.s; i++) {
      iMesh_getEntSetData(mesh,allESH.v[i],nameTag,&name.v,&name.a,&name.s,&err);dICHK(mesh,err);
      err = PetscStrcmp(dBDY_ROOT,name.v,&match);dCHK(err);
      /* printf("checking tag %d/%d: %s\n",i+1,allESH.s,name.v); */
      MeshListFree(name);
      if (match) {
        bdyRoot = allESH.v[i];
        MeshListFree(allESH);
        break;
      }
    }
    if (!bdyRoot) dERROR(1,"%s not found",dBDY_ROOT);
    iMesh_getChldn(mesh,bdyRoot,0,&bdy.v,&bdy.a,&bdy.s,&err);dICHK(mesh,err);
    iMesh_createTag(mesh,dTAG_BDY_NUM,1,iBase_INTEGER,&bdyNum,&err,strlen(dTAG_BDY_NUM));dICHK(mesh,err);
    for (i=0; i<bdy.s; i++) {
      iMesh_getEntSetData(mesh,bdy.v[i],nameTag,&name.v,&name.a,&name.s,&err);dICHK(mesh,err);
      err = dPrintf(PETSC_COMM_WORLD,"Found boundary %d/%d: %s\n",i+1,bdy.s,name.v);dCHK(err);
      MeshListFree(name);
      iMesh_setEntSetIntData(mesh,bdy.v[i],bdyNum,i,&err);dICHK(mesh,err);
      iMesh_getEntities(mesh,bdy.v[i],iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,&e.v,&e.a,&e.s,&err);dICHK(mesh,err);
      iMesh_getEntArrTopo(mesh,e.v,e.s,&type.v,&type.a,&type.s,&err);dICHK(mesh,err);
      for (j=0; j<e.s; j++) {
        if (type.v[j] >= iBase_VERTEX) { /* should we set the boundary on all entities or just faces? */
          iMesh_setIntData(mesh,e.v[j],bdyNum,i,&err);dICHK(mesh,err);
        }
      }
      MeshListFree(e);
    }
    MeshListFree(bdy);
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

  iMesh_save(mesh,0,outfile,outopts,&err,(int)sizeof(outfile),(int)strlen(outopts));dICHK(mesh,err);
  err = PetscFinalize();dCHK(err);
  dFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "dMeshCreateBoundary"
