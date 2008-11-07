static const char help[] = "Test the construction of dFS objects and anisotropic propogation.\n";

#include "dohpfs.h"
#include "dohp.h"
#include "dohpmesh.h"

#define ALEN(a) (dInt)(sizeof(a)/sizeof((a)[0]))

struct Options {
  dTruth constant_degree;
} gopt;

static dErr createHexMesh(iMesh_Instance mi)
{
  double vtx[12*3] = {0,0,0, 1,0,0, 1,1,0, 0,1,0,
                      0,0,1, 1,0,1, 1,1,1, 0,1,1,
                      2,0,0, 2,1,0, 2,1,1, 2,1,0};
  int rconn[16] = {0,1,2,3,4,5,6,7, 1,2,6,5,8,9,10,11};
  int fconn[11*4] = {0,1,2,3, 1,2,6,5, 2,3,7,6, 0,3,7,4, 0,1,5,4, 5,6,7,4,
                     8,11,5,1, 2,1,8,9, 8,9,10,11, 11,10,6,5, 9,10,6,2};
  int econn[20*2] = {0,1,4,5,7,6,3,2, 1,8,5,11,6,10,2,9,
                     0,3,3,7,7,4,4,0, 1,2,2,6,6,5,5,1, 8,11,11,10,10,9,9,8};
  iBase_EntityHandle work[100];
  MeshListEH v=MLZ,e=MLZ,f=MLZ,r=MLZ,tv=MLZ;
  MeshListInt stat=MLZ,off=MLZ;
  dIInt ierr;

  dFunctionBegin;
  iMesh_createVtxArr(mi,ALEN(vtx)/3,iBase_INTERLEAVED,vtx,ALEN(vtx),MLREF(v),&ierr);dICHK(mi,ierr);
  for (int i=0; i<ALEN(rconn); i++) work[i] = v.v[rconn[i]];
  iMesh_createEntArr(mi,iMesh_HEXAHEDRON,work,ALEN(rconn),MLREF(r),MLREF(stat),&ierr);dICHK(mi,ierr);
  {                             /* Check to see if any orientations changed */
    iMesh_getEntArrAdj(mi,r.v,r.s,iMesh_POINT,MLREF(tv),MLREF(off),&ierr);dICHK(mi,ierr);
    if (tv.s != v.s + 4) dERROR(1,"wrong number of vertices returned"); /* interface verts counted twice */
    for (int i=0; i<tv.s; i++) {
      if (v.v[rconn[i]] != tv.v[i]) dERROR(1,"unexpected vertex ordering");
    }
    MeshListFree(tv);
    MeshListFree(off);
  }
  MeshListFree(r); MeshListFree(stat);
  for (int i=0; i<ALEN(fconn); i++) work[i] = v.v[fconn[i]];
  iMesh_createEntArr(mi,iMesh_QUADRILATERAL,work,ALEN(fconn),MLREF(f),MLREF(stat),&ierr);dICHK(mi,ierr);
  MeshListFree(f); MeshListFree(stat);
  for (int i=0; i<ALEN(econn); i++) work[i] = v.v[econn[i]];
  iMesh_createEntArr(mi,iMesh_LINE_SEGMENT,work,ALEN(econn),MLREF(e),MLREF(stat),&ierr);dICHK(mi,ierr);
  MeshListFree(e); MeshListFree(stat);
  MeshListFree(v);
  dFunctionReturn(0);
}

/* Verify the adjacencies against these hand-verified values, should prevent regressions */
static dErr verifyAdjacencies(dMesh mesh)
{
  const dInt nents = 45,toff[5] = {0,12,32,43,45};
  const dInt adjoff[46] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,2,4, 6,8,10,12,14, 16,18,20,22,24, 26,28,30,32,34,
                           36,38,40,44,48, 52,56,60,64,68, 72,76,80,84,90, 96};
  const dInt adjind[96] = {0,1,4,5,7, 6,3,2,1,8, 5,11,6,10,2, 9,0,3,3,7,
                           7,4,4,0,1, 2,2,6,6,5, 5,1,8,11,11, 10,10,9,9,8,
                           12,24,15,20,24, 25,26,27,15,21, 14,25,20,21,22, 23,12,27,13,23,
                           26,14,22,13,28, 17,27,16,24,16, 31,19,31,30,29, 28,29,18,26,17,
                           30,18,25,19,36, 33,34,35,32,37, 39,42,41,38,33, 40};
  const dInt adjperm[96] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
                            0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
                            0,0,1,1,0, 0,0,0,1,0, 0,1,0,0,0, 0,0,1,1,0,
                            1,1,0,0,0, 1,0,0,1,0, 1,1,1,1,1, 1,0,1,0,0,
                            1,1,1,0,0, 0,0,5,4,1, 5,7,2,2,4, 0};
  struct dMeshAdjacency ma;
  dInt i;
  dErr err;

  dFunctionBegin;
  err = dMeshGetAdjacency(mesh,0,&ma);dCHK(err);
  dASSERT(ma.nents == nents);
  for (i=0; i<ALEN(toff); i++) dASSERT(ma.toff[i] == toff[i]);
  for (i=0; i<ALEN(adjoff); i++) dASSERT(ma.adjoff[i] == adjoff[i]);
  for (i=0; i<ALEN(adjind); i++) dASSERT(ma.adjind[i] == adjind[i]);
  for (i=0; i<ALEN(adjperm); i++) dASSERT(ma.adjperm[i] == adjperm[i]);
  err = dMeshRestoreAdjacency(mesh,0,&ma);dCHK(err);
  dFunctionReturn(0);
}

static dErr tagHexes(dMesh mesh,dMeshTag *intag)
{
  dMeshTag tag;
  dMeshEH ents[100],hex[2];
  dInt adata[3*ALEN(ents)],nonconst_data[3*ALEN(hex)] = {3,4,5,6,7,8},const_data[3*ALEN(hex)] = {4,4,4,4,4,4},*data,nhex,nents;
  dErr err;

  dFunctionBegin;
  *intag = 0;
  data = (gopt.constant_degree) ? const_data : nonconst_data;
  err = dMeshTagCreateTemp(mesh,"anisotropic",3,dDATA_INT,&tag);dCHK(err);
  /* tag edges and faces with high values, currently needed to propogate degrees (will overwrite high values) */
  for (dEntType type=dTYPE_VERTEX; type<=dTYPE_FACE; type++) {
    err = dMeshGetEnts(mesh,0,type,dTOPO_ALL,ents,ALEN(ents),&nents);
    for (dInt i=0; i<3*nents; i++) adata[i] = 30;
    err = dMeshTagSetData(mesh,tag,ents,nents,adata,3*nents,dDATA_INT);dCHK(err);
  }
  /* tag the hexes with meaningful values */
  err = dMeshGetEnts(mesh,0,dTYPE_ALL,dTOPO_HEX,hex,ALEN(hex),&nhex);dCHK(err);
  if (nhex != 2) dERROR(1,"wrong number of hexes");
  err = dMeshTagSetData(mesh,tag,hex,nhex,data,3*nhex,dDATA_INT);dCHK(err);
  *intag = tag;
  dFunctionReturn(0);
}

static dErr examine(dMesh mesh,dMeshTag tag)
{
  dMeshEH ents[100];
  dEntTopology topo[100];
  MeshListEH conn=MLZ;
  MeshListInt connoff=MLZ;
  dInt data[3*ALEN(ents)],nents,nv;
  dIInt ierr;
  iMesh_Instance mi;
  dErr err;

  dFunctionBegin;
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  for (dEntType type=dTYPE_EDGE; type<dTYPE_ALL; type++) {
    err = dMeshGetEnts(mesh,0,type,dTOPO_ALL,ents,ALEN(ents),&nents);dCHK(err);
    err = dMeshGetTopo(mesh,nents,ents,topo);dCHK(err);
    iMesh_getEntArrAdj(mi,ents,nents,iBase_VERTEX,MLREF(conn),MLREF(connoff),&ierr);dICHK(mi,ierr);
    err = dMeshTagGetData(mesh,tag,ents,nents,data,ALEN(data),dDATA_INT);dCHK(err);
    for (dInt i=0; i<nents; i++) {
      switch (topo[i]) {
        case dTOPO_HEX: nv = 8; break;
        case dTOPO_QUAD: nv = 4; break;
        case dTOPO_LINE: nv = 2; break;
        default: dERROR(1,"unsupported topology");
      }
      err = dPrintf(PETSC_COMM_WORLD,"%30s %d %d %d   [",iMesh_TopologyName[topo[i]],data[3*i],data[3*i+1],data[3*i+2]);dCHK(err);
      for (dInt j=0; j<nv; j++) { err = dPrintf(PETSC_COMM_WORLD," %p",conn.v[connoff.v[i]+j]);dCHK(err); }
      err = dPrintf(PETSC_COMM_WORLD," ]\n");dCHK(err);
    }
    MeshListFree(conn); MeshListFree(connoff);
  }
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  /* const char pTagName[] = "dohp_partition"; */
  iMesh_Instance mi;
  dJacobi jac;
  dFS fs;
  dMesh mesh;
  dMeshESH domain;
  dMeshTag rtag,dtag;
  MPI_Comm comm;
  PetscViewer viewer;
  dErr err;

  dFunctionBegin;
  err = PetscInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  err = PetscOptionsBegin(comm,NULL,"Test options","ex1");dCHK(err); {
    gopt.constant_degree = dFALSE;
    err = PetscOptionsTruth("-constant_degree","Use constant isotropic degree on all elements",NULL,gopt.constant_degree,&gopt.constant_degree,NULL);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);
  err = dMeshCreate(comm,&mesh);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshGetInstance(mesh,&mi);dCHK(err);
  err = createHexMesh(mi);dCHK(err);
  err = dMeshView(mesh,viewer);dCHK(err);
  err = verifyAdjacencies(mesh);dCHK(err);
  iMesh_getRootSet(mi,&domain,&err);dICHK(mi,err);

  err = dJacobiCreate(comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);
  err = dJacobiSetUp(jac);dCHK(err);

  err = dMeshCreateRuleTagIsotropic(mesh,domain,jac,"ex1_rule",5,&rtag);dCHK(err);
  err = tagHexes(mesh,&dtag);dCHK(err);

  err = dFSCreate(comm,&fs);dCHK(err);
  err = dFSSetMesh(fs,mesh,domain);dCHK(err);
  err = dFSSetRuleTag(fs,jac,rtag);dCHK(err);
  err = dFSSetDegree(fs,jac,dtag);dCHK(err);
  err = dFSSetFromOptions(fs);dCHK(err);

  err = examine(mesh,dtag);dCHK(err);

  err = dFSDestroy(fs);dCHK(err);
  err = dJacobiDestroy(jac);dCHK(err);
  err = dMeshDestroy(mesh);
  err = PetscFinalize();dCHK(err);
  dFunctionReturn(0);
}
