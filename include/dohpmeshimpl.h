#ifndef _DMESHIMPL_H
#define _DMESHIMPL_H

#include "dohpmesh.h"

struct _dMeshOps {
  dErr (*orientfacets)(dMesh);
  dErr (*tagbcast)(dMesh,dMeshTag);
  dErr (*setfromoptions)(dMesh);
  dErr (*load)(dMesh);
  dErr (*view)(dMesh,dViewer);
  dErr (*destroy)(dMesh);
};

struct _p_dMesh {
  PETSCHEADER(struct _dMeshOps);
  char *infile,*inoptions;
  iMesh_Instance mi;
#ifdef dHAVE_ITAPS_REL
  iGeom_Instance igeom;
  iRel_Instance irel;
#endif
  dMeshESH root,emptyset;
  dMeshTag senseTag;
  MeshListEH v,e,f,r;           /* vertices, edges, faces, vertices */
  MeshListEH arf,afe,aev;       /* adjacencies region -> face -> edge -> vertex */
  MeshListData orf,ofe;
  MeshListReal x;
  void *data;
};

#endif
