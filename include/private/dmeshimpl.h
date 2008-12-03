#ifndef _DMESHIMPL_H
#define _DMESHIMPL_H

#include "dohpmesh.h"

struct _dMeshOps {
  dErr (*orientfacets)(dMesh);
  dErr (*tagbcast)(dMesh,dMeshTag);
  dErr (*setfromoptions)(dMesh);
  dErr (*load)(dMesh);
  dErr (*view)(dMesh,PetscViewer);
  dErr (*destroy)(dMesh);
};

struct _p_dMesh {
  PETSCHEADER(struct _dMeshOps);
  char *infile,*inoptions;
  iMesh_Instance mi;
  iBase_EntitySetHandle root;
  MeshListEH v,e,f,r;           /* vertices, edges, faces, vertices */
  MeshListEH arf,afe,aev;       /* adjacencies region -> face -> edge -> vertex */
  MeshListData orf,ofe;
  MeshListReal x;
  void *data;
  dMeshTag manifoldTag,manifoldOrientTag;
  dMeshManifold *manifoldList;
  dInt nManifolds;
};

struct _p_dMeshManifold {
  char name[dNAME_LEN];
  dInt toff[5];
  dMeshEH *ents;
  char *orient;
};

#endif
