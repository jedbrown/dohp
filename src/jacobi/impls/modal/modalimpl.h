#if !defined _MODALIMPL_H
#define _MODALIMPL_H

#include <dohpjacimpl.h>
#include <dohpkhash.h>

typedef struct _ModalBasis *ModalBasis;
struct _ModalBasis {
  dInt P;                       /* number of modes */
  dInt Q;                       /* number of quadrature points */
  dInt dim;                     /* spatial dimension */
  dReal *interp;
  dReal *deriv;
};

typedef struct {
  dEFSHEADER;
  dEntTopology topo;
  ModalBasis basis;
} dEFS_Modal;

KHASH_MAP_INIT_INT(modal,ModalBasis)
KHASH_MAP_INIT_INT64(efs,dEFS_Modal*)

typedef struct {
  khash_t(modal) *topo[dTOPO_ALL];
  khash_t(efs) *efs;
  dJacobiModalFamily family;
  struct _dEFSOps *efsOpsLine,*efsOpsQuad,*efsOpsHex;
} dJacobi_Modal;

extern dErr ModalBasisView(ModalBasis,PetscViewer);
extern dErr dJacobiEFSOpsSetUp_Modal(dJacobi);

#endif
