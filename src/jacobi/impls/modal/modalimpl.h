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

KHASH_MAP_INIT_INT(modal,ModalBasis)

typedef struct {
  dEFSHEADER;
  ModalBasis basis;
  void *unused[2];
} dEFS_Modal;

typedef struct {
  khash_t(modal) *topo[dTOPO_ALL];
  dJacobiModalFamily family;
  struct _dEFSOps *efsOpsLine,*efsOpsQuad,*efsOpsHex;
} dJacobi_Modal;

extern dErr ModalBasisView(ModalBasis,PetscViewer);
extern dErr dJacobiEFSOpsSetUp_Modal(dJacobi);

#endif
