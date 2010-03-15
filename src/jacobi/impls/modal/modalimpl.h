#if !defined _MODALIMPL_H
#define _MODALIMPL_H

#include <dohpjacimpl.h>
#include <dohpkhash.h>

typedef struct {
  dEFSHEADER;
  dEntTopology topo;
  dInt         P;               /* number of modes */
  dInt         Q;               /* number of quadrature points */
  dReal        *interp;
  dReal        *deriv;
} dEFS_Modal;

typedef struct {
  dEntTopology topo;
  dPolynomialOrder degree;
  dRule rule;
} khu_efskey_t;
static inline khint_t khu_efskey_hash_func(khu_efskey_t key)
{ return kh_int_hash_func((khint32_t)key.topo) ^ kh_int_hash_func((khint32_t)key.degree) ^ kh_int64_hash_func((khint64_t)(uintptr_t)key.rule); }
static inline bool khu_efskey_hash_equal(khu_efskey_t a,khu_efskey_t b)
{ return (a.topo == b.topo) && (a.degree == b.degree) && (a.rule == b.rule); }
KHASH_INIT(efs, khu_efskey_t, dEFS_Modal*, 1, khu_efskey_hash_func, khu_efskey_hash_equal)

typedef struct {
  khash_t(efs) *efs;
  dJacobiModalFamily family;
  struct _dEFSOps *efsOpsLine,*efsOpsQuad,*efsOpsHex;
} dJacobi_Modal;

extern dErr dJacobiEFSOpsSetUp_Modal(dJacobi);

#endif
