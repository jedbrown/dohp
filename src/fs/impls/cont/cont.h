#ifndef _CONT_H
#define _CONT_H

#include "private/fsimpl.h"

PETSC_EXTERN_CXX_BEGIN

extern dErr dFSCreate_Cont(dFS);

struct dFS_Cont {
  dBool usecmatrix;
  dInt m;
  Mat F;                        /**< facet constraint matrix */
  Mat Cprimal;                  /**< primal constraint matrix (element dofs to local numbering, full order) */
  Mat Cdual;                    /**< dual constraint matrix (element dofs to local numbering, as sparse as possible) */
};

PETSC_EXTERN_CXX_END

#endif  /* _CONT_H */
