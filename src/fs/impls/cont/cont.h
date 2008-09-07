#ifndef _CONT_H
#define _CONT_H

#include "private/fsimpl.h"

PETSC_EXTERN_CXX_BEGIN

dErr dFSCreate_Cont(dFS);

typedef struct {
  MeshListEH r,f,e,v;
} dFS_Cont;

PETSC_EXTERN_CXX_END

#endif  /* _CONT_H */
