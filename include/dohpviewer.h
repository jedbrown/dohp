#ifndef _DOHPVIEWER_H
#define _DOHPVIEWER_H

#include "dohptype.h"
#include "dohpfs.h"

#define PETSC_VIEWER_DHM "dhm"

PETSC_EXTERN_CXX_BEGIN

EXTERN dErr dViewerDHMSetFS(PetscViewer,dFS);
EXTERN dErr dViewerDHMSetTime(PetscViewer,dReal);
EXTERN dErr dViewerRegisterAll(const char*);

PETSC_EXTERN_CXX_END

#endif
