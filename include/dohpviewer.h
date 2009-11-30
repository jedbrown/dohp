#ifndef _DOHPVIEWER_H
#define _DOHPVIEWER_H

#include "dohptype.h"
#include "dohpfs.h"

#define PETSC_VIEWER_DHM "dhm"

dEXTERN_C_BEGIN

/* For writing */
extern dErr dViewerDHMSetTime(PetscViewer,dReal);
extern dErr dViewerDHMSetTimeUnits(PetscViewer,const char*,dReal);

/* For reading */
extern dErr dViewerDHMSetTimeStep(PetscViewer,dInt);
extern dErr dViewerDHMGetSteps(PetscViewer,dInt *nsteps,dReal **steptimes);
extern dErr dViewerDHMRestoreSteps(PetscViewer,dInt *nsteps,dReal **steptimes);

extern dErr dViewerRegisterAll(const char*);

dEXTERN_C_END

#endif
