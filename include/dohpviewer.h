#ifndef _DOHPVIEWER_H
#define _DOHPVIEWER_H

#include "dohptype.h"
#include "dohpfs.h"

#define PETSC_VIEWER_DHM "dhm"

dEXTERN_C_BEGIN

/* For writing */
extern dErr dViewerDHMSetTime(dViewer,dReal);
extern dErr dViewerDHMSetTimeUnits(dViewer,const char*,dReal);

/* For reading */
extern dErr dViewerDHMSetTimeStep(dViewer,dInt);
extern dErr dViewerDHMGetSteps(dViewer,dInt *nsteps,dReal **steptimes);
extern dErr dViewerDHMRestoreSteps(dViewer,dInt *nsteps,dReal **steptimes);

extern dErr dViewerRegisterAll(const char*);

dEXTERN_C_END

#endif
