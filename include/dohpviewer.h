#ifndef _DOHPVIEWER_H
#define _DOHPVIEWER_H

#include "dohptype.h"
#include "dohpfs.h"

#define PETSCVIEWERDHM "dhm"

dEXTERN_C_BEGIN

/* For writing */
extern dErr dViewerDHMSetTime(dViewer,dReal);
extern dErr dViewerDHMSetTimeUnits(dViewer,const char*,dReal);

/* For reading */
extern dErr dViewerDHMSetTimeStep(dViewer,dInt);
extern dErr dViewerDHMGetSteps(dViewer,dInt *nsteps,dReal **steptimes);
extern dErr dViewerDHMRestoreSteps(dViewer,dInt *nsteps,dReal **steptimes);

extern dErr dViewerRegisterAll(const char*);

/* Sort of a VisIt-specific structure at this point */
struct dViewerDHMSummaryFS {
  char  name[256];
  dInt  nblocks;
  dReal boundingbox[3][2];
};

struct dViewerDHMSummaryField {
  char name[256];
  char fsname[256];
  dInt bs;
};

dErr dViewerDHMGetStepSummary(PetscViewer viewer,dInt *nfs,const struct dViewerDHMSummaryFS **infs,dInt *nfields,const struct dViewerDHMSummaryField **infields);
dErr dViewerDHMRestoreStepSummary(PetscViewer viewer,dInt *nfs,const struct dViewerDHMSummaryFS **infs,dInt *nfields,const struct dViewerDHMSummaryField **infields);

dEXTERN_C_END

#endif
