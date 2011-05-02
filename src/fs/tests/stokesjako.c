#include "stokesimpl.h"
#include <ogr_srs_api.h>

#define dOCHK(err) do {                                                 \
    if (PetscUnlikely((err) != OGRERR_NONE))                             \
      dERROR(PETSC_COMM_SELF,PETSC_ERR_LIB,"OGR library routine failed with error code %d",err); \
  } while (0)

// Trivial gravity model
static void StokesCaseSolution_Gravity(StokesCase dUNUSED scase,const dReal dUNUSED x[3],dScalar u[],dScalar du[],dScalar *p,dScalar dp[])
{                               /* Defines inhomogeneous Dirichlet boundary conditions */
  u[0] = u[1] = u[2] = 0;
  for (dInt i=0; i<9; i++) du[i] = 0;
  *p = 0;
  for (dInt i=0; i<3; i++) dp[i] = 0;
}
static void StokesCaseForcing_Gravity(StokesCase scase,const dReal dUNUSED x[3],dScalar fu[],dScalar *fp)
{
  fu[0] = 0;
  fu[1] = 0;
  fu[2] = scase->gravity;
  fp[0] = 0;
}
static dErr StokesCaseCreate_Gravity(StokesCase scase)
{
  dFunctionBegin;
  scase->reality = dTRUE;
  scase->solution = StokesCaseSolution_Gravity;
  scase->forcing  = StokesCaseForcing_Gravity;
  dFunctionReturn(0);
}

// A real implementation
typedef struct {
  OGRSpatialReferenceH utmref; // UTM zone 22N, the coordinates with Roman's surface elevation
  OGRSpatialReferenceH llref;  // Longitude-Latitude, for the CReSIS bed elevation
  OGRSpatialReferenceH ianref; // Ian Joughin's surface velocity
  OGRCoordinateTransformationH fromll;  // Convert from LonLat to current
  OGRCoordinateTransformationH fromian; // Convert from Ian's projection to current
} StokesCase_Jako;

static dErr StokesCaseSolution_Jako(StokesCase scase,const dReal x[3],dScalar u[],dScalar du[],dScalar *p,dScalar dp[])
{                               /* Defines inhomogeneous Dirichlet boundary conditions */
  dUNUSED StokesCase_Jako *jako = scase->data;

  dFunctionBegin;
  if (x[0] > 580000.) {
    u[0] = -200;
    u[1] = -100;
    u[2] = 0;
  } else {
    u[0] = -400;
    u[1] = -300;
    u[2] = 0;
  }

  for (dInt i=0; i<9; i++) du[i] = 0;
  *p = 0;
  for (dInt i=0; i<3; i++) dp[i] = 0;
  dFunctionReturn(0);
}
static void StokesCaseSolution_Jako_Void(StokesCase scase,const dReal x[3],dScalar u[],dScalar du[],dScalar *p,dScalar dp[])
{
  dErr err;
  err = StokesCaseSolution_Jako(scase,x,u,du,p,dp);CHKERRV(err);
}
static dErr JakoViewWKT(OGRSpatialReferenceH ref,const char *name,PetscViewer viewer)
{
  dErr err;
  OGRErr oerr;
  char *wkt;

  dFunctionBegin;
  oerr = OSRExportToPrettyWkt(ref,&wkt,0);dOCHK(oerr);
  err = PetscViewerASCIIPrintf(viewer,"WKT %s: %s\n\n",name,wkt);dCHK(err);
  OGRFree(wkt);
  dFunctionReturn(0);
}
static dErr StokesCaseView_Jako(StokesCase scase,PetscViewer viewer)
{
  StokesCase_Jako *jako = scase->data;
  dErr err;

  dFunctionBegin;
  err = PetscViewerASCIIPrintf(viewer,"StokesCase: Jakobshavn\n");dCHK(err);
  err = JakoViewWKT(jako->utmref,"UTM",viewer);dCHK(err);
  err = JakoViewWKT(jako->llref,"LonLat",viewer);dCHK(err);
  err = JakoViewWKT(jako->ianref,"Ian",viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseSetUp_Jako(StokesCase scase)
{
  StokesCase_Jako *jako = scase->data;
  dErr err;
  OGRErr oerr;

  dFunctionBegin;
  jako->utmref = OSRNewSpatialReference(NULL);
  oerr = OSRSetProjCS(jako->utmref,"UTM 22N (WGS84)");dOCHK(oerr);
  oerr = OSRSetWellKnownGeogCS(jako->utmref,"WGS84");dOCHK(oerr);
  oerr = OSRSetUTM(jako->utmref,22,1);dOCHK(oerr);

  jako->llref  = OSRNewSpatialReference(NULL);
  oerr = OSRSetWellKnownGeogCS(jako->llref,"WGS84");dOCHK(oerr);

  jako->ianref = OSRNewSpatialReference(NULL);
  oerr = OSRSetProjCS(jako->ianref,"Stereographic Greenland (WGS84)");dOCHK(oerr);
  oerr = OSRSetWellKnownGeogCS(jako->ianref,"WGS84");dOCHK(oerr);
  oerr = OSRSetStereographic(jako->ianref,70.0,-45.0,100.0,-217.75e3,-2302.0e3);dOCHK(oerr);

  err = StokesCaseView_Jako(scase,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
  dFunctionReturn(0);
}

static dErr StokesCaseSetFromOptions_Jako(StokesCase scase)
{
  char fname[256] = "unknown";
  dBool flg;
  dErr err;
  dFunctionBegin;
  err = PetscOptionsHead("StokesCase_Jako options");dCHK(err); {
    err = PetscOptionsString("-jako_surface_velocity","File to read surface velocity from (assume same projection as model, e.g. UTM)","",fname,fname,sizeof(fname),&flg);dCHK(err);
    if (!flg) dERROR(scase->comm,PETSC_ERR_USER,"User must provide surface velocity file with -jako_surface_velocity FILENAME");
  } err = PetscOptionsTail();dCHK(err);
  err = StokesCaseSetUp_Jako(scase);dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseDestroy_Jako(StokesCase scase)
{
  StokesCase_Jako *jako = scase->data;
  dErr err;

  dFunctionBegin;
  OSRDestroySpatialReference(jako->utmref); jako->utmref = NULL;
  OSRDestroySpatialReference(jako->llref);  jako->llref  = NULL;
  err = dFree(scase->data);dCHK(err);
  dFunctionReturn(0);
}
static dErr StokesCaseCreate_Jako(StokesCase scase)
{
  StokesCase_Jako *jako;
  dErr err;

  dFunctionBegin;
  scase->reality = dTRUE;
  scase->solution = StokesCaseSolution_Jako_Void;
  scase->forcing  = StokesCaseForcing_Gravity;
  scase->setfromoptions = StokesCaseSetFromOptions_Jako;
  scase->destroy = StokesCaseDestroy_Jako;

  err = dNew(StokesCase_Jako,&jako);dCHK(err);
  scase->data = jako;
  dFunctionReturn(0);
}

dErr StokesCaseRegisterAll_Jako(void)
{
  dErr err;

  dFunctionBegin;
  err = StokesCaseRegisterAll_Exact();dCHK(err);
  err = StokesCaseRegister("gravity",StokesCaseCreate_Gravity);dCHK(err);
  err = StokesCaseRegister("jako",StokesCaseCreate_Jako);dCHK(err);
  dFunctionReturn(0);
}
