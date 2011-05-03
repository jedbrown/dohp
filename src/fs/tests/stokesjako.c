#include "stokesimpl.h"
#include <ogr_srs_api.h>
#include <gdal.h>
#include <gdal_alg.h>
#include <gdalwarper.h>
#include <cpl_conv.h>   // CPLPrintPointer
#include <cpl_string.h> // CSLSetNameValue

// Most OGR functions return 0 for success
#define dOGRCHK(err) do {                                                 \
    if (PetscUnlikely((err) != OGRERR_NONE))                             \
      dERROR(PETSC_COMM_SELF,PETSC_ERR_LIB,"GDAL/OGR library routine failed with error code %d",(err)); \
  } while (0)

// Most GDAL functions return 1 for success
typedef int GDALErr;
#define dGDALCHK(err) do {                                              \
    if (PetscUnlikely((err) != 1))                                      \
      dERROR(PETSC_COMM_SELF,PETSC_ERR_LIB,"GDAL library routine failed with error code %d",(err)); \
  } while (0)

// Still other functions use this convention
#define dCPLCHK(err) do {                                               \
    if (PetscUnlikely((err) != CE_None))                                \
      dERROR(PETSC_COMM_SELF,PETSC_ERR_LIB,"GDAL/CPL library routine failed with error code %d",(err)); \
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
struct JakoInput {
  char path[PETSC_MAX_PATH_LEN];
};
typedef struct {
  OGRSpatialReferenceH utmref; // UTM zone 22N, the coordinates with Roman's surface elevation
  OGRSpatialReferenceH llref;  // Longitude-Latitude, for the CReSIS bed elevation
  OGRSpatialReferenceH ianref; // Ian Joughin's surface velocity
  OGRSpatialReferenceH myref;  // Computational domain
  OGRCoordinateTransformationH fromutm; // Convert coordinates from LonLat to current
  OGRCoordinateTransformationH fromll;  // Convert coordinates from LonLat to current
  OGRCoordinateTransformationH fromian; // Convert coordinates from Ian's projection to current
  struct JakoInput surface_elevation;
  struct JakoInput bed_elevation;
  struct JakoInput surface_velocity;
  double mygeo[6],myinvgeo[6];
  int nx,ny;
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
  oerr = OSRExportToPrettyWkt(ref,&wkt,0);dOGRCHK(oerr);
  err = PetscViewerASCIIPrintf(viewer,"WKT %s: %s\n\n",name,wkt);dCHK(err);
  OGRFree(wkt);
  dFunctionReturn(0);
}
static dErr JakoGDALMemAddBand(GDALDatasetH dset,GDALDataType dtype,void *memory)
{
  char buf[256] = {0},**bandoptions = NULL;
  int bytes,nx,ny;
  CPLErr cplerr;
  dErr err;

  dFunctionBegin;
  bytes = GDALGetDataTypeSize(dtype);
  nx = GDALGetRasterXSize(dset);
  ny = GDALGetRasterYSize(dset);
  err = dMalloc(nx*ny*bytes,(void**)memory);dCHK(err);

  // This is where the API moves from merely cumbersome to outright demeaning, like some twisted hazing ritual.
  CPLPrintPointer(buf,*(void**)memory,sizeof(buf));
  bandoptions = CSLSetNameValue(bandoptions,"DATAPOINTER",buf);
  cplerr = GDALAddBand(dset,dtype,bandoptions);dCPLCHK(cplerr);
  CSLDestroy(bandoptions);
  dFunctionReturn(0);
}
static dErr JakoGDALDatasetCreateMem(OGRSpatialReferenceH ref,const double geo[6],dInt n,dInt nlines,GDALDataType dtype,GDALDatasetH *dset,void *bandmem)
{
  char *wkt;
  GDALDriverH memdriver;
  CPLErr cplerr;
  OGRErr oerr;
  dErr err;

  dFunctionBegin;
  oerr = OSRExportToWkt(ref,&wkt);dOGRCHK(oerr);
  memdriver = GDALGetDriverByName("MEM");
  *dset = GDALCreate(memdriver,"MEM:::",n,nlines,0,dtype,NULL);
  cplerr = GDALSetProjection(*dset,wkt);dCPLCHK(cplerr);
  cplerr = GDALSetGeoTransform(*dset,(double*)geo);dCPLCHK(cplerr); /* const-incorrect interface */
  OGRFree(wkt);
  if (bandmem) {err = JakoGDALMemAddBand(*dset,GDT_Float64,&bandmem);dCHK(err);}
  dFunctionReturn(0);
}
static dErr JakoFileDataView(GDALDatasetH filedata,const char *name,PetscViewer viewer)
{
  dErr err;
  CPLErr cplerr;
  double geo[6],data[8*12];
  int nx=8,ny=12,snx,sny;
  GDALRasterBandH band;

  dFunctionBegin;
  cplerr = GDALGetGeoTransform(filedata,geo);
  err = dRealTableView(2,3,geo,PETSC_VIEWER_STDOUT_WORLD,name);dCHK(err);
  snx = GDALGetRasterXSize(filedata);
  sny = GDALGetRasterYSize(filedata);
  err = PetscViewerASCIIPrintf(viewer,"%s: nx=%d ny=%d\n",name,snx,sny);dCHK(err);
  band = GDALGetRasterBand(filedata,1);
  cplerr = GDALRasterIO(band,GF_Read,snx/2,sny/2,nx,ny,data,nx,ny,GDT_Float64,0,0);dCPLCHK(cplerr);
  err = dRealTableView(ny,nx,data,PETSC_VIEWER_STDOUT_WORLD,name);dCHK(err);
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
  oerr = OSRSetProjCS(jako->utmref,"UTM 22N (WGS84)");dOGRCHK(oerr);
  oerr = OSRSetWellKnownGeogCS(jako->utmref,"WGS84");dOGRCHK(oerr);
  oerr = OSRSetUTM(jako->utmref,22,1);dOGRCHK(oerr);
  oerr = OSRSetLinearUnits(jako->utmref,SRS_UL_METER,1.0);dOGRCHK(oerr);

  jako->llref  = OSRNewSpatialReference(NULL);
  oerr = OSRSetWellKnownGeogCS(jako->llref,"WGS84");dOGRCHK(oerr);

  jako->ianref = OSRNewSpatialReference(NULL);
  oerr = OSRSetProjCS(jako->ianref,"Stereographic Greenland (WGS84)");dOGRCHK(oerr);
  oerr = OSRSetWellKnownGeogCS(jako->ianref,"WGS84");dOGRCHK(oerr);
  oerr = OSRSetStereographic(jako->ianref,70.0,-45.0,100.0,-217.75e3,-2302.0e3);dOGRCHK(oerr);

  // The computational domain is in UTM 22N
  jako->myref = OSRNewSpatialReference(NULL);
  oerr = OSRSetProjCS(jako->myref,"Computational Domain: UTM 22N (WGS84)");dOGRCHK(oerr);
  oerr = OSRSetWellKnownGeogCS(jako->myref,"WGS84");dOGRCHK(oerr);
  oerr = OSRSetUTM(jako->myref,22,1);dOGRCHK(oerr);
  oerr = OSRSetLinearUnits(jako->myref,SRS_UL_METER,1.0);dOGRCHK(oerr);

  jako->fromutm = OCTNewCoordinateTransformation(jako->utmref,jako->myref);
  jako->fromll  = OCTNewCoordinateTransformation(jako->llref,jako->myref);
  jako->fromian = OCTNewCoordinateTransformation(jako->ianref,jako->myref);

  err = StokesCaseView_Jako(scase,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);

  if (1) {
    double x[] = {-49.858421420039,-48.4344803043268}, y[] = {69.6629932261656,69.4707311345926}, z[] = {373.8359,155.0298};
    int n = 2;
    for (dInt i=0; i<n; i++) printf("lonlat[%d] = %f %f %f\n",i,x[i],y[i],z[i]);
    oerr = !OCTTransform(jako->fromll,n,x,y,z);dOGRCHK(oerr); // This function has the opposite convention for error codes, but is in the OGR package
    for (dInt i=0; i<n; i++) printf("utm22n[%d] = %f %f %f\n",i,x[i],y[i],z[i]);
  }

  if (1) {
    double *hmem,*bmem,x0,y0,Lx,Ly,dx,dy;
    dReal *data;
    GDALDatasetH filedata,mysurf,mybed;
    GDALErr gerr;
    GDALRasterBandH band;
    CPLErr cplerr;

    jako->nx = 8;
    jako->ny = 12;
    x0 = scase->bbox[0][0];
    y0 = scase->bbox[1][0];
    Lx = scase->bbox[0][1] - x0;
    Ly = scase->bbox[1][1] - y0;
    // It appears from gdalwarp.cpp that pixels are cell centered
    dx = Lx / jako->nx;
    dy = Ly / jako->ny;

    // Convert a pixel to a physical coordinate in the current projection
    jako->mygeo[0] = x0 + dx/2;
    jako->mygeo[1] = dx;
    jako->mygeo[2] = 0;
    jako->mygeo[3] = y0 + Lx - dy/2; // Physical coordinates are "y up", pixels are "y down"
    jako->mygeo[4] = 0;
    jako->mygeo[5] = -dy;

    // Convert a physical coordinate in the current projection to a pixel coordinate
    gerr = GDALInvGeoTransform(jako->mygeo,jako->myinvgeo);dGDALCHK(gerr);
    err = dRealTableView(3,2,&scase->bbox[0][0],PETSC_VIEWER_STDOUT_WORLD,"bbox");dCHK(err);
    err = dRealTableView(2,3,jako->mygeo,PETSC_VIEWER_STDOUT_WORLD,"mygeo");dCHK(err);

    err = dMallocA(jako->nx*jako->ny,&data);dCHK(err);

    err = JakoGDALDatasetCreateMem(jako->myref,jako->mygeo,jako->nx,jako->ny,GDT_Float64,&mysurf,&hmem);dCHK(err);
    err = JakoGDALDatasetCreateMem(jako->myref,jako->mygeo,jako->nx,jako->ny,GDT_Float64,&mybed,&bmem);dCHK(err);

    filedata = GDALOpen(jako->surface_elevation.path, GA_ReadOnly); if (!filedata) dERROR(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"GDALOpen(\"%s\")",jako->surface_elevation.path);
    err = JakoFileDataView(filedata,"surface_elevation",PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
    cplerr = GDALReprojectImage(filedata,NULL,mysurf,NULL,GRA_Bilinear,0.0,0.0,GDALTermProgress,NULL,NULL);dCPLCHK(cplerr);
    GDALClose(filedata);

    band = GDALGetRasterBand(mysurf,1); // 1-based indexing
    cplerr = GDALRasterIO(band,GF_Read,0,0,jako->nx,jako->ny,data,jako->nx,jako->ny,GDT_Float64,0,0);dCPLCHK(cplerr);
    err = dRealTableView(jako->ny,jako->nx,data,PETSC_VIEWER_STDOUT_WORLD,"mysurf");dCHK(err);

    filedata = GDALOpen(jako->bed_elevation.path, GA_ReadOnly); if (!filedata) dERROR(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"GDALOpen(\"%s\")",jako->bed_elevation.path);
    err = JakoFileDataView(filedata,"bed_elevation",PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
    cplerr = GDALReprojectImage(filedata,NULL,mybed,NULL,GRA_Bilinear,0.0,0.0,GDALTermProgress,NULL,NULL);dCPLCHK(cplerr);
    GDALClose(filedata);

    band = GDALGetRasterBand(mybed,1); // 1-based indexing
    cplerr = GDALRasterIO(band,GF_Read,0,0,jako->nx,jako->ny,data,jako->nx,jako->ny,GDT_Float64,0,0);dCPLCHK(cplerr);
    err = dRealTableView(jako->ny,jako->nx,data,PETSC_VIEWER_STDOUT_WORLD,"mysurf");dCHK(err);

    GDALClose(mysurf);
    GDALClose(mybed);
    err = dFree(hmem);dCHK(err);
    err = dFree(bmem);dCHK(err);
    err = dFree(data);dCHK(err);
  }

  dFunctionReturn(0);
}

static dErr StokesCaseSetFromOptions_Jako(StokesCase scase)
{
  StokesCase_Jako *jako = scase->data;
  dBool flg;
  dErr err;
  dFunctionBegin;
  err = PetscOptionsHead("StokesCase_Jako options");dCHK(err); {
    err = PetscOptionsString("-jako_surface_elevation","File to read surface elevation from (assume same projection as model, e.g. UTM)","",
                             jako->surface_elevation.path,jako->surface_elevation.path,PETSC_MAX_PATH_LEN,&flg);dCHK(err);
    if (!flg) dERROR(scase->comm,PETSC_ERR_USER,"User must provide surface elevation file with -jako_surface_elevation FILENAME");

    err = PetscOptionsString("-jako_bed_elevation","File to read bed elevation from (assume same projection as model, e.g. UTM)","",
                                 jako->bed_elevation.path,jako->bed_elevation.path,PETSC_MAX_PATH_LEN,&flg);dCHK(err);
    if (!flg) dERROR(scase->comm,PETSC_ERR_USER,"User must provide bed elevation file with -jako_bed_elevation FILENAME");

    err = PetscOptionsString("-jako_surface_velocity","File to read surface velocity from (assume same projection as model, e.g. UTM)","",
                             jako->surface_velocity.path,jako->surface_velocity.path,PETSC_MAX_PATH_LEN,&flg);dCHK(err);
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
  OSRDestroySpatialReference(jako->ianref); jako->ianref = NULL;
  OSRDestroySpatialReference(jako->myref);  jako->myref = NULL;

  OCTDestroyCoordinateTransformation(jako->fromutm); jako->fromutm = NULL;
  OCTDestroyCoordinateTransformation(jako->fromll);  jako->fromll = NULL;
  OCTDestroyCoordinateTransformation(jako->fromian); jako->fromian = NULL;

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
  GDALAllRegister();
  err = StokesCaseRegister("gravity",StokesCaseCreate_Gravity);dCHK(err);
  err = StokesCaseRegister("jako",StokesCaseCreate_Jako);dCHK(err);
  dFunctionReturn(0);
}
