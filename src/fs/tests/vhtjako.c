#include "vhtimpl.h"
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

// Advection model
static dErr VHTCaseSolution_Wind(VHTCase dUNUSED scase,const dReal x[3],dScalar rhou[],dScalar drhou[],dScalar *p,dScalar dp[],dScalar *E,dScalar dE[])
{                               /* Defines inhomogeneous Dirichlet boundary conditions */
  dReal z = (1 + x[2])/2; // in range [0,1]
  dReal depth = 1-z;
  rhou[0] = (1 - dSqr(x[1])) * (1 - pow(depth,2));
  rhou[1] = 0;
  rhou[2] = 0;
  for (dInt i=0; i<9; i++) drhou[i] = 0;
  *p = 0;
  for (dInt i=0; i<3; i++) dp[i] = 0;
  *E = -x[0] * (1 - dSqr(x[1]));
  for (dInt i=0; i<3; i++) dE[i] = 0;
  return 0;
}
static dErr VHTCaseForcing_Wind(VHTCase dUNUSED scase,const dReal dUNUSED x[3],dScalar frhou[],dScalar *fp,dScalar *fE)
{
  for (dInt i=0; i<3; i++) frhou[i] = 0;
  fp[0] = 0;
  fE[0] = 0;
  return 0;
}
static dErr VHTCaseCreate_Wind(VHTCase scase)
{
  dFunctionBegin;
  scase->reality = dTRUE;
  scase->solution = VHTCaseSolution_Wind;
  scase->forcing  = VHTCaseForcing_Wind;
  dFunctionReturn(0);
}

// A real implementation
struct JakoInput {
  char path[PETSC_MAX_PATH_LEN];
  GDALDatasetH dataset;
  double *memory;
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
  double *h,*b;
  int nx,ny;
  dBool verbose;
} VHTCase_Jako;
static dErr JakoFindPixel(VHTCase scase,const dReal x[3],dInt *ix,dInt *iy)
{
  VHTCase_Jako *jako = scase->data;
  double xpixel,ypixel;
  dInt i,j;

  dFunctionBegin;
  *ix = -1;
  *iy = -1;
  if (x[0] < scase->bbox[0][0] || scase->bbox[0][1] < x[0])
    dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"x[0]=%f not in bounding box %f:%f",x[0],scase->bbox[0][0],scase->bbox[0][1]);
  if (x[1] < scase->bbox[1][0] || scase->bbox[1][1] < x[1])
    dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"x[1]=%f not in bounding box %f:%f",x[1],scase->bbox[1][0],scase->bbox[1][1]);
  GDALApplyGeoTransform(jako->myinvgeo,x[0],x[1],&xpixel,&ypixel);
  i = (dInt)xpixel;
  j = (dInt)ypixel;
  if (i < 0 || jako->nx <= i) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Computed xp=%D not in range %D:%D",i,0,jako->nx);
  if (j < 0 || jako->ny <= j) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Computed yp=%D not in range %D:%D",j,0,jako->ny);
  *ix = i;
  *iy = j;
  dFunctionReturn(0);
}
static dErr JakoGradient2(VHTCase scase,const dReal field[],dInt i,dInt j,dReal grad[2])
{
  VHTCase_Jako *jako = scase->data;
  const dReal dx = jako->mygeo[1],dy = jako->mygeo[5];
  const dInt nx = jako->nx,ny = jako->ny;
  dFunctionBegin;
  if (i+1 >= nx) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"X slope would evaluate point (%D,%D) outside of valid range %Dx%D",i+1,j,nx,ny);
  if (j-1 < 0) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Y slope would evaluate point (%D,%D) outside of valid range %Dx%D",i,j-1,nx,ny);
  grad[0] = (field[j*nx + i + 1] - field[j*nx + i]) / dx;
  grad[1] = (field[(j-1)*nx + i] - field[j*nx + i]) / dy; // The raster has y axis flipped
  dFunctionReturn(0);
}
static dErr JakoSIAVelocity(VHTCase scase,dReal b,dReal h,dReal dh[2],dReal z,dScalar u[3])
{
  struct VHTRheology *rheo = &scase->rheo;
  const dReal
    B = rheo->B0,
    p = rheo->pe,
    n = 1./(p-1),
    A = pow(B,-n);
  const dScalar
    slope = dSqrt(dSqr(dh[0]) + dSqr(dh[1])),
    sliding = 0,
    hmz = z > h ? 0 : (z < b ? h-b : h-z),
    siaspeed = sliding + A / (n+1) * pow(slope,n) * (pow(h-b,n+1) - pow(hmz,n+1)),
    speed = 0 + (siaspeed - 0) * (1 + tanh((z-b)/100))/2; // Integrate stress from the bottom
  u[0] = -speed / slope * dh[0];
  u[1] = -speed / slope * dh[1];
  u[2] = 0; // Would need another derivative to evaluate this and it should not be a big deal for this stationary computation
  return 0;
}
static dErr JakoInternalEnergy(VHTCase scase,dReal b,dReal h,const dReal x[],dScalar *e)
{
  struct VHTRheology *rheo = &scase->rheo;
  dReal z,t,T;

  dFunctionBegin;
  z = (x[2] - b) / (h-b); // Normalize to [0,1] plus possible fuzz
  z = dMax(0,dMin(1,z)); // Clip to [0,1]
  t = 90*dSqr(z)*exp(-5*z); // hump with maximum value close to 2
  T = rheo->T3 - (rheo->T3 - rheo->T0)*t; // Maximum value is at the bed and equal to T3, extends past T0
  *e = rheo->c_i * (T - rheo->T0);        // energy/mass
  dFunctionReturn(0);
}
static dErr VHTCaseSolution_Jako(VHTCase scase,const dReal x[3],dScalar rhou[],dScalar drhou[],dScalar *p,dScalar dp[],dScalar *E,dScalar dE[])
{                               /* Defines inhomogeneous Dirichlet boundary conditions */
  VHTCase_Jako *jako = scase->data;
  struct VHTRheology *rheo = &scase->rheo;
  double h,b;
  dInt pixel,xp,yp;
  dReal u[3],dh[2] = {0,0};
  dScalar e;
  dErr err;

  dFunctionBegin;
  err = JakoFindPixel(scase,x,&xp,&yp);dCHK(err);
  pixel = yp*jako->nx + xp;
  h = jako->h[pixel];
  b = jako->b[pixel];
  err = JakoGradient2(scase,jako->h,xp,yp,dh);dCHK(err);
  err = JakoSIAVelocity(scase,b,h,dh,x[2],u);dCHK(err);
  err = JakoInternalEnergy(scase,b,h,x,&e);dCHK(err);

  for (dInt i=0; i<3; i++) rhou[i] = rheo->rhoi * u[i]; // Just assume constant density
  for (dInt i=0; i<9; i++) drhou[i] = 0;
  *p = rheo->rhoi * rheo->gravity[2] * (h - x[2]);
  for (dInt i=0; i<3; i++) dp[i] = -rheo->rhoi * rheo->gravity[i];
  *E = rheo->rhoi * e;
  for (dInt i=0; i<3; i++) dE[i] = 0;
  dFunctionReturn(0);
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
  err = dRealTableView(2,3,geo,PETSC_VIEWER_STDOUT_WORLD,"%s:geo",name);dCHK(err);
  snx = GDALGetRasterXSize(filedata);
  sny = GDALGetRasterYSize(filedata);
  err = PetscViewerASCIIPrintf(viewer,"%s: nx=%d ny=%d\n",name,snx,sny);dCHK(err);
  band = GDALGetRasterBand(filedata,1);
  cplerr = GDALRasterIO(band,GF_Read,snx/2,sny/2,nx,ny,data,nx,ny,GDT_Float64,0,0);dCPLCHK(cplerr);
  err = dRealTableView(ny,nx,data,PETSC_VIEWER_STDOUT_WORLD,name);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTCaseView_Jako(VHTCase scase,PetscViewer viewer)
{
  VHTCase_Jako *jako = scase->data;
  dErr err;

  dFunctionBegin;
  err = PetscViewerASCIIPrintf(viewer,"VHTCase: Jakobshavn\n");dCHK(err);
  err = JakoViewWKT(jako->utmref,"UTM",viewer);dCHK(err);
  err = JakoViewWKT(jako->llref,"LonLat",viewer);dCHK(err);
  err = JakoViewWKT(jako->ianref,"Ian",viewer);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseSetUp_Jako(VHTCase scase)
{
  VHTCase_Jako *jako = scase->data;
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

  if (jako->verbose) {err = VHTCaseView_Jako(scase,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);}

  if (1) {
    GDALProgressFunc gdalprogress = jako->verbose ? GDALTermProgress : GDALDummyProgress;
    double x0,y0,Lx,Ly,dx,dy;
    GDALDatasetH filedata;
    GDALErr gerr;
    GDALRasterBandH band;
    CPLErr cplerr;

    x0 = scase->bbox[0][0];
    y0 = scase->bbox[1][0];
    Lx = scase->bbox[0][1] - x0;
    Ly = scase->bbox[1][1] - y0;
    // It appears from gdalwarp.cpp that pixels are cell centered
    dx = Lx / jako->nx;
    dy = Ly / jako->ny;

    // Extend the domain in the positive direction for slope evaluation
    jako->nx += 2;
    jako->ny += 2;
    Lx += 2*dx;
    Ly += 2*dy;

    // Convert a pixel to a physical coordinate in the current projection
    jako->mygeo[0] = x0 + dx/2;
    jako->mygeo[1] = dx;
    jako->mygeo[2] = 0;
    jako->mygeo[3] = y0 + Ly - dy/2; // Physical coordinates are "y up", pixels are "y down"
    jako->mygeo[4] = 0;
    jako->mygeo[5] = -dy;

    // Convert a physical coordinate in the current projection to a pixel coordinate
    gerr = GDALInvGeoTransform(jako->mygeo,jako->myinvgeo);dGDALCHK(gerr);
    if (0) {
      dInt i,j;
      double a,b;
      err = dRealTableView(3,2,&scase->bbox[0][0],PETSC_VIEWER_STDOUT_WORLD,"bbox");dCHK(err);
      err = dRealTableView(2,3,jako->mygeo,PETSC_VIEWER_STDOUT_WORLD,"mygeo");dCHK(err);
      GDALApplyGeoTransform(jako->mygeo,0,0,&a,&b); printf("geo[0,0] = %f,%f\n",a,b);
      GDALApplyGeoTransform(jako->mygeo,1,1,&a,&b); printf("geo[1,1] = %f,%f\n",a,b);
      GDALApplyGeoTransform(jako->mygeo,0.5,0.5,&a,&b); printf("geo[0.5,0.5] = %f,%f\n",a,b);
      err = JakoFindPixel(scase,(dReal[]){scase->bbox[0][0],scase->bbox[1][0]},&i,&j);dCHK(err);
      printf("xmin,ymin: (%d,%d)\n",i,j);
      err = JakoFindPixel(scase,(dReal[]){scase->bbox[0][1],scase->bbox[1][1]},&i,&j);dCHK(err);
      printf("xmax,ymax: (%d,%d)\n",i,j);
    }

    err = dMallocA2(jako->nx*jako->ny,&jako->h,jako->nx*jako->ny,&jako->b);dCHK(err);

    err = JakoGDALDatasetCreateMem(jako->myref,jako->mygeo,jako->nx,jako->ny,GDT_Float64,&jako->surface_elevation.dataset,&jako->surface_elevation.memory);dCHK(err);
    err = JakoGDALDatasetCreateMem(jako->myref,jako->mygeo,jako->nx,jako->ny,GDT_Float64,&jako->bed_elevation.dataset,&jako->bed_elevation.memory);dCHK(err);

    filedata = GDALOpen(jako->surface_elevation.path, GA_ReadOnly); if (!filedata) dERROR(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"GDALOpen(\"%s\")",jako->surface_elevation.path);
    if (jako->verbose) {err = JakoFileDataView(filedata,"surface_elevation",PETSC_VIEWER_STDOUT_WORLD);dCHK(err);}
    cplerr = GDALReprojectImage(filedata,NULL,jako->surface_elevation.dataset,NULL,GRA_Bilinear,0.0,0.0,gdalprogress,NULL,NULL);dCPLCHK(cplerr);
    GDALClose(filedata);

    band = GDALGetRasterBand(jako->surface_elevation.dataset,1); // 1-based indexing
    cplerr = GDALRasterIO(band,GF_Read,0,0,jako->nx,jako->ny,jako->h,jako->nx,jako->ny,GDT_Float64,0,0);dCPLCHK(cplerr);
    if (jako->verbose) {err = dRealTableView(jako->ny,jako->nx,jako->h,PETSC_VIEWER_STDOUT_WORLD,"mysurf");dCHK(err);}

    filedata = GDALOpen(jako->bed_elevation.path, GA_ReadOnly); if (!filedata) dERROR(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"GDALOpen(\"%s\")",jako->bed_elevation.path);
    if (jako->verbose) {err = JakoFileDataView(filedata,"bed_elevation",PETSC_VIEWER_STDOUT_WORLD);dCHK(err);}
    cplerr = GDALReprojectImage(filedata,NULL,jako->bed_elevation.dataset,NULL,GRA_Bilinear,0.0,0.0,gdalprogress,NULL,NULL);dCPLCHK(cplerr);
    GDALClose(filedata);

    band = GDALGetRasterBand(jako->bed_elevation.dataset,1); // 1-based indexing
    cplerr = GDALRasterIO(band,GF_Read,0,0,jako->nx,jako->ny,jako->b,jako->nx,jako->ny,GDT_Float64,0,0);dCPLCHK(cplerr);
    if (jako->verbose) {err = dRealTableView(jako->ny,jako->nx,jako->b,PETSC_VIEWER_STDOUT_WORLD,"mybed");dCHK(err);}
  }

  dFunctionReturn(0);
}

static dErr VHTCaseSetFromOptions_Jako(VHTCase scase)
{
  VHTCase_Jako *jako = scase->data;
  dBool flg;
  dErr err;
  dFunctionBegin;
  err = PetscOptionsHead("VHTCase_Jako options");dCHK(err); {
    err = PetscOptionsString("-jako_surface_elevation","File to read surface elevation from (assume same projection as model, e.g. UTM)","",
                             jako->surface_elevation.path,jako->surface_elevation.path,PETSC_MAX_PATH_LEN,&flg);dCHK(err);
    if (!flg) dERROR(scase->comm,PETSC_ERR_USER,"User must provide surface elevation file with -jako_surface_elevation FILENAME");

    err = PetscOptionsString("-jako_bed_elevation","File to read bed elevation from (assume same projection as model, e.g. UTM)","",
                                 jako->bed_elevation.path,jako->bed_elevation.path,PETSC_MAX_PATH_LEN,&flg);dCHK(err);
    if (!flg) dERROR(scase->comm,PETSC_ERR_USER,"User must provide bed elevation file with -jako_bed_elevation FILENAME");

    err = PetscOptionsString("-jako_surface_velocity","File to read surface velocity from (assume same projection as model, e.g. UTM)","",
                             jako->surface_velocity.path,jako->surface_velocity.path,PETSC_MAX_PATH_LEN,&flg);dCHK(err);
#if 0 // Not used yet
    if (!flg) dERROR(scase->comm,PETSC_ERR_USER,"User must provide surface velocity file with -jako_surface_velocity FILENAME");
#endif
    err = PetscOptionsBool("-jako_verbose","Turn on verbose output about projections","",jako->verbose,&jako->verbose,NULL);dCHK(err);
    err = PetscOptionsInt("-jako_nx","Number of pixels in x-direction for raster of input fields","",jako->nx,&jako->nx,NULL);dCHK(err);
    err = PetscOptionsInt("-jako_ny","Number of piyels in y-direction for raster of input fields","",jako->ny,&jako->ny,NULL);dCHK(err);
  } err = PetscOptionsTail();dCHK(err);
  err = VHTCaseSetUp_Jako(scase);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseDestroy_Jako(VHTCase scase)
{
  VHTCase_Jako *jako = scase->data;
  dErr err;

  dFunctionBegin;
  OSRDestroySpatialReference(jako->utmref); jako->utmref = NULL;
  OSRDestroySpatialReference(jako->llref);  jako->llref  = NULL;
  OSRDestroySpatialReference(jako->ianref); jako->ianref = NULL;
  OSRDestroySpatialReference(jako->myref);  jako->myref = NULL;

  OCTDestroyCoordinateTransformation(jako->fromutm); jako->fromutm = NULL;
  OCTDestroyCoordinateTransformation(jako->fromll);  jako->fromll = NULL;
  OCTDestroyCoordinateTransformation(jako->fromian); jako->fromian = NULL;

  GDALClose(jako->surface_elevation.dataset);
  GDALClose(jako->bed_elevation.dataset);
  err = dFree(jako->surface_elevation.memory);dCHK(err);
  err = dFree(jako->bed_elevation.memory);dCHK(err);
  err = dFree2(jako->h,jako->b);dCHK(err);
  err = dFree(scase->data);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseForcing_Jako(VHTCase dUNUSED scase,const dReal dUNUSED x[3],dScalar frhou[],dScalar *fp,dScalar *fE)
{
  for (dInt i=0; i<3; i++) frhou[0] = 0;
  fp[0] = 0;
  fE[0] = 0;
  return 0;
}
static dErr VHTCaseCreate_Jako(VHTCase scase)
{
  VHTCase_Jako *jako;
  dErr err;

  dFunctionBegin;
  scase->reality = dTRUE;
  scase->solution = VHTCaseSolution_Jako;
  scase->forcing  = VHTCaseForcing_Jako;
  scase->setfromoptions = VHTCaseSetFromOptions_Jako;
  scase->destroy = VHTCaseDestroy_Jako;

  err = dNew(VHTCase_Jako,&jako);dCHK(err);
  scase->data = jako;

  jako->nx = 8;
  jako->ny = 7;
  dFunctionReturn(0);
}

dErr VHTCaseRegisterAll_Jako(void)
{
  dErr err;

  dFunctionBegin;
  GDALAllRegister();
  err = VHTCaseRegister("wind",VHTCaseCreate_Wind);dCHK(err);
  err = VHTCaseRegister("jako",VHTCaseCreate_Jako);dCHK(err);
  dFunctionReturn(0);
}
