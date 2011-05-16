static const char help[] = "Test nondimensionalization using dUnits\n\n";
#include <dohpunits.h>
#include <dohpsys.h>
#include <dohp.h>

static dErr QuantityView(PetscViewer viewer,const char *name,dScalar value,dUnit unit)
{
  dErr err;

  dFunctionBegin;
  err = PetscViewerASCIIPrintf(viewer,"%-12s: %10.4e nondim    %10.4e %-15s    %10.4e %-15s\n",name,value,dUnitDimensionalize(unit,value),dUnitShortName(unit),dUnitDimensionalizeSI(unit,value),dUnitSIName(unit));dCHK(err);
  dFunctionReturn(0);
}


int main(int argc, char *argv[])
{
  dErr err;
  dUnits units;
  dUnit length,time,mass;
  dUnit forcedensity,force,energy,power,pressure,viscosity,density,accel,velocity,strainrate;
  dScalar H,h3,rho,grav,eta0,eta1,sigma0,sigma1;

  err = dInitialize(&argc,&argv,0,help);dCHK(err);
  err = dUnitsCreate(PETSC_COMM_WORLD,&units);dCHK(err);
  err = dUnitsSetBase(units,dUNITS_LENGTH,"metre","m",1,100,&length);dCHK(err);
  err = dUnitsSetBase(units,dUNITS_TIME,"year","a",31556926,1,&time);dCHK(err);
  err = dUnitsSetBase(units,dUNITS_MASS,"exaton","Et",1e21,1000,&mass);dCHK(err);
  err = dUnitsSetFromOptions(units);dCHK(err);
  // Define a bunch of derived units
  err = dUnitsCreateUnit(units,"DENSITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=-3.0,[dUNITS_MASS]=1.0},&density);dCHK(err);
  err = dUnitsCreateUnit(units,"VELOCITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=1.0,[dUNITS_TIME]=-1.0},&velocity);dCHK(err);
  err = dUnitsCreateUnit(units,"ACCELERATION",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=1.0,[dUNITS_TIME]=-2.0},&accel);dCHK(err);
  err = dUnitsCreateUnit(units,"FORCEDENSITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=-2.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-2.0},&forcedensity);dCHK(err);
  err = dUnitsCreateUnit(units,"FORCE",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=1.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-2.0},&force);dCHK(err);
  err = dUnitsCreateUnit(units,"ENERGY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=2.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-2.0},&energy);dCHK(err);
  err = dUnitsCreateUnit(units,"POWER",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=2.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-1.0},&power);dCHK(err);
  err = dUnitsCreateUnit(units,"PRESSURE",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=-1.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-2.0},&pressure);dCHK(err);
  err = dUnitsCreateUnit(units,"VISCOSITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=-1.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-1.0},&viscosity);dCHK(err);
  err = dUnitsCreateUnit(units,"STRAINRATE",NULL,NULL,4,(dReal[4]){[dUNITS_TIME]=-1.0},&strainrate);dCHK(err);

  err = dUnitsView(units,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
  H = dUnitNonDimensionalizeSI(length,1500);         // typical ice thickness
  h3 = dUnitNonDimensionalizeSI(length,200); h3 = h3*h3*h3; // typical element volume
  rho = dUnitNonDimensionalizeSI(density,910);
  grav = dUnitNonDimensionalizeSI(accel,9.81);
  eta0 = dUnitNonDimensionalizeSI(viscosity,3e15); // Typical large viscosity in Pa s
  eta1 = dUnitNonDimensionalizeSI(viscosity,1e13); // Typical small viscosity in Pa s
  sigma0 = dUnitNonDimensionalizeSI(pressure,1e4); // Typical smallish stress 10 kPa
  sigma1 = dUnitNonDimensionalizeSI(pressure,2e5); // Typical large stress 200 kPa
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"rho",rho,density);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"rhog",rho*grav,forcedensity);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"rhogh3",rho*grav*h3,force);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"eta0",eta0,viscosity);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"eta1",eta1,viscosity);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"sigma0",sigma0,pressure);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"sigma1",sigma1,pressure);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"strainrate0",sigma0/eta0,strainrate);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"strainrate1",sigma1/eta1,strainrate);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"vel0",sigma0/eta0*H,velocity);dCHK(err);
  err = QuantityView(PETSC_VIEWER_STDOUT_WORLD,"vel1",sigma1/eta1*H,velocity);dCHK(err);
  err = dUnitsDestroy(&units);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
