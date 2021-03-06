static const char help[] = "Solve viscous flow coupled to a heat transport problem using dual order elements.\n"
  "The model problem is\n"
  "  -div(eta Du) + grad(p) = f\n"
  "                  div(u) = g\n"
  "    div(u T) - eta Du:Du = h\n"
  "where\n"
  "  D is the symmetric gradient operator\n"
  "  eta(gamma,T) = B(T) (0.5*eps^2 + gamma)^{(p-2)/2}\n"
  "  gamma = Du : Du/2\n"
  "  B(T) = B_0 exp(Q/(n R T))\n"
  "The weak form is\n"
  "  int_Omega eta Dv:Du - p div(v) - q div(u) - v.f - q g -  = 0\n"
  "with Jacobian\n"
  "  int_Omega eta Dv:Du + eta' (Dv:Dw)(Dw:Du) - p div(v) - q div(u) = 0\n"
  "The problem is linear for p=2, an incompressible for g=0\n\n";

#include <dohpstring.h>
#include <dohpviewer.h>
#include <dohpgeom.h>
#include <petscts.h>

#include "vhtimpl.h"

#if 1
#  define VHTAssertRange(val,low,high) do {                             \
    if ((val) < (low)-1e-10 || (high)+1e-10 < (val))                    \
      dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Value %s=%G not in range [%G,%G]",#val,(dReal)(val),(dReal)(low),(dReal)(high)); \
  } while (0)
#else
#  define VHTAssertRange(val,low,high) do {} while (0)
#endif

#if 0
#  define VHTAssertRangeClip(val,low,high) do {                         \
    if ((val) < (low)-1e-10) (val) = (low);                             \
    else if ((val) > (high)+1e-10) (val) = (high);                      \
  } while (0)
#else
#  define VHTAssertRangeClip(val,low,high) VHTAssertRange((val),(low),(high))
#endif

PetscFunctionList VHTCaseList = NULL;

#define VHTCaseType char*

static void VHTStashGetRho(const struct VHTStash *st,const struct VHTRheology *rheo,VHTScalarD *rho);
static void VHTStashGetUH(const struct VHTStash *stash,const dReal dx[9],dReal *uh);
static dErr VHTStashGetStreamlineStabilization(const struct VHTStash *st,const struct VHTRheology *rheo,const dReal dx[9],const dScalar u1[3],dScalar *stab,dScalar *stab1);

dErr VHTCaseRegister(const char *name,VHTCaseCreateFunction screate)
{
  dErr err;
  dFunctionBegin;
  err = PetscFunctionListAdd(PETSC_COMM_WORLD,&VHTCaseList,name,"",(void(*)(void))screate);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseFind(const char *name,VHTCaseCreateFunction *screate)
{
  dErr err;

  dFunctionBegin;
  err = PetscFunctionListFind(PETSC_COMM_WORLD,VHTCaseList,name,PETSC_FALSE,(void(**)(void))screate);dCHK(err);
  if (!*screate) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"VHT Case \"%s\" could not be found",name);
  dFunctionReturn(0);
}
static dErr VHTCaseSetType(VHTCase scase,const VHTCaseType type)
{
  dErr err;
  VHTCaseCreateFunction f;

  dFunctionBegin;
  err = VHTCaseFind(type,&f);dCHK(err);
  err = dStrcpyS(scase->name,sizeof scase->name,type);dCHK(err);
  err = (*f)(scase);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseUnitsSetFromOptions(VHTCase scase)
{
  struct VHTUnitTable *u = &scase->utable;
  dUnits units = scase->units;
  dErr err;
  dInt bprof;
  dBool flg;

  dFunctionBegin;
  err = PetscOptionsInt("-vht_units_profile","Set base units by profile 0=SI, 1=annual, 2=gravitational","",bprof=0,&bprof,NULL);dCHK(err);
  switch (bprof) {
  case 0:
    err = dUnitsGetBase(units,dUNITS_LENGTH,&u->Length);dCHK(err);
    err = dUnitsGetBase(units,dUNITS_TIME,&u->Time);dCHK(err);
    err = dUnitsGetBase(units,dUNITS_MASS,&u->Mass);dCHK(err);
    err = dUnitsGetBase(units,dUNITS_TEMPERATURE,&u->Temperature);dCHK(err);
    break;
  case 1:
    err = dUnitsSetBase(units,dUNITS_LENGTH,"metre","m",1,100,&u->Length);dCHK(err);
    err = dUnitsSetBase(units,dUNITS_TIME,"year","a",31556926,1,&u->Time);dCHK(err);
    err = dUnitsSetBase(units,dUNITS_MASS,"exaton","Et",1e21,1000,&u->Mass);dCHK(err);
    err = dUnitsSetBase(units,dUNITS_TEMPERATURE,"Kelvin","K",1,1,&u->Temperature);dCHK(err);
    break;
  case 2:
    err = dUnitsSetBase(units,dUNITS_LENGTH,"metre","m",1,1,&u->Length);dCHK(err);
    err = dUnitsSetBase(units,dUNITS_TIME,"second","s",1,1e-8,&u->Time);dCHK(err);
    err = dUnitsSetBase(units,dUNITS_MASS,"kilogram","kg",1,1e-9,&u->Mass);dCHK(err);
    err = dUnitsSetBase(units,dUNITS_TEMPERATURE,"Kelvin","K",1,1,&u->Temperature);dCHK(err);
    break;
  case 3:
    err = dUnitsSetBase(units,dUNITS_LENGTH,"metre","m",1,100,&u->Length);dCHK(err);
    err = dUnitsSetBase(units,dUNITS_TIME,"second","s",1,1e-6,&u->Time);dCHK(err);
    err = dUnitsSetBase(units,dUNITS_MASS,"kilogram","kg",1,1e-3,&u->Mass);dCHK(err);
    err = dUnitsSetBase(units,dUNITS_TEMPERATURE,"Kelvin","K",1,1,&u->Temperature);dCHK(err);
    break;
  default: dERROR(scase->comm,PETSC_ERR_ARG_OUTOFRANGE,"Unknown base units profile '%D'",bprof);
  }
  err = dUnitsSetFromOptions(units);dCHK(err);
  err = dUnitsCreateUnit(units,"DENSITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=-3.0,[dUNITS_MASS]=1.0},&u->Density);dCHK(err);
  err = dUnitsCreateUnit(units,"ENERGY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=2.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-2.0},&u->Energy);dCHK(err);
  err = dUnitsCreateUnit(units,"PRESSURE",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=-1.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-2.0},&u->Pressure);dCHK(err);
  u->EnergyDensity = u->Pressure;
  err = dUnitsCreateUnit(units,"MOMENTUM/VOLUME",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=-2.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-1.0},&u->MomentumDensity);dCHK(err);
  err = dUnitsCreateUnit(units,"STRAINRATE",NULL,NULL,4,(dReal[4]){[dUNITS_TIME]=-1.0},&u->StrainRate);dCHK(err);
  err = dUnitsCreateUnit(units,"VELOCITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=1.0,[dUNITS_TIME]=-1.0},&u->Velocity);dCHK(err);
  err = dUnitsCreateUnit(units,"ACCELERATION",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=1.0,[dUNITS_TIME]=-2.0},&u->Acceleration);dCHK(err);
  err = dUnitsCreateUnit(units,"VISCOSITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=-1.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-1.0},&u->Viscosity);dCHK(err);
  err = dUnitsCreateUnit(units,"VOLUME",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=3.0},&u->Volume);dCHK(err);
  err = dUnitsCreateUnit(units,"ENERGY/TEMPERATURE",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=2.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-2.0,[dUNITS_TEMPERATURE]=-1.0},&u->EnergyPerTemperature);dCHK(err);
  err = dUnitsCreateUnit(units,"THERMALCONDUCTIVITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=1.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-3.0,[dUNITS_TEMPERATURE]=-1.0},&u->ThermalConductivity);dCHK(err);
  err = dUnitsCreateUnit(units,"HYDROCONDUCTIVITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=-1.0,[dUNITS_MASS]=1.0,[dUNITS_TIME]=-1.0,[dUNITS_TEMPERATURE]=0},&u->HydroConductivity);dCHK(err);
  err = dUnitsCreateUnit(units,"SPECIFICHEAT",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=2.0,[dUNITS_MASS]=0.0,[dUNITS_TIME]=-2.0,[dUNITS_TEMPERATURE]=-1.0},&u->SpecificHeat);dCHK(err);
  err = dUnitsCreateUnit(units,"DIFFUSIVITY",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=2.0,[dUNITS_TIME]=-1.0},&u->Diffusivity);dCHK(err);
  err = dUnitsCreateUnit(units,"ENERGY/MASS",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=2.0,[dUNITS_MASS]=0.0,[dUNITS_TIME]=-2.0},&u->EnergyPerMass);dCHK(err);
  err = dUnitsCreateUnit(units,"CCGRADIENT",NULL,NULL,4,(dReal[4]){[dUNITS_LENGTH]=1.0,[dUNITS_MASS]=-1.0,[dUNITS_TIME]=2.0,[dUNITS_TEMPERATURE]=1.0},&u->CCGradient);dCHK(err);
  flg = dFALSE;
  err = PetscOptionsGetBool(NULL,"-units_view",&flg,NULL);dCHK(err);
  if (flg) {
    PetscViewer viewer;
    err = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)units),&viewer);dCHK(err);
    err = dUnitsView(units,viewer);dCHK(err);
  }
  dFunctionReturn(0);
}
static dErr VHTCaseProfile_Default(VHTCase scase)
{
  struct VHTRheology *rheo = &scase->rheo;
  dFunctionBegin;
  rheo->B0           = 1;
  rheo->Bomega       = 1;
  rheo->R            = 1;
  rheo->Q            = 1;
  rheo->V            = -0.1;
  rheo->du0          = 1;
  rheo->eps          = 1;
  rheo->pe           = 2;
  rheo->k_T          = 1;
  rheo->kappa_w      = 0.5;
  rheo->c_i          = 1;
  rheo->Latent       = 2;
  rheo->rhoi         = 1;
  rheo->rhow         = 2;
  rheo->beta_CC      = 0.1;
  rheo->T0           = 10;
  rheo->T3           = 10;
  rheo->splice_delta = 1;
  rheo->rscale.momentum = 1;
  rheo->rscale.mass     = 1;
  rheo->rscale.energy   = 1;
  rheo->small         = 1e-10;
  rheo->expstab       = 0;
  rheo->expstab_rate  = 2;
  dFunctionReturn(0);
}
static dErr VHTCaseProfile_Ice(VHTCase scase)
{
  const struct VHTUnitTable *u = &scase->utable;
  struct VHTRheology *rheo = &scase->rheo;
  const dReal
    n = 3,
    refstress_si = 1e5, // 100 kPa
    reftemperature_si = 263.15,
    Asoftness_si = 3.61e-13 * exp(-6e4/(8.314*reftemperature_si)), // Softness parameter
    refstrainrate_si = Asoftness_si * pow(refstress_si,3.0); // about 0.014 / year

  dFunctionBegin;
  // Viscosity at reference strain rate before dimensionless Arrhenius term
  rheo->B0            = dUnitNonDimensionalizeSI(u->Viscosity,pow(Asoftness_si,-1/n) * pow(0.5*dSqr(refstrainrate_si),(1-n)/(2*n)));
  rheo->Bomega        = 181.25; // nondimensional
  rheo->R             = dUnitNonDimensionalizeSI(u->EnergyPerTemperature,8.314); // 8.314 J mol^-1 K^-1
  rheo->Q             = dUnitNonDimensionalizeSI(u->Energy,6.0e4); // 6.0 J mol^-1
  rheo->V             = dUnitNonDimensionalizeSI(u->Volume,-13.0e-6); // m^3 / mol; it always bothered my that the "volume" was negative
  rheo->du0           = dUnitNonDimensionalizeSI(u->StrainRate,refstrainrate_si); // Reference strain rate
  rheo->gamma0        = 0.5*dSqr(rheo->du0); // second invariant of reference strain rate
  rheo->eps           = 1e-3;   // Dimensionless fraction of du0
  rheo->pe            = 1 + 1./n;
  rheo->k_T           = dUnitNonDimensionalizeSI(u->ThermalConductivity,2.1); // W/(m K)
  rheo->kappa_w       = dUnitNonDimensionalizeSI(u->HydroConductivity,1.045e-4); // kg / (m s)
  rheo->c_i           = dUnitNonDimensionalizeSI(u->SpecificHeat,2009);
  rheo->Latent        = dUnitNonDimensionalizeSI(u->EnergyPerMass,3.34e5);
  rheo->rhoi          = dUnitNonDimensionalizeSI(u->Density,910); // Density of ice
  rheo->rhow          = dUnitNonDimensionalizeSI(u->Density,999.8395); // Density of water at 0 degrees C, STP
  rheo->beta_CC       = dUnitNonDimensionalizeSI(u->CCGradient,7.9e-8);
  rheo->T0            = dUnitNonDimensionalizeSI(u->Temperature,260.);
  rheo->T3            = dUnitNonDimensionalizeSI(u->Temperature,273.15);
  rheo->splice_delta  = 1e-3 * rheo->Latent;
  rheo->rscale.momentum = 1;
  rheo->rscale.mass     = 1;
  rheo->rscale.energy   = 1e10;
  rheo->small         = 1e-10;
  rheo->expstab       = 0;
  rheo->expstab_rate  = 2;
  dFunctionReturn(0);
}
static dErr VHTCaseSetFromOptions(VHTCase scase)
{
  struct VHTRheology *rheo = &scase->rheo;
  struct VHTUnitTable *u = &scase->utable;
  PetscFunctionList profiles = NULL;
  char prof[256] = "default";
  dErr (*rprof)(VHTCase);
  dErr err;

  dFunctionBegin;
  err = PetscFunctionListAdd(PETSC_COMM_WORLD,&profiles,"default",NULL,(void(*)(void))VHTCaseProfile_Default);dCHK(err);
  err = PetscFunctionListAdd(PETSC_COMM_WORLD,&profiles,"ice",NULL,(void(*)(void))VHTCaseProfile_Ice);dCHK(err);
  err = PetscOptionsBegin(scase->comm,NULL,"VHTCase options",__FILE__);dCHK(err); {
    err = PetscOptionsList("-rheo_profile","Rheological profile","",profiles,prof,prof,sizeof prof,NULL);dCHK(err);
    err = PetscFunctionListFind(scase->comm,profiles,prof,PETSC_FALSE,(void(**)(void))&rprof);dCHK(err);
    err = (*rprof)(scase);dCHK(err);
    rheo->mask_kinetic  = 0;
    rheo->mask_momtrans = 1;
    rheo->mask_rho      = 0;
    rheo->mask_Ep       = 1;
    rheo->eta_min       = 1e-10 * rheo->B0;
    rheo->eta_max       = 1e5 * rheo->B0;
    err = PetscFunctionListDestroy(&profiles);dCHK(err);
    err = dOptionsRealUnits("-rheo_B0","Viscosity at reference strain rate and temperature","",u->Viscosity,rheo->B0,&rheo->B0,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_Bomega","Softening due to water content","",NULL,rheo->Bomega,&rheo->Bomega,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_R","Ideal gas constant","",u->EnergyPerTemperature,rheo->R,&rheo->R,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_Q","Activation Energy","",u->Energy,rheo->Q,&rheo->Q,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_V","Activation Volume","",u->Volume,rheo->V,&rheo->V,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_du0","Reference strain rate","",u->StrainRate,rheo->du0,&rheo->du0,NULL);dCHK(err);
    rheo->gamma0 = 0.5*dSqr(rheo->du0);
    err = dOptionsRealUnits("-rheo_eps","Strain rate regularization (fraction of reference du0)","",NULL,rheo->eps,&rheo->eps,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_p","Power p=1+1/n where n is Glen exponent","",NULL,rheo->pe,&rheo->pe,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_k_T","Thermal conductivity in the cold part","",u->ThermalConductivity,rheo->k_T,&rheo->k_T,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_kappa_w","Hydraulic conductivity in the warm part","",u->HydroConductivity,rheo->kappa_w,&rheo->kappa_w,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_c_i","Specific heat capacity of cold part","",u->SpecificHeat,rheo->c_i,&rheo->c_i,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_Latent","Latent heat of fusion","",u->EnergyPerMass,rheo->Latent,&rheo->Latent,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_rhoi","Density of cold part","",u->Density,rheo->rhoi,&rheo->rhoi,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_rhow","Density of melted part","",u->Density,rheo->rhow,&rheo->rhow,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_beta_CC","Clausius-Clapeyron gradient","",u->CCGradient,rheo->beta_CC,&rheo->beta_CC,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_T0","Reference temperature (corresponds to solution Energy=0)","",u->Temperature,rheo->T0,&rheo->T0,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_p0","Reference pressure (corresponds to solution pressure=0)","",u->Pressure,rheo->p0,&rheo->p0,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_T3","Triple point temperature","",u->Temperature,rheo->T3,&rheo->T3,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_splice_delta","Characteristic width of split","",u->EnergyPerMass,rheo->splice_delta,&rheo->splice_delta,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_gravity_x","Gravity in x-direction","",u->Acceleration,rheo->gravity[0],&rheo->gravity[0],NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_gravity_y","Gravity in y-direction","",u->Acceleration,rheo->gravity[1],&rheo->gravity[1],NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_gravity_z","Gravity in z-direction","",u->Acceleration,rheo->gravity[2],&rheo->gravity[2],NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_Kstab","Stabilization diffusivity, use to prevent negative diffusivity due to non-monotone splice","",NULL,rheo->Kstab,&rheo->Kstab,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_supg","Multiplier for SU/PG stabilization","",NULL,rheo->supg,&rheo->supg,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_supg_crosswind","Fraction of streamline diffusion to apply in crosswind direction","",NULL,rheo->supg_crosswind,&rheo->supg_crosswind,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_shockhmm","Multiplier for the Hughes-Mallet-Mizukami shock-capturing stabilization","",NULL,rheo->shockhmm,&rheo->shockhmm,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_expstab","Multiplier for the exponential low-temperature stabilization","",NULL,rheo->expstab,&rheo->expstab,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_expstab_rate","Inverse characteristic width of exponential low-temperature stabilization","",NULL,rheo->expstab_rate,&rheo->expstab_rate,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_mask_kinetic","Coefficient of kinetic energy used to determine internal energy","",NULL,rheo->mask_kinetic,&rheo->mask_kinetic,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_mask_momtrans","Multiplier for the transport term in momentum balance","",NULL,rheo->mask_momtrans,&rheo->mask_momtrans,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_mask_rho","Multiplier for density variation due to omega","",NULL,rheo->mask_rho,&rheo->mask_rho,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_mask_Ep","Multiplier for pressure in (E+p) term of energy equation","",NULL,rheo->mask_Ep,&rheo->mask_Ep,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_eta_min","Lower bound for effective viscosity","",u->Viscosity,rheo->eta_min,&rheo->eta_min,NULL);dCHK(err);
    err = dOptionsRealUnits("-rheo_eta_max","Upper bound for effective viscosity","",u->Viscosity,rheo->eta_max,&rheo->eta_max,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_rscale_momentum","Scaling for momentum residual","",rheo->rscale.momentum,&rheo->rscale.momentum,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_rscale_mass","Scaling for mass residual","",rheo->rscale.mass,&rheo->rscale.mass,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_rscale_energy","Scaling for energy residual","",rheo->rscale.energy,&rheo->rscale.energy,NULL);dCHK(err);
    err = PetscOptionsReal("-rheo_small","Small parameter to avoid dividing by zero","",rheo->small,&rheo->small,NULL);dCHK(err);
    if (scase->setfromoptions) {err = (*scase->setfromoptions)(scase);dCHK(err);}
  } err = PetscOptionsEnd();dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseDestroy(VHTCase *scase)
{
  dErr err;

  dFunctionBegin;
  if ((*scase)->destroy) {err = ((*scase)->destroy)(*scase);dCHK(err);}
  err = dUnitsDestroy(&(*scase)->units);dCHK(err);
  err = dFree(*scase);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTCaseRegisterAll(void)
{
  dErr err;

  dFunctionBegin;
  err = VHTCaseRegisterAll_Exact();dCHK(err);
#if defined dHAVE_GDAL
  err = VHTCaseRegisterAll_Jako();dCHK(err);
#endif
  dFunctionReturn(0);
}
static dUNUSED dErr VHTLogEpochViewLine2_Private(PetscViewer viewer,const char *linename,const char *n1,const dReal v1[2],const char *n2,const dReal v2[2])
{return PetscViewerASCIIPrintf(viewer,"%s: %-9s [% 8.2e,% 8.2e]  %-9s [% 8.2e,% 8.2e]\n",linename,n1,v1[0],v1[1],n2,v2[0],v2[1]);}
static dUNUSED dErr VHTLogEpochViewLine3_Private(PetscViewer viewer,const char *linename,const char *n1,const dReal v1[2],const char *n2,const dReal v2[2],const char *n3,const dReal v3[2])
{return PetscViewerASCIIPrintf(viewer,"%s: %-9s [% 8.2e,% 8.2e]  %-9s [% 8.2e,% 8.2e]  %-9s [% 8.2e,% 8.2e]\n",linename,n1,v1[0],v1[1],n2,v2[0],v2[1],n3,v3[0],v3[1]);}
static dErr VHTLogEpochViewLine4_Private(PetscViewer viewer,const char *linename,const char *n1,const dReal v1[2],const char *n2,const dReal v2[2],const char *n3,const dReal v3[2],const char *n4,const dReal v4[2])
{return PetscViewerASCIIPrintf(viewer,"%s: %-9s [% 8.2e,% 8.2e]  %-9s [% 8.2e,% 8.2e]  %-9s [% 8.2e,% 8.2e]  %-9s [% 8.2e,% 8.2e]\n",linename,n1,v1[0],v1[1],n2,v2[0],v2[1],n3,v3[0],v3[1],n4,v4[0],v4[1]);}
#define VHTLogEpochViewLine2(viewer,linename,ep,f1,f2) VHTLogEpochViewLine2_Private((viewer),linename,#f1,(ep)->f1,#f2,(ep)->f2)
#define VHTLogEpochViewLine3(viewer,linename,ep,f1,f2,f3) VHTLogEpochViewLine3_Private((viewer),linename,#f1,(ep)->f1,#f2,(ep)->f2,#f3,(ep)->f3)
#define VHTLogEpochViewLine4(viewer,linename,ep,f1,f2,f3,f4) VHTLogEpochViewLine4_Private((viewer),linename,#f1,(ep)->f1,#f2,(ep)->f2,#f3,(ep)->f3,#f4,(ep)->f4)
static dErr VHTLogEpochView(struct VHTLogEpoch *ep,PetscViewer viewer,const char *fmt,...)
{
  va_list Argp;
  char name[4096];
  size_t fullLen;
  dErr err;

  dFunctionBegin;
  va_start(Argp,fmt);
  err = PetscVSNPrintf(name,sizeof name,fmt,&fullLen,Argp);dCHK(err);
  va_end(Argp);
  err = VHTLogEpochViewLine4(viewer,name,ep,p,E,T,omega);dCHK(err);
  err = VHTLogEpochViewLine4(viewer,name,ep,eta,cPeclet,cReynolds,Prandtl);dCHK(err);
  err = VHTLogEpochViewLine3(viewer,name,ep,speed,upwind,K1);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTLogView(struct VHTLog *vlog,PetscViewer viewer)
{
  dErr err;

  dFunctionBegin;
  err = PetscViewerASCIIPrintf(viewer,"Logged %d epochs\n",vlog->epoch+1);dCHK(err);
  err = VHTLogEpochView(&vlog->global,viewer,"Global");dCHK(err);
  err = VHTLogEpochView(&vlog->epochs[vlog->epoch],viewer,"Latest");dCHK(err);
  dFunctionReturn(0);
}
static void VHTLogEpochRangeReset(dReal a[2]) {a[0] = PETSC_MAX_REAL; a[1] = PETSC_MIN_REAL;}
static void VHTLogEpochRangeUpdate(dReal a[2],dReal b) {a[0] = dMin(a[0],b); a[1] = dMax(a[1],b);}
static void VHTLogEpochRangeUpdate2(dReal a[2],const dReal b[2]) {a[0] = dMin(a[0],b[0]); a[1] = dMax(a[1],b[1]);}
static dErr VHTLogEpochReset(struct VHTLogEpoch *ep)
{
  dReal (*pairs)[2] = (dReal(*)[2])ep;
  dInt n = (dInt)(sizeof(*ep)/sizeof(pairs[0]));
  dFunctionBegin;
  for (dInt i=0; i<n; i++) VHTLogEpochRangeReset(pairs[i]);
  dFunctionReturn(0);
}
static dErr VHTLogEpochStart(struct VHTLog *vlog)
{
  dErr err;

  dFunctionBegin;
  vlog->epoch++;
  if (vlog->epoch >= vlog->alloc) {
    dInt newalloc = vlog->alloc * 2 + 16;
    struct VHTLogEpoch *tmp = vlog->epochs;
    err = dCallocA(newalloc,&vlog->epochs);dCHK(err);
    err = dMemcpy(vlog->epochs,tmp,vlog->alloc*sizeof(tmp[0]));dCHK(err);
    err = dFree(tmp);dCHK(err);
    vlog->alloc = newalloc;
  }
  err = VHTLogEpochReset(&vlog->epochs[vlog->epoch]);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTLogEpochEnd(struct VHTLog *vlog)
{
  struct VHTLogEpoch *g = &vlog->global,*e = &vlog->epochs[vlog->epoch];
  dReal (*gpairs)[2] = (dReal(*)[2])g,
    (*epairs)[2] = (dReal(*)[2])e;
  dInt n = (dInt)(sizeof(*g)/sizeof(gpairs[0]));
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<n; i++) VHTLogEpochRangeUpdate2(gpairs[i],epairs[i]);
  if (vlog->monitor) {err = VHTLogEpochView(e,PETSC_VIEWER_STDOUT_WORLD,"Epoch[% 3d]",vlog->epoch);dCHK(err);}
  dFunctionReturn(0);
}
static dErr VHTLogStash(struct VHTLog *vlog,struct VHTRheology *rheo,const dReal dx[9],const struct VHTStash *stash)
{
  struct VHTLogEpoch *ep = &vlog->epochs[vlog->epoch];
  dScalar u1[3] = {0,0,0};
  dReal speed,cPeclet,cReynolds,uh,Prandtl,upwind,upwind1;
  VHTScalarD rho;
  dErr err;

  dFunctionBegin;
  VHTStashGetRho(stash,rheo,&rho);
  VHTStashGetUH(stash,dx,&uh);
  speed = dSqrt(dDotScalar3(stash->u,stash->u));
  err = VHTStashGetStreamlineStabilization(stash,rheo,dx,u1,&upwind,&upwind1);dCHK(err);
  upwind *= dDotScalar3(stash->u,stash->u);
  cPeclet = uh / stash->K[1];
  cReynolds = rho.x * uh / stash->eta.x;
  Prandtl = rheo->c_i * stash->eta.x / rheo->k_T;
  VHTLogEpochRangeUpdate(ep->eta,stash->eta.x);
  VHTLogEpochRangeUpdate(ep->cPeclet,cPeclet);
  VHTLogEpochRangeUpdate(ep->cReynolds,cReynolds);
  VHTLogEpochRangeUpdate(ep->p,stash->p);
  VHTLogEpochRangeUpdate(ep->E,stash->E);
  VHTLogEpochRangeUpdate(ep->T,stash->T.x);
  VHTLogEpochRangeUpdate(ep->omega,stash->omega.x);
  VHTLogEpochRangeUpdate(ep->Prandtl,Prandtl);
  VHTLogEpochRangeUpdate(ep->K1,stash->K[1]);
  VHTLogEpochRangeUpdate(ep->upwind,upwind);
  VHTLogEpochRangeUpdate(ep->speed,speed);
  dFunctionReturn(0);
}
static dErr VHTLogSetFromOptions(struct VHTLog *vlog)
{
  dErr err;

  dFunctionBegin;
  err = PetscOptionsBool("-vht_log_monitor","View each epoch",NULL,vlog->monitor,&vlog->monitor,NULL);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTLogReset(struct VHTLog *vlog)
{
  dFunctionBegin;
  vlog->epoch = 0;
  dFunctionReturn(0);
}

#define ALEN(a) ((dInt)(sizeof(a)/sizeof(a)[0]))
static PetscLogEvent LOG_VHTShellMult;

static dErr VHTGetNullSpace(VHT vht,MatNullSpace *matnull);
static dErr MatMultXIorA_VHT_stokes(Mat A,Vec gx,Vec gy,Vec gz,InsertMode,VHTMultMode);
static dErr MatMult_Nest_VHT_all(Mat J,Vec gx,Vec gy);
static dErr MatMult_VHT_uu(Mat A,Vec gx,Vec gy) {return MatMultXIorA_VHT_stokes(A,gx,gy,NULL,INSERT_VALUES,VHT_MULT_UU);}
static dErr MatMult_VHT_up(Mat A,Vec gx,Vec gy) {return MatMultXIorA_VHT_stokes(A,gx,gy,NULL,INSERT_VALUES,VHT_MULT_UP);}
static dErr MatMult_VHT_ue(Mat A,Vec gx,Vec gy) {return MatMultXIorA_VHT_stokes(A,gx,gy,NULL,INSERT_VALUES,VHT_MULT_UE);}
static dErr MatMult_VHT_pu(Mat A,Vec gx,Vec gy) {return MatMultXIorA_VHT_stokes(A,gx,gy,NULL,INSERT_VALUES,VHT_MULT_PU);}
static dErr MatMultAdd_VHT_uu(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultXIorA_VHT_stokes(A,gx,gy,gz,ADD_VALUES,VHT_MULT_UU);}
static dErr MatMultAdd_VHT_up(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultXIorA_VHT_stokes(A,gx,gy,gz,ADD_VALUES,VHT_MULT_UP);}
static dErr MatMultAdd_VHT_ue(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultXIorA_VHT_stokes(A,gx,gy,gz,ADD_VALUES,VHT_MULT_UE);}
static dErr MatMultAdd_VHT_pu(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultXIorA_VHT_stokes(A,gx,gy,gz,ADD_VALUES,VHT_MULT_PU);}
static dErr MatMult_VHT_eX(Mat A,Vec gx,Vec gy,VHTMultMode);
static dErr MatMult_VHT_eu(Mat A,Vec gx,Vec gy) {return MatMult_VHT_eX(A,gx,gy,VHT_MULT_EU);}
static dErr MatMult_VHT_ep(Mat A,Vec gx,Vec gy) {return MatMult_VHT_eX(A,gx,gy,VHT_MULT_EP);}
static dErr MatMult_VHT_ee(Mat A,Vec gx,Vec gy) {return MatMult_VHT_eX(A,gx,gy,VHT_MULT_EE);}
static dErr MatMultAdd_VHT_eX(Mat A,Vec gx,Vec gy,Vec gz,VHTMultMode);
static dErr MatMultAdd_VHT_eu(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultAdd_VHT_eX(A,gx,gy,gz,VHT_MULT_EU);}
static dErr MatMultAdd_VHT_ep(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultAdd_VHT_eX(A,gx,gy,gz,VHT_MULT_EP);}
static dErr MatMultAdd_VHT_ee(Mat A,Vec gx,Vec gy,Vec gz) {return MatMultAdd_VHT_eX(A,gx,gy,gz,VHT_MULT_EE);}
static dErr MatGetVecs_VHT_stokes(Mat,Vec*,Vec*);
static dErr MatGetVecs_VHT_ee(Mat,Vec*,Vec*);

static dErr VHTCreate(MPI_Comm comm,VHT *invht)
{
  VHT vht;
  dErr err;

  dFunctionBegin;
  *invht = 0;
  err = dNew(struct _n_VHT,&vht);dCHK(err);
  vht->comm = comm;

  vht->velocityBDeg  = 3;
  vht->pressureCodim = 1;
  vht->energyBDeg  = 3;
  vht->u_dirichlet[0]  = 100;
  vht->u_dirichlet[1]  = 200;
  vht->u_dirichlet[2]  = 300;
  vht->e_dirichlet[0]  = 100;
  vht->e_dirichlet[1]  = 200;
  vht->e_dirichlet[2]  = 300;
  vht->alldirichlet  = dTRUE;
  vht->domain_error  = 2;
  vht->function_qmethod = dQUADRATURE_METHOD_FAST;
  vht->jacobian_qmethod = dQUADRATURE_METHOD_SPARSE;
  vht->clip.p[0] = SNES_VI_NINF;
  vht->clip.p[1] = SNES_VI_INF;
  vht->clip.E[0] = SNES_VI_NINF;
  vht->clip.E[1] = SNES_VI_INF;

  err = dCalloc(sizeof(*vht->scase),&vht->scase);dCHK(err);
  vht->scase->comm = comm;

  vht->log.epoch = -1;
  err = VHTLogEpochReset(&vht->log.global);dCHK(err);

  err = dUnitsCreate(vht->comm,&vht->scase->units);dCHK(err);

  *invht = vht;
  dFunctionReturn(0);
}

static dErr MatGetVecs_VHT_stokes(Mat A,Vec *x,Vec *y)
{
  VHT vht;
  dInt m,n,nu,np;
  dErr err;

  dFunctionBegin;
  err = MatShellGetContext(A,(void**)&vht);dCHK(err);
  err = MatGetLocalSize(A,&m,&n);dCHK(err);
  err = VecGetLocalSize(vht->gvelocity,&nu);dCHK(err);
  err = VecGetLocalSize(vht->gpressure,&np);dCHK(err);
  if (nu==np) dERROR(PETSC_COMM_SELF,1,"Degenerate case, don't know which space to copy");
  if (x) {
    if (n == nu) {
      err = VecDuplicate(vht->gvelocity,x);dCHK(err);
    } else if (n == np) {
      err = VecDuplicate(vht->gpressure,x);dCHK(err);
    } else dERROR(PETSC_COMM_SELF,1,"sizes do not agree with either space");
  }
  if (y) {
    if (n == nu) {
      err = VecDuplicate(vht->gvelocity,y);dCHK(err);
    } else if (n == np) {
      err = VecDuplicate(vht->gpressure,y);dCHK(err);
    } else dERROR(PETSC_COMM_SELF,1,"sizes do not agree with either space");
  }
  dFunctionReturn(0);
}

static dErr MatGetVecs_VHT_ee(Mat A,Vec *x,Vec *y)
{
  VHT vht;
  dErr err;

  dFunctionBegin;
  err = MatShellGetContext(A,(void**)&vht);dCHK(err);
  if (x) {err = VecDuplicate(vht->genergy,x);dCHK(err);}
  if (y) {err = VecDuplicate(vht->genergy,y);dCHK(err);}
  dFunctionReturn(0);
}
static dErr dIntArrayJoin(char buf[],size_t bufsize,const char join[],const dInt values[],dInt n)
{ // Concatenates list of values until a zero value is found, up to a maximum of n
  char *p = buf;

  dFunctionBegin;
  buf[0] = 0;
  for (dInt i=0; i<n && values[i]; i++) {
    p += snprintf(p,bufsize-(p-buf),"%s%d",i ? join : "",values[i]);
  }
  dFunctionReturn(0);
}
static dErr VHTView(VHT vht,PetscViewer viewer)
{
  dErr err;
  dInt nu,np,ne;

  dFunctionBegin;
  err = PetscViewerASCIIGetStdout(vht->comm,&viewer);dCHK(err);
  err = VecGetSize(vht->gvelocity,&nu);dCHK(err);
  err = VecGetSize(vht->gpressure,&np);dCHK(err);
  err = VecGetSize(vht->genergy,&ne);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"VHT: basis degree: u %D  p %D  E %D\n",vht->velocityBDeg,vht->velocityBDeg-vht->pressureCodim,vht->energyBDeg);dCHK(err);
  err = PetscViewerASCIIPushTab(viewer);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Number of nodes: u %D  p %D  E %D\n",nu/3,np,ne);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"Quadrature methods: F:%s  B:%s\n",dQuadratureMethods[vht->function_qmethod],dQuadratureMethods[vht->jacobian_qmethod]);dCHK(err);
  {
    char buf[256];
    err = dIntArrayJoin(buf,sizeof buf,", ",vht->u_dirichlet,ALEN(vht->u_dirichlet));dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Dirichlet velocity at: %s\n",buf);dCHK(err);
    err = dIntArrayJoin(buf,sizeof buf,", ",vht->e_dirichlet,ALEN(vht->e_dirichlet));dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Dirichlet energy at  : %s\n",buf);dCHK(err);
  }
  {
    struct VHTRheology *rheo = &vht->scase->rheo;
    struct VHTUnitTable *u = &vht->scase->utable;
    err = PetscViewerASCIIPrintf(viewer,"Reference strain rate: %G %s\n",dUnitDimensionalize(u->StrainRate,rheo->du0),dUnitName(u->StrainRate));dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Reference temperature: %G %s\n",dUnitDimensionalize(u->Temperature,rheo->T0),dUnitName(u->Temperature));dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Reference viscosity  : %G %s\n",dUnitDimensionalize(u->Viscosity,rheo->B0),  dUnitName(u->Viscosity));dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Momentum diffusivity : %G %s\n",dUnitDimensionalize(u->Diffusivity,rheo->B0/rheo->rhoi),  dUnitName(u->Diffusivity));dCHK(err); // kinematic viscosity
    err = PetscViewerASCIIPrintf(viewer,"Thermal diffusivity  : %G %s\n",dUnitDimensionalize(u->Diffusivity,rheo->k_T/(rheo->rhoi*rheo->c_i)),dUnitName(u->Diffusivity));dCHK(err);
    err = PetscViewerASCIIPrintf(viewer,"Hydraulic diffusivity: %G %s\n",dUnitDimensionalize(u->Diffusivity,rheo->kappa_w/rheo->rhoi),dUnitName(u->Diffusivity));dCHK(err);
    if (rheo->Kstab > 0) {err = PetscViewerASCIIPrintf(viewer,"Artificial diffusivity: %G %s\n",dUnitDimensionalize(u->Diffusivity,rheo->Kstab),dUnitName(u->Diffusivity));dCHK(err);}
  }
  err = PetscViewerASCIIPrintf(viewer,"Clipping: p [%G,%G]  E [%G,%G]\n",vht->clip.p[0],vht->clip.p[1],vht->clip.E[0],vht->clip.E[1]);dCHK(err);
  err = dRealTableView(3,2,(const dReal*)vht->scase->bbox,PETSC_VIEWER_STDOUT_WORLD,"Mesh bounding box");dCHK(err);

  err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  dFunctionReturn(0);
}
static void VHTMeshMorph_NonDimensionalizeSI(void *vctx,double *coords)
{
  VHT vht = vctx;
  dReal scale = dUnitNonDimensionalizeSI(vht->scase->utable.Length,1.0); // nondimensionalize 1 meter
  double x = coords[0],y = coords[1],z = coords[2];
  coords[0] = scale * x;
  coords[1] = scale * y;
  coords[2] = scale * z;
}
static dErr VHTSetFromOptions(VHT vht)
{
  char scasename[256] = "Exact0";
  dMesh mesh;
  dFS fsu,fsp,fse;
  dJacobi jac;
  dMeshESH domain;
  dMeshTag dutag,dptag,detag;
  dBool doview = dFALSE;
  dErr err;

  dFunctionBegin;
  err = VHTCaseUnitsSetFromOptions(vht->scase);dCHK(err);
  err = dStrcpyS(vht->mattype_Buu,sizeof(vht->mattype_Buu),MATBAIJ);dCHK(err);
  err = dStrcpyS(vht->mattype_Bpp,sizeof(vht->mattype_Bpp),MATAIJ);dCHK(err);
  err = dStrcpyS(vht->mattype_Bee,sizeof(vht->mattype_Bee),MATAIJ);dCHK(err);
  err = PetscOptionsBegin(vht->comm,NULL,"Viscous Heat Transport options",__FILE__);dCHK(err); {
    err = PetscOptionsInt("-vht_u_bdeg","Constant isotropic degree to use for velocity","",vht->velocityBDeg,&vht->velocityBDeg,NULL);dCHK(err);
    err = PetscOptionsInt("-vht_p_codim","Reduce pressure space by this factor","",vht->pressureCodim,&vht->pressureCodim,NULL);dCHK(err);
    err = PetscOptionsInt("-vht_e_bdeg","Constant isotropic degree to use for energy","",vht->energyBDeg,&vht->energyBDeg,NULL);dCHK(err);
    err = PetscOptionsBool("-vht_cardinal_mass","Assemble diagonal mass matrix","",vht->cardinalMass,&vht->cardinalMass,NULL);dCHK(err);
    err = PetscOptionsList("-vht_Buu_mat_type","Matrix type for velocity-velocity operator","",MatList,vht->mattype_Buu,vht->mattype_Buu,sizeof(vht->mattype_Buu),NULL);dCHK(err);
    err = PetscOptionsList("-vht_Bpp_mat_type","Matrix type for pressure-pressure operator","",MatList,vht->mattype_Bpp,vht->mattype_Bpp,sizeof(vht->mattype_Bpp),NULL);dCHK(err);
    err = PetscOptionsList("-vht_Bee_mat_type","Matrix type for energy-energy operator","",MatList,vht->mattype_Bee,vht->mattype_Bee,sizeof(vht->mattype_Bee),NULL);dCHK(err);
    err = PetscOptionsBool("-vht_split_recursive","Recursively nest the Stokes block inside the overall matrix instead of using a flat split","",vht->split_recursive,&vht->split_recursive,NULL);dCHK(err);
    err = PetscOptionsEnum("-vht_f_qmethod","Quadrature method for residual evaluation/matrix-free","",dQuadratureMethods,(PetscEnum)vht->function_qmethod,(PetscEnum*)&vht->function_qmethod,NULL);dCHK(err);
    err = PetscOptionsEnum("-vht_jac_qmethod","Quadrature to use for Jacobian assembly","",dQuadratureMethods,(PetscEnum)vht->jacobian_qmethod,(PetscEnum*)&vht->jacobian_qmethod,NULL);dCHK(err);
    {
      dBool flg; dInt n = ALEN(vht->u_dirichlet);
      err = PetscOptionsIntArray("-vht_u_dirichlet","List of boundary sets on which to impose Dirichlet velocity conditions","",vht->u_dirichlet,&n,&flg);dCHK(err);
      if (flg) {
        for (dInt i=n; i<ALEN(vht->u_dirichlet); i++) vht->u_dirichlet[i] = 0; /* Clear out any leftover values */
        if (n < 3) vht->alldirichlet = dFALSE;                             /* @bug More work to determine independent of the mesh whether all the boundaries are Dirichlet */
      }
      n = ALEN(vht->e_dirichlet);
      err = PetscOptionsIntArray("-vht_e_dirichlet","List of boundary sets on which to impose Dirichlet energy conditions","",vht->e_dirichlet,&n,&flg);dCHK(err);
      if (flg) for (dInt i=n; i<ALEN(vht->e_dirichlet); i++) vht->e_dirichlet[i] = 0;

    }
    err = PetscOptionsBool("-vht_alldirichlet","Remove the pressure null space (all Dirichlet momentum boundary conditions)","",vht->alldirichlet,&vht->alldirichlet,NULL);dCHK(err);
    err = PetscOptionsInt("-vht_domain_error","Throw an error on domain error instead of just calling SNESFunctionDomainError()","",vht->domain_error,&vht->domain_error,NULL);dCHK(err);
    {
      dInt two = 2;
      err = PetscOptionsRealArray("-vht_clip_p","Clip pressure to this range","",vht->clip.p,&two,NULL);dCHK(err);
      two = 2;
      err = PetscOptionsRealArray("-vht_clip_E","Clip Energy density to this range","",vht->clip.E,&two,NULL);dCHK(err);
    }
    err = PetscOptionsList("-vht_case","Which sort of case to run","",VHTCaseList,scasename,scasename,sizeof(scasename),NULL);dCHK(err);
    err = PetscOptionsBool("-vht_view","View VHT configuration","",doview,&doview,NULL);dCHK(err);
    err = VHTLogSetFromOptions(&vht->log);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);

  err = dMeshCreate(vht->comm,&mesh);dCHK(err);
  err = dMeshSetInFile(mesh,"dblock.h5m",NULL);dCHK(err);
  err = dMeshSetFromOptions(mesh);dCHK(err);
  err = dMeshLoad(mesh);dCHK(err);
  err = dMeshMorph(mesh,VHTMeshMorph_NonDimensionalizeSI,vht);dCHK(err);
  err = dMeshGetRoot(mesh,&domain);dCHK(err); /* Need a taggable set */
  err = dMeshSetDuplicateEntsOnly(mesh,domain,&domain);dCHK(err);
  err = PetscObjectSetName((PetscObject)mesh,"dMesh_0");dCHK(err);

  err = dJacobiCreate(vht->comm,&jac);dCHK(err);
  err = dJacobiSetFromOptions(jac);dCHK(err);

  err = dMeshCreateRuleTagIsotropic(mesh,domain,"vht_efs_velocity_degree",vht->velocityBDeg,&dutag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,"vht_efs_pressure_degree",vht->velocityBDeg-vht->pressureCodim,&dptag);dCHK(err);
  err = dMeshCreateRuleTagIsotropic(mesh,domain,"vht_efs_energy_degree",vht->energyBDeg,&detag);dCHK(err);

  err = dFSCreate(vht->comm,&fsu);dCHK(err);
  err = dFSSetBlockSize(fsu,3);dCHK(err);
  err = dFSSetMesh(fsu,mesh,domain);dCHK(err);
  err = dFSSetDegree(fsu,jac,dutag);dCHK(err);
  for (dInt i=0; i<ALEN(vht->u_dirichlet) && vht->u_dirichlet[i]>0; i++) {
    err = dFSRegisterBoundary(fsu,vht->u_dirichlet[i],dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  }
  err = PetscObjectSetOptionsPrefix((dObject)fsu,"u");dCHK(err);
  err = dFSSetFromOptions(fsu);dCHK(err);
  err = PetscObjectSetName((PetscObject)fsu,"dFS_U_0");dCHK(err);
  err = dFSSetFieldName(fsu,0,"MomX");dCHK(err);
  err = dFSSetFieldName(fsu,1,"MomY");dCHK(err);
  err = dFSSetFieldName(fsu,2,"MomZ");dCHK(err);
  err = dFSSetFieldUnit(fsu,0,vht->scase->utable.Pressure);dCHK(err); // Energy/Volume is the same as Pressure
  err = dFSSetFieldUnit(fsu,1,vht->scase->utable.Pressure);dCHK(err); // Energy/Volume is the same as Pressure
  err = dFSSetFieldUnit(fsu,2,vht->scase->utable.Pressure);dCHK(err); // Energy/Volume is the same as Pressure
  vht->fsu = fsu;

  err = dFSCreate(vht->comm,&fsp);dCHK(err);
  err = dFSSetMesh(fsp,mesh,domain);dCHK(err);
  err = dFSSetDegree(fsp,jac,dptag);dCHK(err);
  err = PetscObjectSetOptionsPrefix((dObject)fsp,"p");dCHK(err);
  /* No boundaries, the pressure space has Neumann conditions when Dirichlet velocity conditions are applied */
  err = dFSSetFromOptions(fsp);dCHK(err);
  err = PetscObjectSetName((PetscObject)fsp,"dFS_P_0");dCHK(err);
  err = dFSSetFieldName(fsp,0,"Pressure");dCHK(err);
  err = dFSSetFieldUnit(fsp,0,vht->scase->utable.Pressure);dCHK(err);
  vht->fsp = fsp;

  err = dFSCreate(vht->comm,&fse);dCHK(err);
  err = dFSSetMesh(fse,mesh,domain);dCHK(err);
  err = dFSSetDegree(fse,jac,detag);dCHK(err);
  for (dInt i=0; i<ALEN(vht->e_dirichlet) && vht->e_dirichlet[i]>0; i++) {
    err = dFSRegisterBoundary(fse,vht->e_dirichlet[i],dFSBSTATUS_DIRICHLET,NULL,NULL);dCHK(err);
  }
  err = PetscObjectSetOptionsPrefix((dObject)fse,"e");dCHK(err);
  err = dFSSetFromOptions(fse);dCHK(err);
  err = PetscObjectSetName((PetscObject)fse,"dFS_E_0");dCHK(err);
  err = dFSSetFieldName(fse,0,"EnergyDensity");dCHK(err);
  err = dFSSetFieldUnit(fse,0,vht->scase->utable.Pressure);dCHK(err); // Energy/Volume is the same as Pressure
  vht->fse = fse;

  err = dFSCreateExpandedVector(fsu,&vht->xu);dCHK(err);
  err = VecDuplicate(vht->xu,&vht->yu);dCHK(err);

  err = dFSCreateExpandedVector(fsp,&vht->xp);dCHK(err);
  err = VecDuplicate(vht->xp,&vht->yp);dCHK(err);

  err = dFSCreateExpandedVector(fse,&vht->xe);dCHK(err);
  err = VecDuplicate(vht->xe,&vht->ye);dCHK(err);

  {
    dInt nu,np,ne,nul,npl,nel;
    err = dFSCreateGlobalVector(vht->fsu,&vht->gvelocity);dCHK(err);
    err = dFSCreateGlobalVector(vht->fsp,&vht->gpressure);dCHK(err);
    err = dFSCreateGlobalVector(vht->fse,&vht->genergy);dCHK(err);
    err = PetscObjectSetName((PetscObject)vht->gvelocity,"MomentumDensity");dCHK(err);
    err = PetscObjectSetName((PetscObject)vht->gpressure,"Pressure");dCHK(err);
    err = PetscObjectSetName((PetscObject)vht->genergy,"EnergyDensity");dCHK(err);
    err = VecGetLocalSize(vht->gvelocity,&nu);dCHK(err);
    err = VecGetLocalSize(vht->gpressure,&np);dCHK(err);
    err = VecGetLocalSize(vht->genergy,&ne);dCHK(err);

    {                           /* Get local sizes of the closure */
      Vec  Vc,Vgh,Pc,Pgh,Ec,Egh;
      err = VecDohpGetClosure(vht->gvelocity,&Vc);dCHK(err);
      err = VecDohpGetClosure(vht->gpressure,&Pc);dCHK(err);
      err = VecDohpGetClosure(vht->genergy,&Ec);dCHK(err);
      err = VecGhostGetLocalForm(Vc,&Vgh);dCHK(err);
      err = VecGhostGetLocalForm(Pc,&Pgh);dCHK(err);
      err = VecGhostGetLocalForm(Ec,&Egh);dCHK(err);
      err = VecGetLocalSize(Vgh,&nul);dCHK(err);
      err = VecGetLocalSize(Pgh,&npl);dCHK(err);
      err = VecGetLocalSize(Egh,&nel);dCHK(err);
      err = VecGhostRestoreLocalForm(Vc,&Vgh);dCHK(err);
      err = VecGhostRestoreLocalForm(Pc,&Pgh);dCHK(err);
      err = VecGhostRestoreLocalForm(Ec,&Egh);dCHK(err);
      err = VecDohpRestoreClosure(vht->gvelocity,&Vc);dCHK(err);
      err = VecDohpRestoreClosure(vht->gpressure,&Pc);dCHK(err);
      err = VecDohpRestoreClosure(vht->genergy,&Ec);dCHK(err);
    }

    {                           /* Set up the Stokes sub-problem */
      IS   ublock,pblock;
      dInt rstart;
      err = VecCreateMPI(vht->comm,nu+np,PETSC_DETERMINE,&vht->stokes.x);dCHK(err);
      err = VecDuplicate(vht->stokes.x,&vht->stokes.y);dCHK(err);
      err = VecGetOwnershipRange(vht->stokes.x,&rstart,NULL);dCHK(err);
      err = ISCreateStride(vht->comm,nu,rstart,1,&ublock);dCHK(err);
      err = ISCreateStride(vht->comm,np,rstart+nu,1,&pblock);dCHK(err);
      err = ISSetBlockSize(ublock,3);dCHK(err);
      err = VecScatterCreate(vht->stokes.x,ublock,vht->gvelocity,NULL,&vht->stokes.extractVelocity);dCHK(err);
      err = VecScatterCreate(vht->stokes.x,pblock,vht->gpressure,NULL,&vht->stokes.extractPressure);dCHK(err);
      vht->stokes.ublock = ublock;
      vht->stokes.pblock = pblock;
      /* Create local index sets */
      err = ISCreateStride(PETSC_COMM_SELF,nul,0,1,&vht->stokes.lublock);dCHK(err);
      err = ISCreateStride(PETSC_COMM_SELF,npl,nul,1,&vht->stokes.lpblock);dCHK(err);
      err = ISSetBlockSize(vht->stokes.lublock,3);dCHK(err);
    }
    {                           /* Set up the Stokes sub-problem */
      IS   ublock,pblock,eblock,sblock;
      dInt rstart;
      err = VecCreateMPI(vht->comm,nu+np+ne,PETSC_DETERMINE,&vht->gpacked);dCHK(err);
      err = VecGetOwnershipRange(vht->gpacked,&rstart,NULL);dCHK(err);
      err = ISCreateStride(vht->comm,nu,rstart,1,&ublock);dCHK(err);
      err = ISCreateStride(vht->comm,np,rstart+nu,1,&pblock);dCHK(err);
      err = ISCreateStride(vht->comm,ne,rstart+nu+np,1,&eblock);dCHK(err);
      err = ISCreateStride(vht->comm,nu+np,rstart,1,&sblock);dCHK(err);
      err = ISSetBlockSize(ublock,3);dCHK(err);
      err = VecScatterCreate(vht->gpacked,ublock,vht->gvelocity,NULL,&vht->all.extractVelocity);dCHK(err);
      err = VecScatterCreate(vht->gpacked,pblock,vht->gpressure,NULL,&vht->all.extractPressure);dCHK(err);
      err = VecScatterCreate(vht->gpacked,eblock,vht->genergy,NULL,&vht->all.extractEnergy);dCHK(err);
      vht->all.ublock = ublock;
      vht->all.pblock = pblock;
      vht->all.eblock = eblock;
      vht->all.sblock = sblock;
      /* Create local index sets */
      err = ISCreateStride(PETSC_COMM_SELF,nul,0,1,&vht->all.lublock);dCHK(err);
      err = ISCreateStride(PETSC_COMM_SELF,npl,nul,1,&vht->all.lpblock);dCHK(err);
      err = ISCreateStride(PETSC_COMM_SELF,nel,nul+npl,1,&vht->all.leblock);dCHK(err);
      err = ISCreateStride(PETSC_COMM_SELF,nul+npl,0,1,&vht->all.lsblock);dCHK(err);
      err = ISSetBlockSize(vht->all.lublock,3);dCHK(err);
    }
  }
  err = dJacobiDestroy(&jac);dCHK(err);
  err = dMeshDestroy(&mesh);dCHK(err);

  err = VHTCaseSetType(vht->scase,scasename);dCHK(err);
  err = dFSGetBoundingBox(vht->fsu,vht->scase->bbox);dCHK(err);
  err = VHTCaseSetFromOptions(vht->scase);dCHK(err);
  if (doview) {err = VHTView(vht,NULL);dCHK(err);}
  dFunctionReturn(0);
}
static dErr VHTGetRegionIterator(VHT vht,VHTEvaluation eval,InsertMode imode,dRulesetIterator *riter)
{
  dErr err;

  dFunctionBegin;
  if (!vht->regioniter[eval]) {
    dRulesetIterator iter;
    dRuleset ruleset;
    dFS cfs;
    dMeshESH domain;
    dQuadratureMethod qmethod;
    switch (eval) {
    case EVAL_FUNCTION: qmethod = vht->function_qmethod; break;
    case EVAL_JACOBIAN: qmethod = vht->jacobian_qmethod; break;
    default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Unknown evaluation context");
    }
    err = dFSGetDomain(vht->fsu,&domain);dCHK(err);
    err = dFSGetPreferredQuadratureRuleSet(vht->fsu,domain,dTYPE_REGION,dTOPO_ALL,qmethod,&ruleset);dCHK(err);
    err = dFSGetCoordinateFS(vht->fsu,&cfs);dCHK(err);
    err = dRulesetCreateIterator(ruleset,cfs,&iter);dCHK(err);
    err = dRulesetDestroy(&ruleset);dCHK(err); /* Give ownership to iterator */
    err = dRulesetIteratorAddFS(iter,vht->fsu);dCHK(err);
    err = dRulesetIteratorAddFS(iter,vht->fsp);dCHK(err);
    err = dRulesetIteratorAddFS(iter,vht->fse);dCHK(err);
    if (eval == EVAL_FUNCTION) {err = dRulesetIteratorAddStash(iter,0,sizeof(struct VHTStash));dCHK(err);}
    vht->regioniter[eval] = iter;
  }
  *riter = vht->regioniter[eval];
  err = dRulesetIteratorSetMode(*riter,imode);dCHK(err);
  dFunctionReturn(0);
}

static dErr MatNestGetVHT(Mat A,VHT *vht)
{
  dErr err;
  Mat B,C;
  dBool flg;

  dFunctionBegin;
  err = MatNestGetSubMat(A,0,0,&B);dCHK(err);
  err = PetscObjectTypeCompare((PetscObject)B,MATNEST,&flg);dCHK(err);
  if (flg) {err = MatNestGetSubMat(B,0,0,&C);dCHK(err);}
  else C = B;
  err = MatShellGetContext(C,(void**)vht);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTExtractGlobalSplit(VHT vht,Vec X,Vec *Xu,Vec *Xp,Vec *Xe)
{
  dErr err;

  dFunctionBegin;
  if (Xu) {
    *Xu = vht->gvelocity;
    err = VecScatterBegin(vht->all.extractVelocity,X,*Xu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (vht->all.extractVelocity,X,*Xu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  if (Xp) {
    *Xp = vht->gpressure;
    err = VecScatterBegin(vht->all.extractPressure,X,*Xp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (vht->all.extractPressure,X,*Xp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  if (Xe) {
    *Xe = vht->genergy;
    err = VecScatterBegin(vht->all.extractEnergy,X,*Xe,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (vht->all.extractEnergy,X,*Xe,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  dFunctionReturn(0);
}
static dErr VHTExtractStokesSplit(VHT vht,Vec X,Vec *Xu,Vec *Xp)
{
  dErr err;

  dFunctionBegin;
  if (Xu) {
    *Xu = vht->gvelocity;
    err = VecScatterBegin(vht->stokes.extractVelocity,X,*Xu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (vht->stokes.extractVelocity,X,*Xu,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  if (Xp) {
    *Xp = vht->gpressure;
    err = VecScatterBegin(vht->stokes.extractPressure,X,*Xp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
    err = VecScatterEnd  (vht->stokes.extractPressure,X,*Xp,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr VHTCommitGlobalSplit(VHT vht,Vec *gxu,Vec *gxp,Vec *gxe,Vec gy,InsertMode imode)
{
  dErr err;

  dFunctionBegin;
  dASSERT(*gxu == vht->gvelocity);
  dASSERT(*gxp == vht->gpressure);
  dASSERT(*gxe == vht->genergy);
  err = VecScatterBegin(vht->all.extractVelocity,*gxu,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractVelocity,*gxu,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(vht->all.extractPressure,*gxp,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractPressure,*gxp,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(vht->all.extractEnergy,*gxe,gy,imode,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractEnergy,*gxe,gy,imode,SCATTER_REVERSE);dCHK(err);
  *gxu = NULL;
  *gxp = NULL;
  dFunctionReturn(0);
}

static dErr VHTDestroy(VHT *invht)
{
  VHT vht = *invht;
  dErr err;

  dFunctionBegin;
  err = dFSDestroy(&vht->fsu);dCHK(err);
  err = dFSDestroy(&vht->fsp);dCHK(err);
  err = dFSDestroy(&vht->fse);dCHK(err);
  err = VecDestroy(&vht->xu);dCHK(err);
  err = VecDestroy(&vht->yu);dCHK(err);
  err = VecDestroy(&vht->xp);dCHK(err);
  err = VecDestroy(&vht->yp);dCHK(err);
  err = VecDestroy(&vht->xe);dCHK(err);
  err = VecDestroy(&vht->ye);dCHK(err);
  err = VecDestroy(&vht->gvelocity);dCHK(err);
  err = VecDestroy(&vht->gpressure);dCHK(err);
  err = VecDestroy(&vht->genergy);dCHK(err);
  err = VecDestroy(&vht->gpacked);dCHK(err);
  {
    err = ISDestroy(&vht->stokes.ublock);dCHK(err);
    err = ISDestroy(&vht->stokes.pblock);dCHK(err);
    err = ISDestroy(&vht->stokes.lublock);dCHK(err);
    err = ISDestroy(&vht->stokes.lpblock);dCHK(err);
    err = VecScatterDestroy(&vht->stokes.extractVelocity);dCHK(err);
    err = VecScatterDestroy(&vht->stokes.extractPressure);dCHK(err);
    err = VecDestroy(&vht->stokes.x);dCHK(err);
    err = VecDestroy(&vht->stokes.y);dCHK(err);
  }
  {
    err = ISDestroy(&vht->all.ublock);dCHK(err);
    err = ISDestroy(&vht->all.pblock);dCHK(err);
    err = ISDestroy(&vht->all.eblock);dCHK(err);
    err = ISDestroy(&vht->all.sblock);dCHK(err);
    err = ISDestroy(&vht->all.lublock);dCHK(err);
    err = ISDestroy(&vht->all.lpblock);dCHK(err);
    err = ISDestroy(&vht->all.leblock);dCHK(err);
    err = ISDestroy(&vht->all.lsblock);dCHK(err);
    err = VecScatterDestroy(&vht->all.extractVelocity);dCHK(err);
    err = VecScatterDestroy(&vht->all.extractPressure);dCHK(err);
    err = VecScatterDestroy(&vht->all.extractEnergy);dCHK(err);
  }
  err = dFree(vht->log.epochs);dCHK(err);
  for (dInt i=0; i<EVAL_UB; i++) {err = dRulesetIteratorDestroy(&vht->regioniter[i]);dCHK(err);}
  err = VHTCaseDestroy(&vht->scase);dCHK(err);
  err = dFree(*invht);dCHK(err);
  dFunctionReturn(0);
}


static dErr VHTGetMatrices(VHT vht,dBool use_jblock,Mat *J,Mat *B)
{
  dErr err;
  dInt m,nu,np,ne;
  Mat Juu,Jup,Jue,Jpu,Jpp,Jpe,Jeu,Jep,Jee,Buu,Bpp,Bee;
  IS splitis[3];

  dFunctionBegin;
  err = VecGetLocalSize(vht->gpacked,&m);dCHK(err);
  err = VecGetLocalSize(vht->gvelocity,&nu);dCHK(err);
  err = VecGetLocalSize(vht->gpressure,&np);dCHK(err);
  err = VecGetLocalSize(vht->genergy,&ne);dCHK(err);

  /* Create high-order matrix for diagonal velocity block, with context \a vht */
  err = MatCreateShell(vht->comm,nu,nu,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Juu);dCHK(err);
  err = MatShellSetOperation(Juu,MATOP_GET_VECS,(void(*)(void))MatGetVecs_VHT_stokes);dCHK(err);
  err = MatShellSetOperation(Juu,MATOP_MULT,(void(*)(void))MatMult_VHT_uu);dCHK(err);
  err = MatShellSetOperation(Juu,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_VHT_uu);dCHK(err);
  /* This is a lie, it is not strictly symmetric. But it is very close to symmetric at low Reynolds number with small melt fraction */
  err = MatShellSetOperation(Juu,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_VHT_uu);dCHK(err);
  err = MatShellSetOperation(Juu,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))MatMultAdd_VHT_uu);dCHK(err);
  err = MatSetOptionsPrefix(Juu,"Juu_");dCHK(err);

  /* Create off-diagonal high-order matrix, with context \a vht */
  err = MatCreateShell(vht->comm,np,nu,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Jpu);dCHK(err);
  err = MatShellSetOperation(Jpu,MATOP_GET_VECS,(void(*)(void))MatGetVecs_VHT_stokes);dCHK(err);
  err = MatShellSetOperation(Jpu,MATOP_MULT,(void(*)(void))MatMult_VHT_pu);dCHK(err);
  err = MatShellSetOperation(Jpu,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_VHT_pu);dCHK(err);
  err = MatSetOptionsPrefix(Jpu,"Jpu_");dCHK(err);

  err = MatCreateShell(vht->comm,nu,np,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Jup);dCHK(err);
  err = MatShellSetOperation(Jup,MATOP_MULT,(void(*)(void))MatMult_VHT_up);dCHK(err);
  err = MatShellSetOperation(Jup,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_VHT_up);dCHK(err);
  err = MatSetOptionsPrefix(Jup,"Jup_");dCHK(err);

  /* These entries are really zero */
  Jpp = NULL;
  Jpe = NULL;

  err = MatCreateShell(vht->comm,nu,ne,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Jue);dCHK(err);
  err = MatShellSetOperation(Jue,MATOP_MULT,(void(*)(void))MatMult_VHT_ue);dCHK(err);
  err = MatShellSetOperation(Jue,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_VHT_ue);dCHK(err);
  err = MatSetOptionsPrefix(Jue,"Jue_");dCHK(err);

  err = MatCreateShell(vht->comm,ne,nu,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Jeu);dCHK(err);
  err = MatShellSetOperation(Jeu,MATOP_MULT,(void(*)(void))MatMult_VHT_eu);dCHK(err);
  err = MatShellSetOperation(Jeu,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_VHT_eu);dCHK(err);
  err = MatSetOptionsPrefix(Jeu,"Jeu_");dCHK(err);

  err = MatCreateShell(vht->comm,ne,np,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Jep);dCHK(err);
  err = MatShellSetOperation(Jep,MATOP_MULT,(void(*)(void))MatMult_VHT_ep);dCHK(err);
  err = MatShellSetOperation(Jep,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_VHT_ep);dCHK(err);
  err = MatSetOptionsPrefix(Jep,"Jep_");dCHK(err);

  /* Energy-energy coupling */
  err = MatCreateShell(vht->comm,ne,ne,PETSC_DETERMINE,PETSC_DETERMINE,vht,&Jee);dCHK(err);
  err = MatShellSetOperation(Jee,MATOP_GET_VECS,(void(*)(void))MatGetVecs_VHT_ee);dCHK(err);
  err = MatShellSetOperation(Jee,MATOP_MULT,(void(*)(void))MatMult_VHT_ee);dCHK(err);
  err = MatShellSetOperation(Jee,MATOP_MULT_ADD,(void(*)(void))MatMultAdd_VHT_ee);dCHK(err);
  err = MatSetOptionsPrefix(Jee,"Jee_");dCHK(err);

  /* Blocks to be used for preconditioning */
  err = dFSCreateMatrix(vht->fsu,vht->mattype_Buu,&Buu);dCHK(err);
  err = dFSCreateMatrix(vht->fsp,vht->mattype_Bpp,&Bpp);dCHK(err);
  err = dFSCreateMatrix(vht->fse,vht->mattype_Bee,&Bee);dCHK(err);
  err = MatSetOptionsPrefix(Buu,"Buu_");dCHK(err);
  err = MatSetOptionsPrefix(Bpp,"Bpp_");dCHK(err);
  err = MatSetOptionsPrefix(Bee,"Bee_");dCHK(err);
  err = MatSetOption(Bpp,MAT_SYMMETRIC,PETSC_TRUE);dCHK(err);
  err = MatSetFromOptions(Buu);dCHK(err);
  err = MatSetFromOptions(Bpp);dCHK(err);
  err = MatSetFromOptions(Bee);dCHK(err);

  splitis[0] = vht->all.ublock;
  splitis[1] = vht->all.pblock;
  splitis[2] = vht->all.eblock;
  if (vht->split_recursive) { /* Do a nested split with the Stokes block inside the overall thing */
    Mat Jss,Jse,Jes,Bss;
    IS stokesis[] = {vht->stokes.ublock,vht->stokes.pblock},allis[] = {vht->all.sblock,vht->all.eblock};
    ISLocalToGlobalMapping ltogs[2],ltogcat;
    err = MatCreateNest(vht->comm,2,stokesis,2,stokesis,((Mat[]){Juu,Jup,Jpu,Jpp}),&Jss);dCHK(err);
    err = MatCreateNest(vht->comm,2,stokesis,1,NULL,((Mat[]){Jue,Jpe}),&Jse);dCHK(err);
    err = MatCreateNest(vht->comm,1,NULL,2,stokesis,((Mat[]){Jeu,Jep}),&Jes);dCHK(err);
    err = MatCreateNest(vht->comm,2,allis,2,allis,((Mat[]){Jss,Jse,Jes,Jee}),J);dCHK(err);
    err = MatSetOptionsPrefix(Jss,"Jss_");dCHK(err);
    err = MatSetOptionsPrefix(Jse,"Jse_");dCHK(err);
    err = MatSetOptionsPrefix(Jes,"Jes_");dCHK(err);
    err = MatCreateNest(vht->comm,2,stokesis,2,stokesis,((Mat[]){Buu,NULL,NULL,Bpp}),&Bss);dCHK(err);
    err = MatSetOptionsPrefix(Bss,"Bss_");dCHK(err);
    err = MatGetLocalToGlobalMapping(Buu,&ltogs[0],NULL);dCHK(err);
    err = MatGetLocalToGlobalMapping(Bpp,&ltogs[1],NULL);dCHK(err);
    err = ISLocalToGlobalMappingConcatenate(PetscObjectComm((PetscObject)Bss),2,ltogs,&ltogcat);dCHK(err);
    err = MatSetLocalToGlobalMapping(Bss,ltogcat,ltogcat);dCHK(err);
    err = ISLocalToGlobalMappingDestroy(&ltogcat);dCHK(err);
    err = MatCreateNest(vht->comm,2,allis,2,allis,((Mat[]){Bss,NULL,NULL,Bee}),B);dCHK(err);
    err = MatDestroy(&Jss);dCHK(err);
    err = MatDestroy(&Jse);dCHK(err);
    err = MatDestroy(&Jes);dCHK(err);
    err = MatDestroy(&Bss);dCHK(err);
  } else {       /* A single top-level 3x3 split */
    err = MatCreateNest(vht->comm,3,splitis,3,splitis,((Mat[]){Juu,Jup,Jue, Jpu,Jpp,Jpe, Jeu,Jep,Jee}),J);dCHK(err);
    err = MatCreateNest(vht->comm,3,splitis,3,splitis,((Mat[]){Buu,NULL,NULL, NULL,Bpp,NULL, NULL,NULL,Bee}),B);dCHK(err);
    //err = MatNestSetSubMat(*J,1,1,Bpp);dCHK(err); // This is inconsistent, activate for certain multiplicative fieldsplit experiments
  }
  err = MatSetOptionsPrefix(*J,"J_");dCHK(err);
  err = MatSetFromOptions(*J);dCHK(err);
  if (!use_jblock) {
    err = MatShellSetOperation(*J,MATOP_MULT,(void(*)(void))MatMult_Nest_VHT_all);dCHK(err);
  }
  err = MatSetOptionsPrefix(*B,"B_");dCHK(err);
  err = MatSetFromOptions(*B);dCHK(err);

  err = MatDestroy(&Juu);dCHK(err);
  err = MatDestroy(&Jup);dCHK(err);
  err = MatDestroy(&Jue);dCHK(err);
  err = MatDestroy(&Jpu);dCHK(err);
  err = MatDestroy(&Jpp);dCHK(err);
  err = MatDestroy(&Jpe);dCHK(err);
  err = MatDestroy(&Jeu);dCHK(err);
  err = MatDestroy(&Jep);dCHK(err);
  err = MatDestroy(&Jee);dCHK(err);
  err = MatDestroy(&Buu);dCHK(err);
  err = MatDestroy(&Bpp);dCHK(err);
  err = MatDestroy(&Bee);dCHK(err);
  dFunctionReturn(0);
}

// The "physics" functions below perform forward-mode derivative propagation.
// Every argument depending on model state U is accompanied by a dual component U1.
static inline void VHTRheoSplice(dScalar a,dScalar a1,dScalar a1x,dScalar b,dScalar b1,dScalar b1x,dReal x0,dReal x01,dReal width,dScalar x,dScalar x1,dScalar *y,dScalar *y1,dScalar *y1x,dScalar *y1x1)
{ // Smooth transition from state a to state b at x0 over width
  // Propagates two derivatives:
  //   a1,b1,x01,x1 is a standard perturbation
  //   a1x,b1x are derivatives with respect to x
  dScalar
    arg = (x-x0)/width,
    arg_x = 1/width,
    f   = 1 + tanh(arg),
    f_x = (1 - dSqr(tanh(arg))) * arg_x,
    f_xx = -2 * tanh(arg) * f_x * arg_x * arg_x;
  *y = a + (b-a)/2 * f;
  *y1 = a1 + (b1-a1)/2 * f + (b - a)/2 * f_x * (x1 - x01);
  *y1x = a1x + (b1x-a1x)/2 * f + (b - a)/2 * f_x;
  // For the derivative of y1x with moment, we simplify since currently a1x, b1x are independent of x
  *y1x1 = ((b1x-a1x)/2 * f_x * (x1-x01)
           + (b1-a1)/2 * f_x + (b-a)/2 * f_xx * (x1 - x01));
}
static inline void VHTRheoSplice0(dScalar a,dScalar b,dReal width,dScalar x,dScalar *y,dScalar *y1a,dScalar *y1b,dScalar *y1x,dScalar *y2ax,dScalar *y2bx,dScalar *y2xx)
{ // Smooth transition from state a to state b at x0 over width
  dScalar
    arg = x/width,
    arg_x = 1/width,
    f   = 1 + tanh(arg),
    f_x = (1 - dSqr(tanh(arg))) * arg_x,
    f_xx = -2 * tanh(arg) * f_x * arg_x * arg_x;
  *y = a + (b-a)*0.5*f;
  *y1a = 1 - 0.5*f;
  *y1b = 0.5*f;
  *y1x = (b-a)*0.5*f_x;
  *y2ax = -0.5*f_x;
  *y2bx = 0.5*f_x;
  *y2xx = (b-a)*0.5*f_xx;
}
static dErr VHTRheoSeparate(dScalar x,dScalar w,dScalar y[2],dScalar y1x[2],dScalar y2xx[2])
{ // Separate the function f(x)=x into a smooth negative part and a positive part.
  // The characteristic width of the transition is w
  const dScalar
    w2 = w*w,
    x2 = x*x,
    theerf = erf(x/(sqrt(2*w2))),
    arg = x2/(2*w2),
    theexp = arg < 50 ? exp(-arg) : 0; // To avoid underflow
  // The function values satisfy: xneg+xpos = x
  y[0] = 0.5*x*(1 - theerf) - w*theexp/sqrt(2*PETSC_PI);
  y[1] = x - y[0]; //0.5*x*(1 + theerf) + w*theexp/sqrt(2*PETSC_PI);
  // The first derivatives satisfy: xneg1+xpos1 = 1
  y1x[0] = 0.5*(1 - theerf);
  y1x[1] = 1 - y1x[0]; //0.5*(1 + theerf);
  // The second derivatives satisy: xneg2+xpos2 = 0
  y2xx[0] = -theexp/sqrt(2*PETSC_PI*w2);
  y2xx[1] = -y2xx[0]; //theexp/sqrt(2*PETSC_PI*w2);
  return 0;
}
dUNUSED
static dErr VHTRheoSeparate2(dScalar x,dScalar dUNUSED w,dScalar y[2],dScalar y1x[2],dScalar y2xx[2])
{  // Non-positive, creates melt fraction only
  y[0] = 0;
  y[1] = x;
  y1x[0] = 0;
  y1x[1] = 1;
  y2xx[0] = 0;
  y2xx[1] = 0;
  return 0;
}
dUNUSED
static dErr VHTRheoSolveEqStateAdjoint_splice(struct VHTRheology *rheo,const dScalar rhou[3],dScalar p,dScalar E,
                                       VHTScalarD *T,VHTScalarD *omega, dScalar K[2],dScalar K1[2][2])
{ // This version provides adjoints. Note the scalar derivatives T1E and omega1E can be used to perturb the gradient
  const dScalar
    rhotmp = rheo->rhoi, // cheat
    Tm = rheo->T3 - rheo->beta_CC * p,
    Tm1p = -rheo->beta_CC,
    em = rheo->c_i * (Tm - rheo->T0),
    em1p = rheo->c_i * Tm1p;
  dScalar e,e1E,T1a,T1b,T1x,T2ax,T2bx,T2xx,T2pp,T2pE,T2EE,o1a,o1b,o1x,o2ax,o2bx,o2xx,o2pp,o2pE,o2EE;
  dScalar a,a1p,a1E,b,b1p,b1E,x,x1E,x1p;

  dFunctionBegin;
  e = (E - rheo->mask_kinetic*1/(2*rhotmp) * dDotScalar3(rhou,rhou)) / rhotmp; // Derivatives are not propagated through kinetic energy
  e1E = 1/rhotmp;
  a = rheo->T0+e/rheo->c_i;
  a1p = 0;
  a1E = e1E / rheo->c_i;
  b = Tm;
  b1p = Tm1p;
  b1E = 0;
  x = e-em;
  x1E = e1E;
  x1p = -em1p;
  VHTRheoSplice0(a, b, rheo->splice_delta, x, &T->x,&T1a,&T1b,&T1x,&T2ax,&T2bx,&T2xx);
  T->dp = T1a*a1p + T1b*b1p + T1x*x1p; // Adjoint back to p
  T->dE = T1a*a1E + T1b*b1E + T1x*x1E; // Adjoint back to E
  T2pp = T2ax*x1p*a1p + T2bx*x1p*b1p + T2xx*x1p*x1p;
  T2pE = T2ax*x1E*a1p + T2bx*x1E*b1p + T2xx*x1E*x1p;
  T2EE = T2ax*x1E*a1E + T2bx*x1E*b1E + T2xx*x1E*x1E;
  VHTAssertRange(T->x,dMax(0,rheo->T0-40),rheo->T3);

  a = 0;
  a1p = 0;
  a1E = 0;
  b = (e-em)/rheo->Latent;
  b1p = -em1p/rheo->Latent;
  b1E = e1E/rheo->Latent;
  VHTRheoSplice0(a, b, rheo->splice_delta, x, &omega->x,&o1a,&o1b,&o1x,&o2ax,&o2bx,&o2xx);
  omega->dp = o1a*a1p + o1b*b1p + o1x*x1p;
  omega->dE = o1a*a1E + o1b*b1E + o1x*x1E;
  o2pp = o2ax*x1p*a1p + o2bx*x1p*b1p + o2xx*x1p*x1p;
  o2pE = o2ax*x1E*a1p + o2bx*x1E*b1p + o2xx*x1E*x1p;
  o2EE = o2ax*x1E*a1E + o2bx*x1E*b1E + o2xx*x1E*x1E;
  VHTAssertRangeClip(omega->x,0,1);

  // Diffusivity with respect to dp and dE, flux is: -Kp dp - KE dE
  K[0] = rheo->k_T * T->dp + rheo->Latent * rheo->kappa_w * omega->dp;
  K[1] = rheo->k_T * T->dE + rheo->Latent * rheo->kappa_w * omega->dE;
  K1[0][0] = rheo->k_T * T2pp + rheo->Latent * rheo->kappa_w * o2pp;
  K1[0][1] = rheo->k_T * T2pE + rheo->Latent * rheo->kappa_w * o2pE;
  K1[1][0] = rheo->k_T * T2pE + rheo->Latent * rheo->kappa_w * o2pE;
  K1[1][1] = rheo->k_T * T2EE + rheo->Latent * rheo->kappa_w * o2EE;
  K[1] += rheo->Kstab; // Stabilization because non-monotone splice functions can create negative diffusivity
  if (K[1] < 0) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Computed negative diffusivity %G, perhaps -rheo_Kstab is needed",K[1]);
  dFunctionReturn(0);
}
static dErr VHTRheoSolveEqStateAdjoint_separate(struct VHTRheology *rheo,const dScalar rhou[3],dScalar p,dScalar E,
                                                VHTScalarD *T,VHTScalarD *omega,dScalar omega2[3],VHTScalarD *rho, dScalar K[2],dScalar K1[2][2])
{ // This version provides adjoints. Note the scalar derivatives T1E and omega1E can be used to perturb the gradient
  const dScalar
    rhotmp = rheo->rhoi, // cheat
    Tm = rheo->T3 - rheo->beta_CC * p,
    Tm1p = -rheo->beta_CC,
    em = rheo->c_i * (Tm - rheo->T0),
    em1p = rheo->c_i * Tm1p;
  dScalar e,T2pp,T2pE,T2Ep,T2EE,o2pp,o2pE,o2Ep,o2EE;
  dScalar G,G1p,G1E,Y[2],Yg[2],Ygg[2],MeltFractionFromY,r,r1p,r1E;
  VHTScalarD Kexpstabmem,*Kexpstab = &Kexpstabmem;
  dErr err;

  dFunctionBegin;
  e = (E - rheo->mask_kinetic*1/(2*rhotmp) * dDotScalar3(rhou,rhou)) / rhotmp; // Derivatives are not propagated through kinetic energy
  G = e - em;
  G1p = rheo->c_i * rheo->beta_CC;
  G1E = 1/rhotmp;
  err = VHTRheoSeparate(G,rheo->splice_delta,Y,Yg,Ygg);dCHK(err); // Separate the energy into the thermal and moisture parts

  T->x = rheo->T0 + (em+Y[0])/rheo->c_i;
  T->dp = (em1p + Yg[0]*G1p)/rheo->c_i;
  T->dE = Yg[0]*G1E/rheo->c_i;
  T2pp = Ygg[0]*G1p*G1p/rheo->c_i;
  T2pE = Ygg[0]*G1E*G1p/rheo->c_i;
  T2Ep = Ygg[0]*G1p*G1E/rheo->c_i;
  T2EE = Ygg[0]*G1E*G1E/rheo->c_i;
  VHTAssertRangeClip(T->x,dMax(0,rheo->T0-40),rheo->T3+1);

  MeltFractionFromY = rheo->rhoi / (rheo->rhow * rheo->Latent); // kg/Joule, should be: rhoi*Y[1]/(rhow*L + (rhow-rhoi)*Y[1])
  omega->x = MeltFractionFromY * Y[1];
  omega->dp = MeltFractionFromY * Yg[1]*G1p;
  omega->dE = MeltFractionFromY * Yg[1]*G1E;
  o2pp = MeltFractionFromY * Ygg[1]*G1p*G1p;
  o2pE = MeltFractionFromY * Ygg[1]*G1p*G1E;
  o2Ep = MeltFractionFromY * Ygg[1]*G1E*G1p;
  o2EE = MeltFractionFromY * Ygg[1]*G1E*G1E;
  omega2[0] = o2pp;
  omega2[1] = o2pE;
  omega2[2] = o2EE;
  VHTAssertRangeClip(omega->x,0,1);

  rho->x = (1 - rheo->mask_rho*omega->x) * rheo->rhoi + rheo->mask_rho*omega->x * rheo->rhow;
  rho->dp = (rheo->rhow - rheo->rhoi) * rheo->mask_rho*omega->dp;
  rho->dE = (rheo->rhow - rheo->rhoi) * rheo->mask_rho*omega->dE;

  r = rheo->rhoi * (1-omega->x)/rho->x; // Factor to make moisture flux relative to total velocity instead of ice velocity
  r1p = rheo->rhoi * (-omega->dp/rho->x - (1-omega->x)/dSqr(rho->x)*rho->dp);
  r1E = rheo->rhoi * (-omega->dE/rho->x - (1-omega->x)/dSqr(rho->x)*rho->dE);

  // Diffusivity with respect to dp and dE, flux is: -Kp dp - KE dE
  K[0] = rheo->k_T * T->dp + rheo->Latent * rheo->kappa_w * r * omega->dp;
  K[1] = rheo->k_T * T->dE + rheo->Latent * rheo->kappa_w * r * omega->dE;
  K1[0][0] = rheo->k_T * T2pp + rheo->Latent * rheo->kappa_w * (r1p * omega->dp + r * o2pp);
  K1[0][1] = rheo->k_T * T2pE + rheo->Latent * rheo->kappa_w * (r1E * omega->dp + r * o2pE);
  K1[1][0] = rheo->k_T * T2Ep + rheo->Latent * rheo->kappa_w * (r1p * omega->dE + r * o2Ep);
  K1[1][1] = rheo->k_T * T2EE + rheo->Latent * rheo->kappa_w * (r1E * omega->dE + r * o2EE);

  // Extra stabilizing thermal diffusion, grows exponentially for non-physical temperatures
  Kexpstab->x = rheo->expstab * exp(rheo->expstab_rate*(2*rheo->T0 - rheo->T3 - T->x) / (rheo->T3 - rheo->T0));
  Kexpstab->dp = Kexpstab->x * -rheo->expstab_rate / (rheo->T3 - rheo->T0) * T->dp;
  Kexpstab->dE = Kexpstab->x * -rheo->expstab_rate / (rheo->T3 - rheo->T0) * T->dE;
  K[0] += Kexpstab->x * T->dp;
  K[1] += Kexpstab->x * T->dE;
  K1[0][0] = Kexpstab->dp * T->dp + Kexpstab->x * T2pp;
  K1[0][1] = Kexpstab->dE * T->dp + Kexpstab->x * T2pE;
  K1[1][0] = Kexpstab->dp * T->dE + Kexpstab->x * T2Ep;
  K1[1][1] = Kexpstab->dE * T->dE + Kexpstab->x * T2EE;

  K[1] += rheo->Kstab; // Stabilization because non-monotone splice functions can create negative diffusivity
  if (K[1] < 0) dERROR(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Computed negative diffusivity %G, perhaps -rheo_Kstab is needed",K[1]);
  dFunctionReturn(0);
}
static dErr VHTRheoArrhenius(struct VHTRheology *rheo,dScalar p,const VHTScalarD *T,const VHTScalarD *omega,VHTScalarD *B)
{
  dScalar
    n        = 1./(rheo->pe-1),
    Q        = rheo->Q,
    V        = rheo->V,
    R        = rheo->R,
    T0       = rheo->T0,
    Tstar    = T->x - rheo->beta_CC*(p+rheo->p0),
    Tstar1p  = T->dp - rheo->beta_CC,
    Tstar1E  = T->dE,
    exparg   = (Q - (p+rheo->p0)*V) / (n*R*Tstar) - Q/(n*R*T0),
    exparg1p = -V/(n*R*Tstar) - (Q-(p+rheo->p0)*V)/dSqr(n*R*Tstar) * n*R*Tstar1p,
    exparg1E = - (Q-(p+rheo->p0)*V)/dSqr(n*R*Tstar) * n*R*Tstar1E,
    warg     = 1 + rheo->Bomega * omega->x,
    warg1p   = rheo->Bomega * omega->dp,
    warg1E   = rheo->Bomega * omega->dE,
    wpow     = pow(warg, -1/n),
    wpow1p   = -1/n * wpow / warg * warg1p,
    wpow1E   = -1/n * wpow / warg * warg1E;
  dFunctionBegin;
  B->x  = rheo->B0 * exp(exparg) * wpow;
  B->dp = rheo->B0 * exp(exparg) * (exparg1p*wpow + wpow1p);
  B->dE = rheo->B0 * exp(exparg) * (exparg1E*wpow + wpow1E);
  VHTAssertRangeClip(exparg,-10,10);
  dFunctionReturn(0);
}
static dErr VHTRheoViscosity(struct VHTRheology *rheo,dScalar p,VHTScalarD *T,VHTScalarD *omega,const dScalar Du[6],VHTScalarD *eta,dScalar *eta1gamma)
{
  const dScalar
    pe = rheo->pe,
    gamma_reg = dSqr(rheo->eps) + 0.5*dColonSymScalar3(Du,Du)/rheo->gamma0,
    power = pow(gamma_reg,0.5*(pe-2)),
    power1gamma = 0.5*(pe-2) * power / gamma_reg;
  VHTScalarD B;
  dErr err;

  dFunctionBegin;
  err = VHTRheoArrhenius(rheo,p,T,omega,&B);dCHK(err);
  VHTAssertRange(gamma_reg,dSqr(rheo->eps),1e20);
  eta->x = B.x * power;
  VHTAssertRangeClip(eta->x,rheo->eta_min,rheo->eta_max);
  eta->dp = B.dp * power;
  eta->dE = B.dE * power;
  *eta1gamma = B.x * power1gamma / rheo->gamma0;
  dFunctionReturn(0);
}
static void VHTStashGetRho(const struct VHTStash *st,const struct VHTRheology dUNUSED *rheo,VHTScalarD *rho)
{
  memcpy(rho,&st->rho,sizeof(*rho));
}
static void VHTStashGetDui(const struct VHTStash *st,const struct VHTRheology *rheo,const dScalar drhou[9],dScalar Dui[6])
{
  dScalar du[9];
  VHTScalarD rho;
  VHTStashGetRho(st,rheo,&rho);
  for (dInt i=0; i<9; i++) du[i] = drhou[i] / rho.x; // cheat: should include rhou . grad(1/rho)
  dTensorSymCompress3(du,Dui);
}
static void VHTStashGetDui1(const struct VHTStash *st,const struct VHTRheology *rheo,const dScalar drhou1[9],dScalar p1,dScalar E1,dScalar Dui1[6])
{
  dScalar dui[9],du1[9];
  VHTScalarD rho;
  VHTStashGetRho(st,rheo,&rho);
  dTensorSymUncompress3(st->Dui,dui); // Have Dui, want drhou[] =~ rho*dui[]
  for (dInt i=0; i<3; i++)
    for (dInt j=0; j<3; j++)
      du1[i*3+j] = (drhou1[i*3+j] / rho.x
                    - (rho.x*dui[i*3+j]) / dSqr(rho.x) * (rho.dp*p1 + rho.dE*E1));
  dTensorSymCompress3(du1,Dui1);
}
static void VHTStashGetStress1(const struct VHTStash *st,const struct VHTRheology *rheo,const dScalar drhou1[9],dScalar p1,dScalar E1,dScalar Stress1[6],dScalar *Sigma1) {
  dScalar SymStress1[6],deta1,Dui1[6];
  VHTStashGetDui1(st,rheo,drhou1,p1,E1,Dui1);
  deta1 = st->eta1gamma * dColonSymScalar3(st->Dui,Dui1) + st->eta.dp*p1 + st->eta.dE*E1; // Perturbation of eta in current direction
  for (dInt i=0; i<6; i++) SymStress1[i] = st->eta.x * Dui1[i] + deta1*st->Dui[i] - p1*(i<3);
  dTensorSymUncompress3(SymStress1,Stress1);
  if (Sigma1) *Sigma1 = 2*st->eta.x*dColonSymScalar3(st->Dui,Dui1) + deta1*dColonSymScalar3(st->Dui,st->Dui);
}
static void VHTGeomTransformVecToRef(const dReal dx[9],const dScalar u[3],const dScalar u1[3],dReal ur[3],dReal ur1[3])
{                               // Transform a physical vector (usually streamline velocity) to reference coordinates
  dReal jinv[3][3],jdet;
  dGeomInvert3(dx,&jinv[0][0],&jdet);
  for (dInt i=0; i<3; i++) {
    ur[i] = jinv[i][0]*u[0] + jinv[i][1]*u[1] + jinv[i][2]*u[2];
    if (u1) ur1[i] = jinv[i][0]*u1[0] + jinv[i][1]*u1[1] + jinv[i][2]*u1[2];
  }
}
static void VHTStashGetUH(const struct VHTStash *stash,const dReal dx[9],dReal *uh)
{                               // Estimate the cell size in the streamline direction, multiplied by the velocity.
  const dScalar *u = stash->u;
  dReal ur[3];
  VHTGeomTransformVecToRef(dx,u,NULL,ur,NULL);
  *uh = dDotScalar3(u,u) / dNormScalar3p(ur,4);
}
static dErr VHTStashGetStreamlineStabilization(const struct VHTStash *st,const struct VHTRheology *rheo,const dReal dx[9],const dScalar u1[3],dScalar *stab,dScalar *stab1)
{
  const dScalar *u = st->u;
  const dReal p = 4.;
  dReal ur[3],ur1[3],urp,urp1;
  dFunctionBegin;
  VHTGeomTransformVecToRef(dx,u,u1,ur,ur1);
  urp = rheo->small + dNormScalar3p(ur,p);
  urp1 = 1./p * urp / (rheo->small + dPowReal(dAbs(ur[0]),p) + dPowReal(dAbs(ur[1]),p) + dPowReal(dAbs(ur[2]),p))
    * (p*dPowReal(dAbs(ur[0]),p-1)*dSign(ur[0])*ur1[0] + p*dPowReal(dAbs(ur[1]),p-1)*dSign(ur[1])*ur1[1] + p*dPowReal(dAbs(ur[2]),p)*dSign(ur[2])*ur1[2]);
  *stab = rheo->supg / urp;
  *stab1 = - *stab / urp * urp1;
  dFunctionReturn(0);
}
static dErr VHTStashGetStreamlineFlux(const struct VHTStash *st,const struct VHTRheology *rheo,const dReal dx[9],dScalar flux[3])
{
  const dScalar *u = st->u,*dE = st->dE,u1[3] = {0,0,0};
  dReal stab,stab1,cross;
  dErr err;

  dFunctionBegin;
  for (dInt i=0; i<3; i++) flux[i] = 0; // Initialize here due to compiler warning bug
  err = VHTStashGetStreamlineStabilization(st,rheo,dx,u1,&stab,&stab1);dCHK(err);
  cross = rheo->supg_crosswind * stab * dDotScalar3(u,u);
  for (dInt i=0; i<3; i++) {
    for (dInt j=0; j<3; j++)
      flux[i] -= (cross + stab * u[i] * u[j]) * dE[j];
  }
  dFunctionReturn(0);
}
static dErr VHTStashGetStreamlineFlux1(const struct VHTStash *st,const struct VHTRheology *rheo,const dScalar dx[9],const dScalar u1[3],const dScalar dE1[3],dScalar flux1[3])
{
  const dScalar *u = st->u,*dE = st->dE;
  dReal cross,cross1,stab,stab1;
  dErr err;

  dFunctionBegin;
  err = VHTStashGetStreamlineStabilization(st,rheo,dx,u1,&stab,&stab1);dCHK(err);
  cross = rheo->supg_crosswind * stab * dDotScalar3(u,u);
  cross1 = rheo->supg_crosswind * (stab1 * dDotScalar3(u,u)
                                   + stab * 2 * dDotScalar3(u,u1));
  for (dInt i=0; i<3; i++) {
    flux1[i] = 0;
    for (dInt j=0; j<3; j++)
      flux1[i] -= (stab1 * u[i] * u[j] * dE[j]
                   + stab * u1[i] * u[j] * dE[j]
                   + stab * u[i] * u1[j] * dE[j]
                   + stab * u[i] * u[j] * dE1[j]
                   + cross1 * dE[j]
                   + cross * dE1[j]);
  }
  dFunctionReturn(0);
}
static dErr VHTStashGetShockStabilization(const struct VHTStash *st,const struct VHTRheology *rheo,const dReal dx[9],dScalar *kstab,dScalar *kstab1)
{ // Computes nonlinear rank-1 diffusion tensor: kstab \outer kstab
  // kstab = u \cdot dE dE / |dE|^2         interpret as projection of u in direction dE
  const dScalar *dE = st->dE,*u = st->u;
  const dReal p = 4.;
  dScalar u_par[3],ur_par[3];
  dReal urp_par;

  dFunctionBegin;
  for (dInt i=0; i<3; i++) u_par[i] = dDotScalar3(u,dE) * dE[i] / (rheo->small + dDotScalar3(dE,dE));
  VHTGeomTransformVecToRef(dx,u_par,NULL,ur_par,NULL);
  urp_par = rheo->small + dNormScalar3p(ur_par,p);
  for (dInt i=0; i<3; i++) kstab[i] = dSqrt(rheo->shockhmm / urp_par) * ur_par[i];
  if (kstab1) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"No derivative of shock stabilization");
  dFunctionReturn(0);
}
static dErr VHTStashGetShockFlux(const struct VHTStash *st,const struct VHTRheology *rheo,const dReal dx[9],dScalar flux[3])
{
  const dScalar *dE = st->dE;
  dReal kstab[3];
  dErr err;

  dFunctionBegin;
  err = VHTStashGetShockStabilization(st,rheo,dx,kstab,NULL);dCHK(err);
  for (dInt i=0; i<3; i++) {
    flux[i] = 0;
    for (dInt j=0; j<3; j++)
      flux[i] -= kstab[i] * kstab[j] * dE[j];
  }
  dFunctionReturn(0);
}
static dErr VHTPointwiseComputeStash(struct VHTRheology *rheo,const dScalar rhou[3],const dScalar drhou[9],const dScalar p[1],const dScalar dp[3],const dScalar E[1],const dScalar dE[3],struct VHTStash *st)
{
  dErr err;
  VHTScalarD rho;

  dFunctionBegin;
  memset(st,0xff,sizeof(*st));
  dMakeMemUndefined(st,sizeof(*st));
  err = VHTRheoSolveEqStateAdjoint_separate(rheo,rhou,p[0],E[0], &st->T,&st->omega,st->omega2,&st->rho,st->K,st->K1);dCHK(err);
  VHTStashGetRho(st,rheo,&rho);
  for (dInt i=0; i<3; i++) st->u[i] = rhou[i] / rho.x;
  for (dInt i=0; i<3; i++) st->dp[i] = dp[i];
  for (dInt i=0; i<3; i++) st->dE[i] = dE[i];
  VHTStashGetDui(st,rheo,drhou,st->Dui);
  st->p = p[0];
  st->E = E[0];
  err = VHTRheoViscosity(rheo,p[0],&st->T,&st->omega,st->Dui,&st->eta,&st->eta1gamma);dCHK(err);
  //err = dRealTableView(sizeof(*st)/sizeof(dReal),1,(dReal*)st,PETSC_VIEWER_STDOUT_WORLD,"stash");dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTPointwiseFunction(VHTCase scase,const dReal x[3],const dScalar dx[9],dReal weight,
                                 const dScalar rhou[3],const dScalar drhou[9],const dScalar p[1],const dScalar dp[3],const dScalar E[1],const dScalar dE[3],
                                 struct VHTStash *st, dScalar rhou_[3],dScalar drhou_[6],dScalar p_[1],dScalar E_[1],dScalar dE_[3])
{
  struct VHTRheology *rheo = &scase->rheo;
  dScalar frhou[3],fp[1],fE[1],supgflux[3],shockflux[3],heatflux[3],Sigma,symstress[6],stress[9];
  dErr err;
  VHTScalarD rho;

  dFunctionBegin;
  err = VHTPointwiseComputeStash(rheo,rhou,drhou,p,dp,E,dE,st);dCHK(err);
  VHTStashGetRho(st,rheo,&rho);
  err = VHTStashGetStreamlineFlux(st,rheo,dx,supgflux);dCHK(err);
  err = VHTStashGetShockFlux(st,rheo,dx,shockflux);dCHK(err);
  scase->forcing(scase,x,frhou,fp,fE);
  for (dInt i=0; i<3; i++) heatflux[i] = -st->K[0]*dp[i] - st->K[1]*dE[i] + supgflux[i] + shockflux[i];
  for (dInt i=0; i<6; i++) symstress[i] = st->eta.x * st->Dui[i] - (i<3)*(p[0]+rheo->p0); // eta Du - p I
  dTensorSymUncompress3(symstress,stress);
  Sigma = st->eta.x*dColonSymScalar3(st->Dui,st->Dui);                         // Strain heating
  for (dInt i=0; i<3; i++) rhou_[i] = -rheo->rscale.momentum*weight * (rho.x * rheo->gravity[i] + frhou[i]); // Momentum source terms
  for (dInt i=0; i<3; i++) for (dInt j=0; j<3; j++)
                             drhou_[i*3+j] = -rheo->rscale.momentum*weight * (rheo->mask_momtrans*rhou[i]*st->u[j] - stress[i*3+j]);
  p_[0] = -rheo->rscale.mass*weight * (drhou[0]+drhou[4]+drhou[8] + fp[0]);                        // -q tr(drhou) - forcing, note tr(drhou) = div(rhou)
  E_[0] = -rheo->rscale.energy*weight * (Sigma + fE[0]);                                             // Strain heating and thermal forcing
  E_[0] -= rheo->rscale.energy*weight * dDotScalar3(rhou,rheo->gravity);                             // Gravitational source term for energy
  for (dInt i=0; i<3; i++) dE_[i] = -rheo->rscale.energy*weight * (st->u[i]*(E[0]+rheo->mask_Ep*(p[0]+rheo->p0)) + heatflux[i]);        // Transport and diffusion
  dFunctionReturn(0);
}
static dErr VHTPointwiseJacobian(const struct VHTStash *st,const struct VHTRheology *rheo,const dScalar dx[9],dReal weight,
                                 const dScalar rhou1[3],const dScalar drhou1[9],const dScalar p1[1],const dScalar dp1[3],const dScalar E1[1],const dScalar dE1[3],
                                 dScalar rhou_[3],dScalar drhou_[9],dScalar p_[1],dScalar E_[1],dScalar dE_[3])
{
  dScalar Stress1[9],Sigma1,u1[3],rho1,supgflux1[3];
  VHTScalarD rho;

  dFunctionBegin;
  VHTStashGetRho(st,rheo,&rho);
  VHTStashGetStress1(st,rheo,drhou1,p1[0],E1[0],Stress1,&Sigma1);
  rho1 = rho.dp*p1[0] + rho.dE*E1[0];
  for (dInt i=0; i<3; i++) u1[i] = rhou1[i] / rho.x - (rho.x*st->u[i]) / dSqr(rho.x) * rho1; // perturb u = rhou/rho
  VHTStashGetStreamlineFlux1(st,rheo,dx,u1,dE1,supgflux1);

  for (dInt i=0; i<3; i++) rhou_[i] = -rheo->rscale.momentum*weight * rho1 * rheo->gravity[i];
  for (dInt i=0; i<3; i++) for (dInt j=0; j<3; j++)
                             drhou_[i*3+j] = -rheo->rscale.momentum*weight
                               * (rheo->mask_momtrans*(rhou1[i]*st->u[j] + rho.x*st->u[i]*u1[j]) - Stress1[i*3+j]);
  p_[0] = -rheo->rscale.mass*weight*(drhou1[0]+drhou1[4]+drhou1[8]);                 // -q tr(Du)
  E_[0] = -rheo->rscale.energy*weight*(Sigma1 + dDotScalar3(rhou1,rheo->gravity));
  for (dInt i=0; i<3; i++) dE_[i] = -rheo->rscale.energy*weight
                             * (st->u[i] * (E1[0]+rheo->mask_Ep*p1[0]) + u1[i] * (st->E + rheo->mask_Ep*(st->p+rheo->p0))
                                - st->K[0]*dp1[i] - st->K[1]*dE1[i]
                                - (st->K1[0][0]*p1[0] + st->K1[0][1]*E1[0])*st->dp[i]
                                - (st->K1[1][0]*p1[0] + st->K1[1][1]*E1[0])*st->dE[i]
                                + supgflux1[i]);
  dFunctionReturn(0);
}

static void VHTPointwiseJacobian_uu(const struct VHTStash *st,const struct VHTRheology *rheo,dReal weight,const dScalar rhou1[3],const dScalar drhou1[9],dScalar drhou_[9])
{
  dScalar p1 = 0,E1 = 0,Stress1[9],u1[3];
  VHTScalarD rho;
  VHTStashGetRho(st,rheo,&rho);
  VHTStashGetStress1(st,rheo,drhou1,p1,E1,Stress1,NULL);
  for (dInt i=0; i<3; i++) u1[i] = rhou1[i] / rho.x; // perturb u = rhou/rho
  for (dInt i=0; i<3; i++) for (dInt j=0; j<3; j++)
                             drhou_[i*3+j] = -rheo->rscale.momentum*weight
                               * (rheo->mask_momtrans*(rhou1[i]*st->u[j] + rho.x*st->u[i]*u1[j]) - Stress1[i*3+j]);
}

static void VHTPointwiseJacobian_pu(const struct VHTRheology *rheo,dReal weight,const dScalar drhou1[9],dScalar p_[1])
{
  p_[0] = -rheo->rscale.mass*weight*(drhou1[0]+drhou1[4]+drhou1[8]);
}

static void VHTPointwiseJacobian_up(const struct VHTStash *st,const struct VHTRheology *rheo,dReal weight,dScalar p1,dScalar rhou_[3],dScalar drhou_[9])
{
  dScalar drhou1[9] = {0},E1 = 0,Stress1[9],rho1;
  VHTScalarD rho;
  VHTStashGetStress1(st,rheo,drhou1,p1,E1,Stress1,NULL);
  VHTStashGetRho(st,rheo,&rho);
  rho1 = rho.dp*p1 + rho.dE*E1;
  for (dInt i=0; i<3; i++) rhou_[i] = -rheo->rscale.momentum*weight * rho1 * rheo->gravity[i];
  for (dInt i=0; i<9; i++) drhou_[i] = rheo->rscale.momentum*weight * Stress1[i];
}
static dErr VHTPointwiseJacobian_ue(const struct VHTStash *st,const struct VHTRheology *rheo,const dScalar dx[9],dReal weight,const dScalar E1[1],const dScalar dE1[3],dScalar rhou_[3],dScalar drhou_[9])
{
  const dScalar rhou1[3] = {0,0,0},drhou1[9] = {0,0,0,0,0,0,0,0,0},p1[1] = {0},dp1[3] = {0,0,0};
  dScalar E_[3],dE_[9],p_[1];
  dErr err;
  dFunctionBegin;
  err = VHTPointwiseJacobian(st,rheo,dx,weight,rhou1,drhou1,p1,dp1,E1,dE1,rhou_,drhou_,p_,E_,dE_);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTPointwiseJacobian_eu(const struct VHTStash *st,const struct VHTRheology *rheo,const dScalar dx[9],dReal weight,const dScalar rhou1[3],const dScalar drhou1[9],dScalar E_[1],dScalar dE_[3])
{
  const dScalar p1[1] = {0},dp1[3] = {0,0,0},E1[1] = {0},dE1[3] = {0,0,0};
  dScalar rhou_[3],drhou_[9],p_[1];
  dErr err;
  dFunctionBegin;
  err = VHTPointwiseJacobian(st,rheo,dx,weight,rhou1,drhou1,p1,dp1,E1,dE1,rhou_,drhou_,p_,E_,dE_);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTPointwiseJacobian_ep(const struct VHTStash *st,const struct VHTRheology *rheo,const dScalar dx[9],dReal weight,const dScalar p1[1],const dScalar dp1[3],dScalar E_[1],dScalar dE_[3])
{
  const dScalar rhou1[3] = {0,0,0},drhou1[9] = {0,0,0,0,0,0,0,0,0},E1[1] = {0},dE1[3] = {0,0,0};
  dScalar rhou_[3],drhou_[9],p_[1];
  dErr err;
  dFunctionBegin;
  err = VHTPointwiseJacobian(st,rheo,dx,weight,rhou1,drhou1,p1,dp1,E1,dE1,rhou_,drhou_,p_,E_,dE_);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTPointwiseJacobian_ee(const struct VHTStash *st,const struct VHTRheology *rheo,const dScalar dx[9],dReal weight,const dScalar E1[1],const dScalar dE1[3],dScalar E_[1],dScalar dE_[3])
{
  dScalar u1[3],Stress1[9],Sigma1,drhou1[9] = {0},p1=0,rho1,supgflux1[3];
  VHTScalarD rho;
  dErr err;
  dFunctionBegin;
  VHTStashGetRho(st,rheo,&rho);
  rho1 = rho.dE*E1[0];
  VHTStashGetStress1(st,rheo,drhou1,p1,E1[0],Stress1,&Sigma1);
  for (dInt i=0; i<3; i++) u1[i] = - (rho.x*st->u[i]) / dSqr(rho.x) * rho1; // perturb u = rhou/rho
  err = VHTStashGetStreamlineFlux1(st,rheo,dx,u1,dE1,supgflux1);dCHK(err);

  E_[0] = -rheo->rscale.energy*weight * Sigma1;
  for (dInt i=0; i<3; i++) dE_[i] = -rheo->rscale.energy*weight
                             * (st->u[i]*(E1[0]+p1) + u1[i]*(st->E + rheo->mask_Ep*(st->p+rheo->p0))
                                - st->K[1]*dE1[i]
                                - st->K1[0][1]*E1[0]*st->dp[i]
                                - st->K1[1][1]*E1[0]*st->dE[i]
                                + supgflux1[i]);
  dFunctionReturn(0);
}

static dErr VHTCheckDomain_Private(VHT vht,SNES snes,Vec Xp,dBool *in) {
  dErr err;
  dInt minloc;
  dReal pmin;

  PetscFunctionBegin;
  *in = dTRUE;
  err = VecMin(Xp,&minloc,&pmin);dCHK(err);
  pmin += vht->scase->rheo.p0; // Convert from perturbation to absolute
  if (pmin < 0) {
    err = PetscInfo2(snes,"Negative pressure at %D value %G\n",minloc,pmin);dCHK(err);
    if (vht->domain_error >= 1) {
      *in = PETSC_FALSE;
      err = SNESSetFunctionDomainError(snes);dCHK(err);
    }
    if (vht->domain_error == 2) dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Negative pressure at %D value %G",minloc,pmin);
  }
  dFunctionReturn(0);
}
#define VHTCheckDomain(vht,snes,Xp) do {        \
    dErr  _err;                                 \
    dBool _in;                                  \
    _err = VHTCheckDomain_Private(vht,snes,Xp,&_in);dCHK(_err); \
    if (!_in) dFunctionReturn(0);               \
  } while (0)
static dErr VecClip(Vec X,PetscReal lb,PetscReal ub)
{
  dErr err;
  dInt n;
  dScalar *x;

  dFunctionBegin;
  err = VecGetLocalSize(X,&n);dCHK(err);
  err = VecGetArray(X,&x);dCHK(err);
  for (dInt i=0; i<n; i++) {
    x[i] = PetscMax(x[i],lb);
    x[i] = PetscMin(x[i],ub);
  }
  err = VecRestoreArray(X,&x);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTProjectOntoAdmissibleSet(VHT vht,Vec dUNUSED Xu,Vec Xp,Vec Xe)
{
  dErr err;

  dFunctionBegin;
  err = VecClip(Xp,vht->clip.p[0],vht->clip.p[1]);dCHK(err);
  err = VecClip(Xe,vht->clip.E[0],vht->clip.E[1]);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTViewDHM(VHT vht,const char *fname,Vec X)
{
  Vec Xu,Xp,Xe;
  dViewer view;
  dMesh mesh;
  dErr err;

  dFunctionBegin;
  err = PetscViewerCreate(vht->comm,&view);dCHK(err);
  err = PetscViewerSetType(view,PETSCVIEWERDHM);dCHK(err);
  err = PetscViewerFileSetName(view,fname);dCHK(err);
  err = PetscViewerFileSetMode(view,FILE_MODE_WRITE);dCHK(err);
  err = dFSGetMesh(vht->fsu,&mesh);dCHK(err);dCHK(err);
  err = VHTExtractGlobalSplit(vht,X,&Xu,&Xp,&Xe);dCHK(err);
  err = dFSDirichletProject(vht->fsu,Xu,dFS_INHOMOGENEOUS);dCHK(err);
  err = dFSDirichletProject(vht->fsp,Xp,dFS_INHOMOGENEOUS);dCHK(err);
  err = dFSDirichletProject(vht->fse,Xe,dFS_INHOMOGENEOUS);dCHK(err);
  err = VecView(Xu,view);dCHK(err);
  err = VecView(Xp,view);dCHK(err);
  err = VecView(Xe,view);dCHK(err);
  err = PetscViewerDestroy(&view);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTFunctionMonitorDHM(VHT vht,SNES snes,Vec X)
{
  dErr err;
  dInt it;
  char fname[dMAX_PATH_LEN];

  dFunctionBegin;
  err = SNESGetIterationNumber(snes,&it);dCHK(err);
  err = PetscSNPrintf(fname,sizeof fname,"vht-monitor-%03d.dhm",it);dCHK(err);
  err = VHTViewDHM(vht,fname,X);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTFunction(SNES snes,Vec X,Vec Y,void *ctx)
{
  VHT              vht = ctx;
  dErr             err;
  Vec              Coords,Xu,Xp,Xe;
  dRulesetIterator iter;

  dFunctionBegin;
  err = VHTLogEpochStart(&vht->log);dCHK(err);
  err = VHTFunctionMonitorDHM(vht,snes,X);dCHK(err);
  err = VHTExtractGlobalSplit(vht,X,&Xu,&Xp,&Xe);dCHK(err);
  err = VHTProjectOntoAdmissibleSet(vht,Xu,Xp,Xe);dCHK(err);
  VHTCheckDomain(vht,snes,Xp);
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,INSERT_VALUES,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,Xu,dFS_INHOMOGENEOUS, Xp,dFS_INHOMOGENEOUS,Xp,dFS_INHOMOGENEOUS, Xe,dFS_INHOMOGENEOUS,Xe,dFS_INHOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*p)[1],(*dp)[3],(*e)[1],(*de)[3];
    dScalar (*u_)[3],(*du_)[9],(*p_)[1],(*e_)[1],(*de_)[3];
    dInt Q;
    struct VHTStash *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,&u_,&du_, &p,&dp,&p_,NULL, &e,&de,&e_,&de_);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      err = VHTPointwiseFunction(vht->scase,x[i],dx[i],jw[i], u[i],du[i],p[i],dp[i],e[i],de[i], &stash[i], u_[i],du_[i],p_[i],e_[i],de_[i]);dCHK(err);
      err = VHTLogStash(&vht->log,&vht->scase->rheo,dx[i],&stash[i]);dCHK(err);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, u_,du_, p_,NULL, e_,de_);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = VHTCommitGlobalSplit(vht,&Xu,&Xp,&Xe,Y,INSERT_VALUES);dCHK(err);
  err = VHTLogEpochEnd(&vht->log);dCHK(err);
  dFunctionReturn(0);
}

static dErr MatMult_Nest_VHT_all(Mat J,Vec X,Vec Y)
{
  VHT              vht;
  Vec              Coords,Xu,Xp,Xe;
  dRulesetIterator iter;
  dErr             err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_VHTShellMult,J,X,Y,0);dCHK(err);
  err = MatNestGetVHT(J,&vht);dCHK(err);
  err = VHTExtractGlobalSplit(vht,X,&Xu,&Xp,&Xe);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,INSERT_VALUES,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_HOMOGENEOUS,Xu,dFS_HOMOGENEOUS, Xp,dFS_HOMOGENEOUS,Xp,dFS_HOMOGENEOUS, Xe,dFS_HOMOGENEOUS,Xe,dFS_HOMOGENEOUS);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*p)[1],(*dp)[3],(*e)[1],(*de)[3];
    dScalar (*u_)[3],(*du_)[9],(*p_)[1],(*e_)[1],(*de_)[3];
    dInt Q;
    struct VHTStash *stash;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,&u_,&du_, &p,&dp,&p_,NULL, &e,&de,&e_,&de_);dCHK(err);
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      err = VHTPointwiseJacobian(&stash[i],&vht->scase->rheo,dx[i],jw[i],u[i],du[i],p[i],dp[i],e[i],de[i],u_[i],du_[i],p_[i],e_[i],de_[i]);dCHK(err);
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, u_,du_, p_,NULL, e_,de_);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = VHTCommitGlobalSplit(vht,&Xu,&Xp,&Xe,Y,INSERT_VALUES);dCHK(err);
  err = PetscLogEventEnd(LOG_VHTShellMult,J,X,Y,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr MatMultXIorA_VHT_stokes(Mat A,Vec X,Vec Y,Vec Z,InsertMode imode,VHTMultMode mmode)
{
  VHT           vht;
  dRulesetIterator iter;
  Vec              Coords;
  dErr             err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_VHTShellMult,A,X,Y,Z);dCHK(err);
  err = MatShellGetContext(A,(void**)&vht);dCHK(err);
  {  /* Check that we have correct sizes */
    dInt nu,np,ne,nx,ny;
    err = VecGetSize(vht->gvelocity,&nu);dCHK(err);
    err = VecGetSize(vht->gpressure,&np);dCHK(err);
    err = VecGetSize(vht->genergy,&ne);dCHK(err);
    err = VecGetSize(X,&nx);dCHK(err);
    err = VecGetSize(Y,&ny);dCHK(err);
    switch (mmode) {
    case VHT_MULT_UU: dASSERT(nx==nu && ny==nu); break;
    case VHT_MULT_UP: dASSERT(nx==np && ny==nu); break;
    case VHT_MULT_UE: dASSERT(nx==ne && ny==nu); break;
    case VHT_MULT_PU: dASSERT(nx==nu && ny==np); break;
    default: dERROR(PETSC_COMM_SELF,1,"Sizes do not match, unknown mult operation");
    }
  }

  switch (imode) {
  case INSERT_VALUES:
    if (Z) dERROR(vht->comm,PETSC_ERR_ARG_INCOMP,"Cannot use INSERT_VALUES and set gz");
    Z = Y;
    break;
  case ADD_VALUES:
    err = VecCopy(Y,Z);dCHK(err);
    break;
  default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"unsupported imode");
  }

  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,imode,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  switch (mmode) {
  case VHT_MULT_UU:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, X,dFS_HOMOGENEOUS,Z,dFS_HOMOGENEOUS, NULL,             NULL,              NULL,NULL);dCHK(err);
    break;
  case VHT_MULT_UP:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,             Z,dFS_HOMOGENEOUS, X,dFS_HOMOGENEOUS,NULL,              NULL,NULL);dCHK(err);
    break;
  case VHT_MULT_UE:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,             Z,dFS_HOMOGENEOUS, NULL,             NULL,              X,dFS_HOMOGENEOUS,NULL);dCHK(err);
    break;
  case VHT_MULT_PU:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, X,dFS_HOMOGENEOUS,NULL,              NULL,             Z,dFS_HOMOGENEOUS, NULL,NULL);dCHK(err);
    break;
  default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmode");
  }
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*u_)[3],(*du_)[9],*p,*p_,(*e)[1],(*de)[3];
    dInt Q;
    struct VHTStash *stash;
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    switch (mmode) {
    case VHT_MULT_UU:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,&du_, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {VHTPointwiseJacobian_uu(&stash[i],&vht->scase->rheo,jw[i],u[i],du[i],du_[i]);}
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,du_, NULL,NULL, NULL,NULL);dCHK(err);
      break;
    case VHT_MULT_UP:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,NULL,&u_,&du_, &p,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {VHTPointwiseJacobian_up(&stash[i],&vht->scase->rheo,jw[i],p[i],u_[i],du_[i]);}
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, u_,du_, NULL,NULL, NULL,NULL);dCHK(err);
      break;
    case VHT_MULT_UE:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,NULL,&u_,&du_, NULL,NULL,NULL,NULL, &e,&de,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {VHTPointwiseJacobian_ue(&stash[i],&vht->scase->rheo,dx[i],jw[i],e[i],de[i],u_[i],du_[i]);}
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, u_,du_, NULL,NULL, NULL,NULL);dCHK(err);
      break;
    case VHT_MULT_PU:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,&du,NULL,NULL, NULL,NULL,&p_,NULL, NULL,NULL,NULL,NULL);dCHK(err);
      for (dInt i=0; i<Q; i++) {VHTPointwiseJacobian_pu(&vht->scase->rheo,jw[i],du[i],&p_[i]);}
      err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,NULL, p_,NULL, NULL,NULL);dCHK(err);
      break;
    default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmode");
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = PetscLogEventEnd(LOG_VHTShellMult,A,X,Y,Z);dCHK(err);
  dFunctionReturn(0);
}

static dErr MatMult_VHT_eX(Mat A,Vec X,Vec Y,VHTMultMode mmode)
{
  dErr err;

  dFunctionBegin;
  err = VecZeroEntries(Y);dCHK(err);
  err = MatMultAdd_VHT_eX(A,X,Y,Y,mmode);dCHK(err);
  dFunctionReturn(0);
}
static dErr MatMultAdd_VHT_eX(Mat A,Vec X,Vec Y,Vec Z,VHTMultMode mmode)
{
  VHT              vht;
  dRulesetIterator iter;
  Vec              Coords;
  dErr             err;

  dFunctionBegin;
  err = PetscLogEventBegin(LOG_VHTShellMult,A,X,Z,0);dCHK(err);
  err = MatShellGetContext(A,(void**)&vht);dCHK(err);
  err = VecCopy(Y,Z);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,ADD_VALUES,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  switch (mmode) {
  case VHT_MULT_EU:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, X,dFS_HOMOGENEOUS,NULL, NULL,NULL, NULL,Z,dFS_HOMOGENEOUS);dCHK(err);
    break;
  case VHT_MULT_EP:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,NULL, X,dFS_HOMOGENEOUS,NULL, NULL,Z,dFS_HOMOGENEOUS);dCHK(err);
    break;
  case VHT_MULT_EE:
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,NULL, NULL,NULL, X,dFS_HOMOGENEOUS,Z,dFS_HOMOGENEOUS);dCHK(err);
    break;
  default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmode");
  }
  while (dRulesetIteratorHasPatch(iter)) {
    const dScalar *jw;
    dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*p)[1],(*dp)[3],(*e)[1],(*de)[3],(*e_)[1],(*de_)[3];
    dInt Q;
    struct VHTStash *stash;
    err = dRulesetIteratorGetStash(iter,NULL,&stash);dCHK(err);
    switch (mmode) {
    case VHT_MULT_EU:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,NULL, NULL,NULL,NULL,NULL, NULL,NULL,&e_,&de_);dCHK(err);
      for (dInt i=0; i<Q; i++) {err = VHTPointwiseJacobian_eu(&stash[i],&vht->scase->rheo,dx[i],jw[i],u[i],du[i],e_[i],de_[i]);dCHK(err);}
      break;
    case VHT_MULT_EP:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,NULL,NULL,NULL, &p,&dp,NULL,NULL, NULL,NULL,&e_,&de_);dCHK(err);
      for (dInt i=0; i<Q; i++) {err = VHTPointwiseJacobian_ep(&stash[i],&vht->scase->rheo,dx[i],jw[i],p[i],dp[i],e_[i],de_[i]);dCHK(err);}
      break;
    case VHT_MULT_EE:
      err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL, &e,&de,&e_,&de_);dCHK(err);
      for (dInt i=0; i<Q; i++) {err = VHTPointwiseJacobian_ee(&stash[i],&vht->scase->rheo,dx[i],jw[i],e[i],de[i],e_[i],de_[i]);dCHK(err);}
      break;
    default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Invalid mmode");
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, NULL,NULL, NULL,NULL, e_,de_);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = PetscLogEventEnd(LOG_VHTShellMult,A,X,Z,0);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTJacobianAssemble_Velocity(VHT vht,Mat Buu,Vec Mdiag,Vec X)
{
  dRulesetIterator iter;
  Vec Coords,Xu;
  dScalar *Kflat;
  dErr err;

  dFunctionBegin;
  err = VHTExtractGlobalSplit(vht,X,&Xu,NULL,NULL);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_JACOBIAN,NOT_SET_VALUES,&iter);dCHK(err);
  if (Mdiag) {
    err = VecZeroEntries(Mdiag);dCHK(err);
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,Mdiag,dFS_HOMOGENEOUS, NULL,NULL, NULL,NULL);dCHK(err);
  } else {
    err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,NULL, NULL,NULL, NULL,NULL);dCHK(err);
  }
  err = dRulesetIteratorGetMatrixSpaceSplit(iter, NULL,NULL,NULL,NULL, NULL,&Kflat,NULL,NULL, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw,*interp_flat,*deriv_flat;
    const dInt *rowcol;
    dScalar (*x)[3],(*dx)[3][3],(*u)[3],(*du)[9],(*v)[3],(*p)[1],(*dp)[3],(*e)[1],(*de)[3];
    dInt Q,P;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,&v,NULL, &p,&dp,NULL,NULL, &e,&de,NULL,NULL);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
    {                           /* Scope so that we can declare new VLA pointers for convenient assembly */
      const dReal (*interp)[P] = (const dReal(*)[P])interp_flat;
      const dReal (*deriv)[P][3] = (const dReal(*)[P][3])deriv_flat;
      dScalar (*K)[3][P][3] = (dScalar(*)[3][P][3])Kflat;
      err = PetscMemzero(K,P*3*P*3*sizeof(K[0][0][0][0]));dCHK(err);
      for (dInt q=0; q<Q; q++) {
        struct VHTStash stash;
        VHTPointwiseComputeStash(&vht->scase->rheo,u[q],du[q],p[q],dp[q],e[q],de[q],&stash);
        for (dInt j=0; j<P; j++) { /* trial functions */
          for (dInt fj=0; fj<3; fj++) {
            dScalar uu[3]={0},duu[3][3] = {{0},{0},{0}},du_[3][3];
            uu[fj] = interp[q][j];
            duu[fj][0] = deriv[q][j][0];
            duu[fj][1] = deriv[q][j][1];
            duu[fj][2] = deriv[q][j][2];
            VHTPointwiseJacobian_uu(&stash,&vht->scase->rheo,jw[q],&uu[0],&duu[0][0],&du_[0][0]);
            for (dInt i=0; i<P; i++) {
              for (dInt fi=0; fi<3; fi++) {
                K[i][fi][j][fj] += (+ deriv[q][i][0] * du_[fi][0]
                                    + deriv[q][i][1] * du_[fi][1]
                                    + deriv[q][i][2] * du_[fi][2]);
              }
            }
          }
        }
      }
      err = dFSMatSetValuesBlockedExpanded(vht->fsu,Buu,P,rowcol,P,rowcol,&K[0][0][0][0],ADD_VALUES);dCHK(err);
      for (dInt i=0; i<P; i++) {
        dScalar Mentry = 0;
        for (dInt q=0; q<Q; q++) Mentry += interp[q][i] * jw[q] * interp[q][i]; /* Integrate the diagonal entry over this element */
        v[i][0] += Mentry;
        v[i][1] += Mentry;
        v[i][2] += Mentry;
      }
    }
    err = dRulesetIteratorCommitPatchApplied(iter,INSERT_VALUES, NULL,NULL, v,NULL, NULL,NULL, NULL,NULL);dCHK(err);
    err = dRulesetIteratorRestorePatchAssembly(iter, NULL,NULL,NULL,NULL, &P,&rowcol,&interp_flat,&deriv_flat, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTJacobianAssemble_PressureEnergy(VHT vht,Mat Bpp,Mat Daux,Mat Bee,Vec X)
{
  dRulesetIterator iter;
  Vec              Coords,Xu,Xe;
  dScalar          *Kpp_flat,*Kppaux_flat,*Kee_flat;
  const dInt       *Ksizes;
  dErr             err;

  dFunctionBegin;
  /* It might seem weird to be getting velocity and energy in the pressure assembly.  The reason is that this preconditioner
  * (indeed the entire problem) is always linear in pressure.  It \e might be nonlinear in velocity and energy. */
  err = VHTExtractGlobalSplit(vht,X,&Xu,NULL,&Xe);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_JACOBIAN,NOT_SET_VALUES,&iter);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,NULL, NULL,NULL, Xe,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSplit(iter, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL, NULL,NULL,&Kpp_flat,NULL, NULL,NULL,NULL,&Kee_flat);dCHK(err);
  err = dRulesetIteratorGetMatrixSpaceSizes(iter,NULL,NULL,&Ksizes);dCHK(err);
  err = dMallocA(Ksizes[2*4+2],&Kppaux_flat);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw,*interpp_flat,*derivp_flat,*interpe_flat,*derive_flat;
    const dInt *rowcolp,*rowcole;
    dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*p)[1],(*dp)[3],(*e)[1],(*de)[3];
    dInt Q,Pp,Pe;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,NULL, &p,&dp,NULL,NULL, &e,&de,NULL,NULL);dCHK(err);
    err = dRulesetIteratorGetPatchAssembly(iter, NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL, &Pp,&rowcolp,&interpp_flat,&derivp_flat, &Pe,&rowcole,&interpe_flat,&derive_flat);dCHK(err);
    {
      const dReal (*interpp)[Pp] = (const dReal(*)[Pp])interpp_flat;
      const dReal (*derivp)[Pp][3] = (const dReal(*)[Pp][3])derivp_flat;
      const dReal (*interpe)[Pe] = (const dReal(*)[Pe])interpe_flat;
      const dReal (*derive)[Pe][3] = (const dReal(*)[Pe][3])derive_flat;
      dScalar (*Kpp)[Pp] = (dScalar(*)[Pp])Kpp_flat,(*Kppaux)[Pp] = (dScalar(*)[Pp])Kppaux_flat;
      dScalar (*Kee)[Pe] = (dScalar(*)[Pe])Kee_flat;
      err = PetscMemzero(Kpp,Pp*Pp*sizeof(Kpp[0][0]));dCHK(err);
      err = PetscMemzero(Kppaux,Pp*Pp*sizeof(Kppaux[0][0]));dCHK(err);
      err = PetscMemzero(Kee,Pe*Pe*sizeof(Kee[0][0]));dCHK(err);
      for (dInt q=0; q<Q; q++) {
        struct VHTStash stash;
        err = VHTPointwiseComputeStash(&vht->scase->rheo,u[q],du[q],p[q],dp[q],e[q],de[q],&stash);dCHK(err);
        /* Pressure-pressure Jacobians  */
        for (dInt j=0; j<Pp; j++) { /* trial functions */
          for (dInt i=0; i<Pp; i++) {
            /* Scaled mass matrx */
            Kpp[i][j] += interpp[q][i] * jw[q] * (stash.rho.x/stash.eta.x) * interpp[q][j];
            /* Neumann Laplacian */
            Kppaux[i][j] += (+ derivp[q][i][0] * jw[q] * derivp[q][j][0]
                             + derivp[q][i][1] * jw[q] * derivp[q][j][1]
                             + derivp[q][i][2] * jw[q] * derivp[q][j][2]);
          }
        }
        /* Energy-energy Jacobian */
        for (dInt j=0; j<Pe; j++) {
          const dScalar ez[1] = {interpe[q][j]},dez[3] = {derive[q][j][0],derive[q][j][1],derive[q][j][2]};
          dScalar e_[1],de_[3];
          VHTPointwiseJacobian_ee(&stash,&vht->scase->rheo,dx[q],jw[q],ez,dez,e_,de_);
          for (dInt i=0; i<Pe; i++) {
            Kee[i][j] += (interpe[q][i] * e_[0]
                          + derive[q][i][0] * de_[0]
                          + derive[q][i][1] * de_[1]
                          + derive[q][i][2] * de_[2]);
          }
        }
      }
      err = dFSMatSetValuesBlockedExpanded(vht->fsp,Bpp,Pp,rowcolp,Pp,rowcolp,&Kpp[0][0],ADD_VALUES);dCHK(err);
      if (Daux) {err = dFSMatSetValuesBlockedExpanded(vht->fsp,Daux,Pp,rowcolp,Pp,rowcolp,&Kppaux[0][0],ADD_VALUES);dCHK(err);}
      err = dFSMatSetValuesBlockedExpanded(vht->fse,Bee,Pe,rowcole,Pe,rowcole,&Kee[0][0],ADD_VALUES);dCHK(err);
    }
    err = dRulesetIteratorRestorePatchAssembly(iter, NULL,NULL,NULL,NULL, &Pp,&rowcolp,&interpp_flat,&derivp_flat, &Pe,&rowcole,&interpe_flat,&derive_flat);dCHK(err);
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = dFree(Kppaux_flat);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTJacobian(SNES dUNUSED snes,Vec X,Mat *J,Mat *B,MatStructure *structure,void *ctx)
{
  VHT vht = ctx;
  dErr err;
  Mat Bss,Buu,Bpp,Daux,Bee;
  Vec Mdiag;

  dFunctionBegin;
  if (vht->split_recursive) {
    err = MatGetLocalSubMatrix(*B,vht->all.lsblock,vht->all.lsblock,&Bss);dCHK(err);
    err = MatGetLocalSubMatrix(Bss,vht->stokes.lublock,vht->stokes.lublock,&Buu);dCHK(err);
    err = MatGetLocalSubMatrix(Bss,vht->stokes.lpblock,vht->stokes.lpblock,&Bpp);dCHK(err);
    err = MatGetLocalSubMatrix(*B,vht->all.leblock,vht->all.leblock,&Bee);dCHK(err);
  } else {
    err = MatGetLocalSubMatrix(*B,vht->all.lublock,vht->all.lublock,&Buu);dCHK(err);
    err = MatGetLocalSubMatrix(*B,vht->all.lpblock,vht->all.lpblock,&Bpp);dCHK(err);
    err = MatGetLocalSubMatrix(*B,vht->all.leblock,vht->all.leblock,&Bee);dCHK(err);
  }
  err = PetscObjectQuery((dObject)Bpp,"LSC_M_diag",(dObject*)&Mdiag);dCHK(err);
  err = PetscObjectQuery((dObject)Bpp,"LSC_L",(dObject*)&Daux);dCHK(err);
  err = MatZeroEntries(*B);dCHK(err);
  if (Daux) {err = MatZeroEntries(Daux);dCHK(err);}
  err = VHTJacobianAssemble_Velocity(vht,Buu,Mdiag,X);dCHK(err);
  err = VHTJacobianAssemble_PressureEnergy(vht,Bpp,Daux,Bee,X);dCHK(err);
  if (Daux) {
    err = MatAssemblyBegin(Daux,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd  (Daux,MAT_FINAL_ASSEMBLY);dCHK(err);
  }
  if (vht->split_recursive) {
    err = MatRestoreLocalSubMatrix(Bss,vht->stokes.lublock,vht->stokes.lublock,&Buu);dCHK(err);
    err = MatRestoreLocalSubMatrix(Bss,vht->stokes.lpblock,vht->stokes.lpblock,&Bpp);dCHK(err);
    err = MatRestoreLocalSubMatrix(*B,vht->all.lsblock,vht->all.lsblock,&Bss);dCHK(err);
    err = MatRestoreLocalSubMatrix(*B,vht->all.leblock,vht->all.leblock,&Bee);dCHK(err);
  } else {
    err = MatRestoreLocalSubMatrix(*B,vht->all.lublock,vht->all.lublock,&Buu);dCHK(err);
    err = MatRestoreLocalSubMatrix(*B,vht->all.lpblock,vht->all.lpblock,&Bpp);dCHK(err);
    err = MatRestoreLocalSubMatrix(*B,vht->all.leblock,vht->all.leblock,&Bee);dCHK(err);
  }

  /* MatNest calls assembly on the constituent pieces */
  err = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);dCHK(err);
  err = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);dCHK(err);
  if (*J != *B) {
    err = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
    err = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);dCHK(err);
  }
  *structure = SAME_NONZERO_PATTERN;
  {
    Mat Jss,Juu,Juu_ex,Jss_ex;
    PetscViewer viewer,ascviewer;
    dBool flg = dFALSE,flg_stokes = dFALSE;
    err = PetscOptionsGetBool(NULL,"-explicit_nest_mat_view_draw",&flg,NULL);dCHK(err);
    err = PetscOptionsGetBool(NULL,"-explicit_stokes_nest_mat_view_draw",&flg_stokes,NULL);dCHK(err);
    if (!flg) dFunctionReturn(0);
    err = MatNestGetSubMat(*J,0,0,&Jss);dCHK(err);
    err = MatNestGetSubMat(Jss,0,0,&Juu);dCHK(err);
    err = MatNestGetSubMat(*B,0,0,&Bss);dCHK(err);
    err = MatNestGetSubMat(Bss,0,0,&Buu);dCHK(err);
    viewer = PETSC_VIEWER_DRAW_(PetscObjectComm((PetscObject)Juu));
    ascviewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)Juu));
    err = MatComputeExplicitOperator(Juu,&Juu_ex);dCHK(err);
    err = PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_CONTOUR);dCHK(err);
    err = PetscViewerASCIIPrintf(ascviewer,"Matrix-free operator Juu\n");dCHK(err);
    err = MatView(Juu_ex,viewer);dCHK(err);
    err = PetscViewerASCIIPrintf(ascviewer,"Assembled Buu\n");dCHK(err);
    err = MatView(Buu,viewer);dCHK(err);
    err = MatAXPY_Basic(Juu_ex,-1.0,Buu,DIFFERENT_NONZERO_PATTERN);dCHK(err);
    err = PetscViewerASCIIPrintf(ascviewer,"Difference Juu - Buu\n");dCHK(err);
    err = MatView(Juu_ex,viewer);dCHK(err);
    if (flg_stokes) {
      err = PetscViewerASCIIPrintf(ascviewer,"Stokes block Jss\n");dCHK(err);
      err = MatComputeExplicitOperator(Jss,&Jss_ex);dCHK(err);
      err = MatView(Jss_ex,viewer);dCHK(err);
      err = MatDestroy(&Jss_ex);dCHK(err);
    }
    err = PetscViewerPopFormat(viewer);dCHK(err);
    err = MatDestroy(&Juu_ex);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr VHTGetPressureShift(VHT vht,Vec Xp,dScalar *pressureshift)
{
  dErr             err;
  dScalar          volume = 0,shift = 0;
  dRulesetIterator iter;
  Vec              Coords;

  dFunctionBegin;
  *pressureshift = 0;
  if (!vht->alldirichlet) dFunctionReturn(0);
   // Do a volume integral of the exact solution to that we can remove the constant pressure mode
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,NOT_SET_VALUES,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, NULL,NULL, Xp,dFS_INHOMOGENEOUS,NULL, NULL,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw;
    const dScalar (*x)[3],*p;
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,NULL,NULL,NULL, NULL,NULL,NULL,NULL, &p,NULL,NULL,NULL, NULL,NULL,NULL,NULL);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[3],duu[9],pp[1],dpp[3],ee[1],dee[3];
      err = vht->scase->solution(vht->scase,x[i],uu,duu,pp,dpp,ee,dee);dCHK(err);
      volume += jw[i];
      shift  += (pp[0] - p[i]) * jw[i]; // The computed pressure sum is zero, but the continuous integral may not be
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  *pressureshift = shift/volume;
  dFunctionReturn(0);
}

static dErr VHTErrorNorms(VHT vht,Vec X,dReal N0u[3],dReal N1u[3],dReal N0p[3],dReal N1p[3],dReal N0e[3],dReal N1e[3])
{
  dErr             err;
  Vec              Coords,Xu,Xp,Xe;
  PetscScalar      pressureshift;
  dRulesetIterator iter;

  dFunctionBegin;
  err = dNormsStart(N0u,N1u);dCHK(err);
  err = dNormsStart(N0p,N1p);dCHK(err);
  err = dNormsStart(N0e,N1e);dCHK(err);
  err = VHTExtractGlobalSplit(vht,X,&Xu,&Xp,&Xe);dCHK(err);
  err = VHTGetRegionIterator(vht,EVAL_FUNCTION,NOT_SET_VALUES,&iter);dCHK(err);
  err = dFSGetGeometryVectorExpanded(vht->fsu,&Coords);dCHK(err);
  err = VHTGetPressureShift(vht,Xp,&pressureshift);dCHK(err);
  err = dRulesetIteratorStart(iter, Coords,dFS_INHOMOGENEOUS,NULL, Xu,dFS_INHOMOGENEOUS,NULL, Xp,dFS_INHOMOGENEOUS,NULL, Xe,dFS_INHOMOGENEOUS,NULL);dCHK(err);
  while (dRulesetIteratorHasPatch(iter)) {
    const dReal *jw;
    const dScalar (*x)[3],(*dx)[9],(*u)[3],(*du)[9],(*p)[1],(*dp)[3],(*e)[1],(*de)[3];
    dInt Q;
    err = dRulesetIteratorGetPatchApplied(iter,&Q,&jw, (dScalar**)&x,(dScalar**)&dx,NULL,NULL, &u,&du,NULL,NULL, &p,&dp,NULL,NULL, &e,&de,NULL,NULL);dCHK(err);
    for (dInt i=0; i<Q; i++) {
      dScalar uu[3],duu[9],pp[1],dpp[3],ee[1],dee[3];
      err = vht->scase->solution(vht->scase,x[i],uu,duu,pp,dpp,ee,dee);dCHK(err);
      pp[0] -= pressureshift;
      err = dNormsUpdate(N0u,N1u,jw[i],3,uu,u[i],duu,du[i]);dCHK(err);
      err = dNormsUpdate(N0p,N1p,jw[i],1,pp,p[i],dpp,dp[i]);dCHK(err);
      err = dNormsUpdate(N0e,N1e,jw[i],1,ee,e[i],dee,de[i]);dCHK(err);
    }
    err = dRulesetIteratorNextPatch(iter);dCHK(err);
  }
  err = dRulesetIteratorFinish(iter);dCHK(err);
  err = dNormsFinish(N0u,N1u);dCHK(err);
  err = dNormsFinish(N0p,N1p);dCHK(err);
  err = dNormsFinish(N0e,N1e);dCHK(err);
  dFunctionReturn(0);
}

static dErr VHTVecNormsSplit(VHT vht,Vec X,NormType ntype,dReal norms[3])
{
  dErr err;
  Vec Xu,Xp,Xe;

  dFunctionBegin;
  err = VHTExtractGlobalSplit(vht,X,&Xu,&Xp,&Xe);dCHK(err);
  err = VecNorm(Xu,ntype,&norms[0]);dCHK(err);
  err = VecNorm(Xp,ntype,&norms[1]);dCHK(err);
  err = VecNorm(Xe,ntype,&norms[2]);dCHK(err);
  dFunctionReturn(0);
}
static dErr VHTVecNormsSplitStokes(VHT vht,Vec X,NormType ntype,dReal norms[2])
{
  dErr err;
  Vec Xu,Xp;

  dFunctionBegin;
  dPragmaGCC(diagnostic ignored "-Wuninitialized") // spurious
  err = VHTExtractStokesSplit(vht,X,&Xu,&Xp);dCHK(err);
  err = VecNorm(Xu,ntype,&norms[0]);dCHK(err);
  err = VecNorm(Xp,ntype,&norms[1]);dCHK(err);
  dFunctionReturn(0);
}
static dErr SNESMonitorVHTSplit(SNES snes,PetscInt its,PetscReal norm2,void *ctx)
{
  PetscViewer viewer = ctx;
  VHT vht;
  dErr err;
  Vec F;
  dReal norms[3];
  dInt tablevel;

  dFunctionBegin;
  err = SNESGetFunction(snes,&F,NULL,(void**)&vht);dCHK(err);
  err = VHTVecNormsSplit(vht,F,NORM_2,norms);dCHK(err);
  err = PetscObjectGetTabLevel((PetscObject)snes,&tablevel);dCHK(err);
  err = PetscViewerASCIIAddTab(viewer,tablevel);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"%3D SNES Function norm %12.6e  VHT norms % 12.6e % 12.6e % 12.6e\n",its,(double)norm2,norms[0],norms[1],norms[2]);dCHK(err);
  err = PetscViewerASCIISubtractTab(viewer,tablevel);dCHK(err);
  dFunctionReturn(0);
}
static dErr KSPMonitorVHTSplit(KSP ksp,PetscInt its,PetscReal norm2,void *ctx)
{
  PetscViewer viewer = ctx;
  VHT vht;
  dErr err;
  Vec rhs,work,resid;
  Mat A,B;
  dReal norms[3],scnorm,bnorm;
  char normtype[dSTR_LEN];
  KSPNormType ntype;
  PC pc;
  dInt tablevel;

  dFunctionBegin;
  err = KSPGetApplicationContext(ksp,&vht);dCHK(err);
  err = KSPGetRhs(ksp,&rhs);dCHK(err);
  err = VecDuplicate(rhs,&work);dCHK(err);
  err = KSPBuildResidual(ksp,0,work,&resid);dCHK(err);
  err = VecCopy(resid,work);dCHK(err);
  err = KSPGetPC(ksp,&pc);dCHK(err);
  err = PCGetOperators(pc,&A,&B,PETSC_NULL);dCHK(err);
  err = VecNorm(work,NORM_2,&scnorm);dCHK(err);
  err = VHTVecNormsSplit(vht,work,NORM_2,norms);dCHK(err);
  err = VecDestroy(&work);dCHK(err);
  err = VecNorm(rhs,NORM_2,&bnorm);dCHK(err);

  err = KSPGetNormType(ksp,&ntype);dCHK(err);
  err = PetscStrncpy(normtype,KSPNormTypes[ntype],sizeof normtype);dCHK(err);
  err = PetscStrtolower(normtype);dCHK(err);
  err = PetscObjectGetTabLevel((PetscObject)ksp,&tablevel);dCHK(err);
  err = PetscViewerASCIIAddTab(viewer,tablevel);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"%3D KSP %s norm % 12.6e true % 12.6e rel % 12.6e  VHT % 12.6e % 12.6e % 12.6e\n",its,normtype,norm2,scnorm,scnorm/bnorm,norms[0],norms[1],norms[2]);dCHK(err);
  err = PetscViewerASCIISubtractTab(viewer,tablevel);dCHK(err);
  dFunctionReturn(0);
}
static dErr KSPMonitorVHTStokesSplit(KSP ksp,PetscInt its,PetscReal norm2,void *ctx)
{
  PetscViewer viewer = ctx;
  VHT vht;
  dErr err;
  Vec rhs,work,resid;
  Mat A,B;
  dReal norms[2],scnorm,bnorm;
  char normtype[dSTR_LEN];
  KSPNormType ntype;
  PC pc;
  dInt tablevel;

  dFunctionBegin;
  err = KSPGetApplicationContext(ksp,&vht);dCHK(err);
  err = KSPGetRhs(ksp,&rhs);dCHK(err);
  err = VecDuplicate(rhs,&work);dCHK(err);
  err = KSPBuildResidual(ksp,0,work,&resid);dCHK(err);
  err = VecCopy(resid,work);dCHK(err);
  err = KSPGetPC(ksp,&pc);dCHK(err);
  err = PCGetOperators(pc,&A,&B,PETSC_NULL);dCHK(err);
  err = VecNorm(work,NORM_2,&scnorm);dCHK(err);
  err = VHTVecNormsSplitStokes(vht,work,NORM_2,norms);dCHK(err);
  err = VecDestroy(&work);dCHK(err);
  err = VecNorm(rhs,NORM_2,&bnorm);dCHK(err);

  err = KSPGetNormType(ksp,&ntype);dCHK(err);
  err = PetscStrncpy(normtype,KSPNormTypes[ntype],sizeof normtype);dCHK(err);
  err = PetscStrtolower(normtype);dCHK(err);
  err = PetscObjectGetTabLevel((PetscObject)ksp,&tablevel);dCHK(err);
  err = PetscViewerASCIIAddTab(viewer,tablevel);dCHK(err);
  err = PetscViewerASCIIPrintf(viewer,"%3D KSP %s norm % 12.6e true % 12.6e rel % 12.6e  VHTStokes % 12.6e % 12.6e\n",its,normtype,norm2,scnorm,scnorm/bnorm,norms[0],norms[1]);dCHK(err);
  err = PetscViewerASCIISubtractTab(viewer,tablevel);dCHK(err);
  dFunctionReturn(0);
}
static dErr SNESComputeVariableBounds_VHT(SNES snes,Vec L,Vec U)
{
  VHT vht;
  dErr err;
  Vec Lu,Lp,Le,Uu,Up,Ue;

  dFunctionBegin;
  err = SNESGetApplicationContext(snes,&vht);dCHK(err);
  err = VHTExtractGlobalSplit(vht,L,&Lu,&Lp,&Le);dCHK(err);
  err = VHTExtractGlobalSplit(vht,U,&Uu,&Up,&Ue);dCHK(err);
  err = VecSet(Lu,SNES_VI_NINF);dCHK(err);
  err = VecSet(Uu,SNES_VI_INF);dCHK(err);
  err = VecSet(Lp,SNES_VI_NINF);dCHK(err);
  err = VecSet(Up,SNES_VI_INF);dCHK(err);
  err = VecSet(Le,SNES_VI_NINF);dCHK(err);
  err = VecSet(Ue,SNES_VI_INF);dCHK(err);
  dFunctionReturn(0);
}

// This function cannot runs separately for each field because the nodal basis may be different for each field
static dErr VHTGetSolutionField_All(VHT vht,dFS fs,dInt fieldnumber,Vec *insoln)
{
  Vec      Sol,Xc,Cvecg,Cvec;
  dScalar *x;
  const dScalar *coords;
  dInt     n,bs;
  dErr     err;

  dFunctionBegin;
  *insoln = 0;
  err = dFSCreateGlobalVector(fs,&Sol);dCHK(err);
  err = VecDohpGetClosure(Sol,&Xc);dCHK(err);
  err = dFSGetNodalCoordinatesGlobal(fs,&Cvecg);dCHK(err);
  err = VecDohpGetClosure(Cvecg,&Cvec);dCHK(err);
  err = VecGetLocalSize(Xc,&n);dCHK(err);
  err = VecGetBlockSize(Xc,&bs);dCHK(err);
  {
    dInt nc;
    err = VecGetLocalSize(Cvec,&nc);dCHK(err);
    if (nc*bs != n*3) dERROR(PETSC_COMM_SELF,1,"Coordinate vector has inconsistent size");
  }
  err = VecGetArray(Xc,&x);dCHK(err);
  err = VecGetArrayRead(Cvec,&coords);dCHK(err);
  for (dInt i=0; i<n/bs; i++) {
    dScalar u_unused[3],p_unused[1],du_unused[3*3],dp_unused[3],e_unused[1],de_unused[3];
    switch (fieldnumber) {
    case 0:
      err = vht->scase->solution(vht->scase,&coords[3*i],&x[i*bs],du_unused,p_unused,dp_unused,e_unused,de_unused);dCHK(err);
      break;
    case 1:
      err = vht->scase->solution(vht->scase,&coords[3*i],u_unused,du_unused,&x[i*bs],dp_unused,e_unused,de_unused);dCHK(err);
      break;
    case 2:
      err = vht->scase->solution(vht->scase,&coords[3*i],u_unused,du_unused,p_unused,dp_unused,&x[i*bs],de_unused);dCHK(err);
      break;
    default: dERROR(vht->comm,PETSC_ERR_ARG_OUTOFRANGE,"Requested field number %D",fieldnumber);
    }
  }
  err = VecRestoreArray(Xc,&x);dCHK(err);
  err = VecRestoreArrayRead(Cvec,&coords);dCHK(err);
  err = VecDohpRestoreClosure(Cvecg,&Cvec);dCHK(err);
  err = dFSInhomogeneousDirichletCommit(fs,Xc);dCHK(err);
  err = VecDohpRestoreClosure(Sol,&Xc);dCHK(err);
  *insoln = Sol;
  dFunctionReturn(0);
}

/** Creates a solution vector, commits the closure to each FS, returns packed solution vector */
static dErr VHTGetSolutionVector(VHT vht,Vec *insoln)
{
  dErr err;
  Vec Xu,Xp,Xe,spacked;

  dFunctionBegin;
  *insoln = 0;
  err = VHTGetSolutionField_All(vht,vht->fsu,0,&Xu);dCHK(err);
  err = VHTGetSolutionField_All(vht,vht->fsp,1,&Xp);dCHK(err);
  err = VHTGetSolutionField_All(vht,vht->fse,2,&Xe);dCHK(err);
  err = VecDuplicate(vht->gpacked,&spacked);dCHK(err);
  err = VecScatterBegin(vht->all.extractVelocity,Xu,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractVelocity,Xu,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(vht->all.extractPressure,Xp,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractPressure,Xp,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterBegin(vht->all.extractEnergy,Xe,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractEnergy,Xe,spacked,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecDestroy(&Xu);dCHK(err);
  err = VecDestroy(&Xp);dCHK(err);
  err = VecDestroy(&Xe);dCHK(err);
  *insoln = spacked;
  dFunctionReturn(0);
}

static dErr VHTGetNullSpace(VHT vht,MatNullSpace *matnull)
{
  dErr err;
  Vec r;

  dFunctionBegin;
  err = VecDuplicate(vht->gpacked,&r);dCHK(err);
  err = VecZeroEntries(r);dCHK(err);
  err = VecSet(vht->gpressure,1);dCHK(err);
  err = VecScatterBegin(vht->all.extractPressure,vht->gpressure,r,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd  (vht->all.extractPressure,vht->gpressure,r,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecNormalize(r,PETSC_NULL);dCHK(err);
  err = MatNullSpaceCreate(vht->comm,dFALSE,1,&r,matnull);dCHK(err);
  err = VecDestroy(&r);dCHK(err);
  dFunctionReturn(0);
}

static dErr CheckNullSpace(SNES snes,Vec U,Vec F)
{
  Mat          mffd,J,Jp;
  dBool        isnull;
  MatStructure mstruct;
  MatNullSpace matnull;
  KSP          ksp;
  dErr         err;

  dFunctionBegin;
  err = SNESGetKSP(snes,&ksp);dCHK(err);
  err = KSPGetNullSpace(ksp,&matnull);dCHK(err);
  err = MatCreateSNESMF(snes,&mffd);dCHK(err);
  err = MatSetFromOptions(mffd);dCHK(err);
  err = SNESGetJacobian(snes,&J,&Jp,NULL,NULL);dCHK(err);
  err = MatMFFDSetBase(mffd,U,F);dCHK(err);
  err = MatNullSpaceTest(matnull,mffd,&isnull);dCHK(err);
  if (!isnull) dERROR(PETSC_COMM_SELF,1,"Vector is not in the null space of the MFFD operator");dCHK(err);
  err = MatNullSpaceTest(matnull,J,&isnull);dCHK(err);
  if (!isnull) dERROR(PETSC_COMM_SELF,1,"Vector is not in the null space of J");dCHK(err);
  err = SNESComputeJacobian(snes,U,&J,&Jp,&mstruct);dCHK(err); /* To assemble blocks of Jp */
  err = MatNullSpaceTest(matnull,Jp,&isnull);dCHK(err);
  // At present, Jp intentionally contains an auxilliary matrix in the (p,p) block. It does not have the same null space as the Jacobian so we disable the error below.
  if (false && !isnull) dERROR(PETSC_COMM_SELF,1,"Vector is not in the null space of Jp");dCHK(err);
  err = MatDestroy(&mffd);dCHK(err);
  dFunctionReturn(0);
}
static dErr ComputeExplicit(SNES snes,Vec U)
{
  Mat J,Jp,expmat,expmat_fd;
  dBool contour = dFALSE;
  MatStructure mstruct;
  dErr err;

  dFunctionBegin;
  err = SNESGetJacobian(snes,&J,&Jp,NULL,NULL);dCHK(err);
  err = SNESComputeJacobian(snes,U,&J,&Jp,&mstruct);dCHK(err);
  err = MatComputeExplicitOperator(J,&expmat);dCHK(err);
  err = MatDuplicate(expmat,MAT_DO_NOT_COPY_VALUES,&expmat_fd);dCHK(err);
  err = SNESComputeJacobianDefault(snes,U,&expmat_fd,&expmat_fd,&mstruct,NULL);dCHK(err);
  err = MatSetOptionsPrefix(expmat,"explicit_");dCHK(err);
  err = MatSetOptionsPrefix(expmat_fd,"explicit_fd_");dCHK(err);
  err = MatSetFromOptions(expmat);dCHK(err);
  err = MatSetFromOptions(expmat_fd);dCHK(err);

  err = PetscOptionsGetBool(NULL,"-mat_view_contour",&contour,NULL);dCHK(err);
  if (contour) {err = PetscViewerPushFormat(PETSC_VIEWER_DRAW_WORLD,PETSC_VIEWER_DRAW_CONTOUR);dCHK(err);}
  {
    dBool flg = dFALSE;
    err = PetscOptionsGetBool(NULL,"-explicit_mat_view",&flg,NULL);dCHK(err);
    if (flg) {
      err = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD,"###  Explicit matrix using mat-free implementation of J\n");dCHK(err);
      err = MatView(expmat,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
    }
    flg = dFALSE;
    err = PetscOptionsGetBool(NULL,"-explicit_mat_view_draw",&flg,NULL);dCHK(err);
    if (flg) {err = MatView(expmat,PETSC_VIEWER_DRAW_WORLD);dCHK(err);}
  }

  {
    dBool flg = dFALSE;
    err = PetscOptionsGetBool(NULL,"-explicit_fd_mat_view",&flg,NULL);dCHK(err);
    if (flg) {
      err = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD,"###  Explicit matrix using FD\n");dCHK(err);
      err = MatView(expmat_fd,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
    }
    flg = dFALSE;
    err = PetscOptionsGetBool(NULL,"-explicit_fd_mat_view_draw",&flg,NULL);dCHK(err);
    if (flg) {err = MatView(expmat_fd,PETSC_VIEWER_DRAW_WORLD);dCHK(err);}
  }

  err = MatAXPY(expmat,-1,expmat_fd,SAME_NONZERO_PATTERN);dCHK(err);
  {
    dBool flg = dFALSE;
    err = PetscOptionsGetBool(NULL,"-explicit_diff_mat_view",&flg,NULL);dCHK(err);
    if (flg) {
      err = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_WORLD,"###  Difference between mat-free implementation of J and FD\n");dCHK(err);
      err = MatView(expmat,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
    }
    flg = dFALSE;
    err = PetscOptionsGetBool(NULL,"-explicit_diff_mat_view_draw",&flg,NULL);dCHK(err);
    if (flg) {err = MatView(expmat,PETSC_VIEWER_DRAW_WORLD);dCHK(err);}
  }
  if (contour) {err = PetscViewerPopFormat(PETSC_VIEWER_DRAW_WORLD);dCHK(err);}
  err = MatDestroy(&expmat);dCHK(err);
  err = MatDestroy(&expmat_fd);dCHK(err);
  dFunctionReturn(0);
}

int main(int argc,char *argv[])
{
  VHT vht;
  MPI_Comm comm;
  Mat J,B;
  Vec R,X,Xsoln = NULL;
  SNES snes;
  KSP ksp;
  dBool flg,check_error,check_null,compute_explicit,use_jblock,viewdhm,snes_monitor_vht,ksp_monitor_vht,sksp_monitor_vht;
  char snesmonfilename[PETSC_MAX_PATH_LEN],kspmonfilename[PETSC_MAX_PATH_LEN],skspmonfilename[PETSC_MAX_PATH_LEN],dhmfilename[PETSC_MAX_PATH_LEN];
  dErr err;

  err = dInitialize(&argc,&argv,NULL,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  err = PetscLogEventRegister("VHTShellMult",MAT_CLASSID,&LOG_VHTShellMult);dCHK(err);

  err = VHTCaseRegisterAll();dCHK(err);
  err = VHTCreate(comm,&vht);dCHK(err);
  err = VHTSetFromOptions(vht);dCHK(err);

  err = VecDuplicate(vht->gpacked,&R);dCHK(err);
  err = VecDuplicate(R,&X);dCHK(err);

  err = PetscOptionsBegin(vht->comm,NULL,"VHT solver options",__FILE__);dCHK(err); {
    check_error = vht->scase->reality ? dFALSE : dTRUE;
    err = PetscOptionsBool("-check_error","Compute errors","",check_error,&check_error,NULL);dCHK(err);
    err = PetscOptionsBool("-use_jblock","Use blocks to apply Jacobian instead of unified (more efficient) version","",use_jblock=dFALSE,&use_jblock,NULL);dCHK(err);
    err = PetscOptionsString("-viewdhm","View the solution (filename optional)","","vht.dhm",dhmfilename,sizeof dhmfilename,&viewdhm);dCHK(err);
    if (viewdhm && !dhmfilename[0]) {err = PetscSNPrintf(dhmfilename,sizeof dhmfilename,"vht-%s.dhm",vht->scase->name);dCHK(err);}
    err = PetscOptionsBool("-check_null","Check that constant pressure really is in the null space","",check_null=dFALSE,&check_null,NULL);dCHK(err);
    err = PetscOptionsBool("-compute_explicit","Compute explicit Jacobian (only very small sizes)","",compute_explicit=dFALSE,&compute_explicit,NULL);dCHK(err);
    err = PetscOptionsString("-snes_monitor_vht","Monitor norm of function split into components","SNESMonitorSet","stdout",snesmonfilename,sizeof snesmonfilename,&snes_monitor_vht);dCHK(err);
    err = PetscOptionsString("-ksp_monitor_vht","Monitor norm of function split into components","KSPMonitorSet","stdout",kspmonfilename,sizeof kspmonfilename,&ksp_monitor_vht);dCHK(err);
    err = PetscOptionsString("-fieldsplit_s_ksp_monitor_vht","Monitor norm of function split into components","KSPMonitorSet","stdout",skspmonfilename,sizeof skspmonfilename,&sksp_monitor_vht);dCHK(err);
    err = PetscOptionsBool("-vht_function_monitor_dhm","Monitor VHT residual evaluations by writing to a file","VHTFunctionMonitorDHM",vht->function_monitor,&vht->function_monitor,NULL);dCHK(err);
  } err = PetscOptionsEnd();dCHK(err);
  err = VHTGetMatrices(vht,use_jblock,&J,&B);dCHK(err);
  err = SNESCreate(comm,&snes);dCHK(err);
  err = SNESSetFunction(snes,R,VHTFunction,vht);dCHK(err);
  err = SNESSetJacobian(snes,J,B,VHTJacobian,vht);dCHK(err);
  err = SNESSetFromOptions(snes);dCHK(err);
  err = PetscObjectTypeCompare((PetscObject)snes,SNESVINEWTONRSLS,&flg);dCHK(err);
  if (!flg) {  err = PetscObjectTypeCompare((PetscObject)snes,SNESVINEWTONSSLS,&flg);dCHK(err);}
  if (flg) {err = SNESVISetComputeVariableBounds(snes,SNESComputeVariableBounds_VHT);dCHK(err);}
  err = SNESGetKSP(snes,&ksp);dCHK(err);
  err = SNESSetApplicationContext(snes,vht);dCHK(err);
  err = KSPSetApplicationContext(ksp,vht);dCHK(err);
  if (snes_monitor_vht) {
    PetscViewer monviewer;
    err = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)snes),snesmonfilename,&monviewer);dCHK(err);
    err = SNESMonitorSet(snes,SNESMonitorVHTSplit,monviewer,(PetscErrorCode (*)(void**))PetscViewerDestroy);dCHK(err);
  }
  if (ksp_monitor_vht) {
    PetscViewer monviewer;
    err = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)ksp),kspmonfilename,&monviewer);dCHK(err);
    err = KSPMonitorSet(ksp,KSPMonitorVHTSplit,monviewer,(PetscErrorCode (*)(void**))PetscViewerDestroy);dCHK(err);
  }
  {
    PC     pc;
    err = KSPGetPC(ksp,&pc);dCHK(err);
    if (vht->split_recursive) {
      KSP *subksp,stk_ksp;
      PC stk_pc;
      dInt nsub;
      err = PCFieldSplitSetIS(pc,"s",vht->all.sblock);dCHK(err);
      err = PCFieldSplitSetIS(pc,"e",vht->all.eblock);dCHK(err);
      err = PCFieldSplitGetSubKSP(pc,&nsub,&subksp);dCHK(err);
      stk_ksp = subksp[0];
      err = KSPGetPC(stk_ksp,&stk_pc);dCHK(err);
      err = PCSetType(stk_pc,PCFIELDSPLIT);dCHK(err);
      err = PCFieldSplitSetIS(stk_pc,"u",vht->stokes.ublock);dCHK(err);
      err = PCFieldSplitSetIS(stk_pc,"p",vht->stokes.pblock);dCHK(err);
      err = KSPSetApplicationContext(stk_ksp,vht);dCHK(err);
      if (sksp_monitor_vht) {
        PetscViewer monviewer;
        err = PetscViewerASCIIOpen(PetscObjectComm((PetscObject)stk_ksp),kspmonfilename,&monviewer);dCHK(err);
        err = KSPMonitorSet(stk_ksp,KSPMonitorVHTStokesSplit,monviewer,(PetscErrorCode (*)(void**))PetscViewerDestroy);dCHK(err);
      }
      err = PetscFree(subksp);dCHK(err);
    } else {
      err = PCFieldSplitSetIS(pc,"u",vht->all.ublock);dCHK(err);
      err = PCFieldSplitSetIS(pc,"p",vht->all.pblock);dCHK(err);
      err = PCFieldSplitSetIS(pc,"e",vht->all.eblock);dCHK(err);
    }
  }
  err = VHTGetSolutionVector(vht,&Xsoln);dCHK(err);
  if (!vht->scase->reality) {
    dReal nrm;
    MatStructure mstruct;
    Vec b;
    err = VecDuplicate(X,&b);dCHK(err);
    err = VecZeroEntries(X);dCHK(err);
    err = SNESComputeFunction(snes,X,b);dCHK(err); /* -f */
    err = SNESComputeFunction(snes,Xsoln,R);dCHK(err);
    err = VHTLogReset(&vht->log);dCHK(err); // Exclude the evaluations above from the log
    err = VecNorm(R,NORM_2,&nrm);dCHK(err);
    err = dPrintf(comm,"Norm of discrete residual for exact solution %g\n",nrm);dCHK(err);
    err = SNESComputeJacobian(snes,Xsoln,&J,&B,&mstruct);dCHK(err);
    err = MatMult(J,Xsoln,R);dCHK(err);
    err = VecAXPY(R,1,b);dCHK(err); /* Jx - f */
    err = VecNorm(R,NORM_2,&nrm);dCHK(err);
    err = dPrintf(comm,"Norm of discrete linear residual at exact solution %g\n",nrm);dCHK(err);
    err = VecDestroy(&b);dCHK(err);
  }

  if (vht->alldirichlet) {                             /* Set null space */
    MatNullSpace matnull;
    err = VHTGetNullSpace(vht,&matnull);dCHK(err);
    err = KSPSetNullSpace(ksp,matnull);dCHK(err);
    if (Xsoln) {err = MatNullSpaceRemove(matnull,Xsoln,NULL);dCHK(err);}
    err = MatNullSpaceDestroy(&matnull);dCHK(err);
  }
  if (check_null || compute_explicit) {
    Vec U,F;
    err = VecDuplicate(R,&U);dCHK(err);
    err = VecDuplicate(R,&F);dCHK(err);
    err = VecSet(U,1.0);dCHK(err);
    err = SNESComputeFunction(snes,U,F);dCHK(err); /* Need base for MFFD and analytic matrix-free */
    if (check_null) {err = CheckNullSpace(snes,U,F);dCHK(err);}
    if (compute_explicit) {err = ComputeExplicit(snes,U);dCHK(err);}
    err = VecDestroy(&U);dCHK(err);
    err = VecDestroy(&F);dCHK(err);
  }
  err = VecZeroEntries(R);dCHK(err);
  err = VecZeroEntries(X);dCHK(err);
  err = SNESSolve(snes,NULL,X);dCHK(err); /* ###  SOLVE  ### */
  err = VHTLogView(&vht->log,PETSC_VIEWER_STDOUT_WORLD);dCHK(err);
  if (vht->alldirichlet) {
    MatNullSpace matnull;
    err = KSPGetNullSpace(ksp,&matnull);dCHK(err); /* does not reference */
    err = MatNullSpaceRemove(matnull,X,NULL);dCHK(err);
  }
  if (check_error) {
    dReal NAu[3],NIu[3],N0u[3],N1u[3],N0p[3],N1p[3],N0e[3],N1e[3];
    err = VHTErrorNorms(vht,X,N0u,N1u,N0p,N1p,N0e,N1e);dCHK(err);
    err = dNormsAlgebraicScaled(NAu,R);dCHK(err);
    err = VecWAXPY(R,-1,Xsoln,X);dCHK(err);
    err = dNormsAlgebraicScaled(NIu,R);dCHK(err);
    err = dPrintf(comm,"Algebraic residual        |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",NAu[0],NAu[1],NAu[2]);dCHK(err);
    err = dPrintf(comm,"Interpolation residual    |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",NIu[0],NIu[1],NIu[2]);dCHK(err);
    err = dPrintf(comm,"Integral velocity error 0 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N0u[0],N0u[1],N0u[2]);dCHK(err);
    err = dPrintf(comm,"Integral velocity error 1 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N1u[0],N1u[1],N1u[2]);dCHK(err);
    err = dPrintf(comm,"Integral pressure error 0 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N0p[0],N0p[1],N0p[2]);dCHK(err);
    err = dPrintf(comm,"Integral pressure error 1 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N1p[0],N1p[1],N1p[2]);dCHK(err);
    err = dPrintf(comm,"Integral energy error   0 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N0e[0],N0e[1],N0e[2]);dCHK(err);
    err = dPrintf(comm,"Integral energy error   1 |x|_1 %8.2e  |x|_2 %8.2e  |x|_inf %8.2e\n",N1e[0],N1e[1],N1e[2]);dCHK(err);
  }
  if (viewdhm) {err = VHTViewDHM(vht,dhmfilename,X);dCHK(err);}

  err = VecDestroy(&R);dCHK(err);
  err = VecDestroy(&X);dCHK(err);
  err = VecDestroy(&Xsoln);dCHK(err);
  err = SNESDestroy(&snes);dCHK(err);
  if (J != B) {err = MatDestroy(&J);dCHK(err);}
  err = MatDestroy(&B);dCHK(err);
  err = VHTDestroy(&vht);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
