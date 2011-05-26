#if !defined _vhtcase_h
#define _vhtcase_h

#include <dohpfs.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohpunits.h>
#include <dohp.h>

struct VHTRheology {
  dReal B0;                     /* Viscosity at reference strain rate and temperature */
  dReal Bomega;                 /* Softening due to water content */
  dReal R;                      /* Ideal gas constant, for Arrhenius relation */
  dReal Q;                      /* "Activation energy" for creep */
  dReal V;                      /* "Activation volume" for creep */
  dReal du0;                    /* Reference strain rate */
  dReal gamma0;                 /* Second invariant of reference strain rate */
  dReal eps;                    /* Fraction of reference strain rate at which to regularize */
  dReal pe;                     /* Rheological exponent */
  dReal k_T;                    /* Thermal conductivity in the cold part */
  dReal kappa_w;                /* Hydraulic conductivity in the warm part */
  dReal c_i;                    /* Specific heat capacity of cold part */
  dReal Latent;                 /* Latent heat of fusion */
  dReal rhoi;                   /* Density of cold part */
  dReal rhow;                   /* Density of melt */
  dReal beta_CC;                /* Clausius-Capeyron gradient */
  dReal T0;                     /* Reference temperature for definition of enthalpy and viscosity */
  dReal T3;                     /* Triple point temperature */
  dReal splice_delta;           /* Characteristic width of splice */
  dReal gravity;                /* Strength of gravity in z direction (probably negative) */
};
struct VHTUnitTable {
  dUnit Length;
  dUnit Mass;
  dUnit Time;
  dUnit Temperature;
  dUnit Density;
  dUnit Energy;
  dUnit Pressure;
  dUnit StrainRate;
  dUnit Velocity;
  dUnit Acceleration;
  dUnit Viscosity;
  dUnit Volume;
  dUnit EnergyPerTemperature;
  dUnit ThermalConductivity;
  dUnit HydroConductivity;
  dUnit Diffusivity;
  dUnit SpecificHeat;
  dUnit EnergyPerMass;
  dUnit CCGradient;
};

typedef struct _n_VHTCase *VHTCase;
struct _n_VHTCase {
  MPI_Comm comm;
  dErr (*solution)(VHTCase,const dReal x[3],dScalar u[],dScalar du[],dScalar *p,dScalar dp[],dScalar *e,dScalar *de);
  dErr (*forcing)(VHTCase,const dReal x[3],dScalar fu[],dScalar *fp,dScalar *fe);
  dErr (*setfromoptions)(VHTCase);
  dErr (*destroy)(VHTCase);
  struct VHTRheology rheo;
  dReal gravity;
  dReal bbox[3][2];
  dBool reality; // The "solution" is just a guess or boundary conditions
  dUnits units;
  struct VHTUnitTable utable;
  void *data;
};
typedef dErr (*VHTCaseCreateFunction)(VHTCase);

struct VHTPack {
  dScalar rhou[3],p,E;
  dScalar drhou[9],dp[3],dE[3];
};
struct VHTStash {
  dReal eta;                    /* Viscosity */
  dReal eta1gamma;              /* Derivative of eta with respect to gamma=Du:Du/2 */
  dReal eta1E;                  /* Derivative of eta with respect to energy */
  dReal Dui[6];                 /* Symmetrized velocity gradient */
  dReal u[3];                   /* Total Velocity */
  dReal T1E;                    /* Derivative of temperature with respect to total energy */
  dReal omega1E;                /* Derivative of melt fraction with respect to total energy */
  dReal rho;                    /* Density */
  dReal E;                      /* Energy */
  dReal dT[3];
  dReal wmom[3];                /* Momentum of the water fraction in reference frame of ice */
};

struct VHTLogEpoch {
  dReal eta[2];
  dReal cPeclet[2];
  dReal cReynolds[2];
  dReal E[2];
};
struct VHTLog {
  dInt epoch,alloc;
  dBool monitor;
  struct VHTLogEpoch *epochs;
  struct VHTLogEpoch global;
};

typedef enum {EVAL_FUNCTION,EVAL_JACOBIAN, EVAL_UB} VHTEvaluation;
typedef enum {VHT_MULT_UU,VHT_MULT_UP,VHT_MULT_PU} VHTMultMode;

typedef struct _n_VHT *VHT;
struct _n_VHT {
  MPI_Comm comm;
  VHTCase  scase;
  dFS      fsu,fsp,fse;
  Vec      xu,xp,xe,yu,yp,ye;
  Vec      gvelocity,gpressure,genthalpy;
  Vec      gpacked;
  struct {
    IS ublock,pblock;
    IS lublock,lpblock;
    VecScatter extractVelocity,extractPressure;
    Vec x,y;
  } stokes;
  struct {
    IS         ublock,pblock,eblock; /* Global index sets for each block */
    IS         lublock,lpblock,leblock; /* Local index sets for each block */
    VecScatter extractVelocity,extractPressure,extractEnthalpy,extractStokes;
  } all;
  dInt              velocityBDeg,pressureCodim,enthalpyBDeg;
  dBool             cardinalMass;
  char              mattype_Buu[256],mattype_Bpp[256],mattype_Bee[256];
  dQuadratureMethod function_qmethod,jacobian_qmethod;
  dRulesetIterator  regioniter[EVAL_UB];
  dInt              dirichlet[16]; /* Set numbers for Dirichlet conditions, 0 means unused */
  dBool             alldirichlet;
  struct VHTLog     log;
};

extern PetscFList VHTCaseList;
dErr VHTCaseRegister(const char *name,VHTCaseCreateFunction screate);

dErr VHTCaseRegisterAll_Exact(void); /* Defined in generated vhtexact.c */

#endif
