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
  dReal p0;                     /* Reference pressure, model pressure is the deviation from reference pressure */
  dReal T3;                     /* Triple point temperature */
  dReal splice_delta;           /* Characteristic width of splice */
  dReal gravity;                /* Strength of gravity in z direction (probably negative) */
  dReal Kstab;                  /* Stabilization for energy diffusion */
  dReal supg;                   /* Multiplier for SU/PG stabilization */
  dReal mask_kinetic;           /* Parameter to turn on the use of kinetic energy when computing velocity */
  dReal mask_momtrans;          /* Multiplier for the transport term in momentum balance */
  dReal mask_rho;               /* Multiplier for the true rho */
  dReal mask_Ep;                /* Multiplier for p in (E+p) term in energy equation */
  dReal eta_min,eta_max;
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
  dReal bbox[3][2];
  dBool reality; // The "solution" is just a guess or boundary conditions
  dUnits units;
  struct VHTUnitTable utable;
  char name[dNAME_LEN];
  void *data;
};
typedef dErr (*VHTCaseCreateFunction)(VHTCase);

struct VHTPack {
  dScalar rhou[3],p,E;
  dScalar drhou[9],dp[3],dE[3];
};
typedef struct {
  dReal x;
  dReal dp,dE;
} VHTScalarD;
struct VHTStash {
  dReal Dui[6];                 /* Symmetrized velocity gradient */
  dReal u[3];                   /* Total Velocity */
  dReal K[2];                   /* Total pressure and energy diffusivities, -K[0]*dp - K[1]*dE */
  dReal K1[2][2];               /* Derivatives */
  VHTScalarD T;
  VHTScalarD omega;
  dScalar omega2[3];            /* Second derivatives: o2pp, o2pE=o2Ep, o2EE */
  VHTScalarD eta;
  dReal eta1gamma;
  VHTScalarD rho;
  dReal p,E;
  dReal dE[3];
  dReal dp[3];
};

struct VHTLogEpoch {
  dReal eta[2];
  dReal cPeclet[2];
  dReal cReynolds[2];
  dReal p[2];
  dReal E[2];
  dReal T[2];
  dReal K1[2];
  dReal omega[2];
  dReal Prandtl[2];
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
  Vec      gvelocity,gpressure,genergy;
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
    VecScatter extractVelocity,extractPressure,extractEnergy,extractStokes;
  } all;
  dInt              velocityBDeg,pressureCodim,energyBDeg;
  dBool             cardinalMass;
  char              mattype_Buu[256],mattype_Bpp[256],mattype_Bee[256];
  dQuadratureMethod function_qmethod,jacobian_qmethod;
  dRulesetIterator  regioniter[EVAL_UB];
  dInt              dirichlet[16]; /* Set numbers for Dirichlet conditions, 0 means unused */
  dBool             alldirichlet;
  dInt              domain_error;
  struct VHTLog     log;
};

extern PetscFList VHTCaseList;
dErr VHTCaseRegister(const char *name,VHTCaseCreateFunction screate);

dErr VHTCaseRegisterAll_Exact(void); /* Defined in generated vhtexact.c */
dErr VHTCaseRegisterAll_Jako(void);

#endif
