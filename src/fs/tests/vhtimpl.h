#if !defined _vhtcase_h
#define _vhtcase_h

#include <dohpfs.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>

struct VHTRheology {
  dReal B0;                     /* Leading factor term in hardness parameter */
  dReal R;                      /* Ideal gas constant, for Arrhenius relation */
  dReal Q;                      /* "Activation energy" for creep */
  dReal eps;                    /* Strain rate for which to regularize */
  dReal pe;                     /* Rheological exponent */
  dReal kappa0;                 /* Thermal conductivity in the cold part */
  dReal kappa1;                 /* Thermal conductivity in the warm part */
  dReal T0;                     /* Reference temperature corresponding to enthalpy=0, melting temperature for ice */
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
  void *data;
};
typedef dErr (*VHTCaseCreateFunction)(VHTCase);

struct VHTStash {
  dReal eta;                    /* Viscosity */
  dReal eta_gamma;              /* Derivative of eta with respect to gamma=Du:Du/2 */
  dReal eta_e;                  /* Derivative of eta with respect to enthalpy */
  dReal Du[6];                  /* Symmetrized velocity gradient */
  dReal u[3];                   /* Velocity */
  dReal kappa;                  /* Thermal conductivity at the current enthalpy */
  dReal kappa_e;                /* Derivative of thermal conductivity with respect to enthalpy */
  dReal e;                      /* Enthalpy */
  dReal de[3];
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
};

extern PetscFList VHTCaseList;
dErr VHTCaseRegister(const char *name,VHTCaseCreateFunction screate);

dErr VHTCaseRegisterAll_Exact(void); /* Defined in generated vhtexact.c */

#endif
