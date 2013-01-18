#if !defined _stokescase_h
#define _stokescase_h

#include <dohpfs.h>
#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>

struct StokesRheology {
  dReal B,eps,p,nu;
};

typedef struct _n_StokesCase *StokesCase;
struct _n_StokesCase {
  MPI_Comm comm;
  dErr (*solution)(StokesCase,const dReal x[3],dScalar u[],dScalar du[],dScalar *p,dScalar dp[]);
  dErr (*forcing)(StokesCase,const dReal x[3],dScalar fu[],dScalar *fp);
  dErr (*setfromoptions)(StokesCase);
  dErr (*destroy)(StokesCase);
  struct StokesRheology rheo;
  dReal gravity;
  dReal bbox[3][2];
  dBool reality; // The "solution" is just a guess or boundary conditions
  void *data;
};
typedef dErr (*StokesCaseCreateFunction)(StokesCase);

struct StokesStore {
  dReal eta,deta;
  dReal Du[6];
};

typedef enum {EVAL_FUNCTION,EVAL_JACOBIAN, EVAL_UB} StokesEvaluation;
typedef enum {STOKES_MULT_A,STOKES_MULT_Bt,STOKES_MULT_B} StokesMultMode;

typedef struct _n_Stokes *Stokes;
struct _n_Stokes {
  MPI_Comm               comm;
  StokesCase             scase;
  dFS                    fsu,fsp;
  Vec                    xu,xp,yu,yp;
  Vec                    gvelocity,gpressure,gpacked;
  IS                     ublock,pblock;   /* Global index sets for each block */
  IS                     lublock,lpblock; /* Local index sets for each block */
  VecScatter             extractVelocity,extractPressure;
  dInt                   constBDeg,pressureCodim;
  dBool                  cardinalMass;
  char                   mattype_Ap[256],mattype_Dp[256];
  dQuadratureMethod      function_qmethod,jacobian_qmethod;
  dRulesetIterator       regioniter[EVAL_UB];
  dInt                   dirichlet[16]; /* Set numbers for Dirichlet conditions, 0 means unused */
  dBool                  alldirichlet;
};

extern PetscFunctionList StokesCaseList;
dErr StokesCaseRegister(const char *name,StokesCaseCreateFunction screate);

dErr StokesCaseRegisterAll_Exact(void); /* Defined in generated stokesexact.c */
dErr StokesCaseRegisterAll_Jako(void);  /* Defined in stokesjako.c */

#endif
