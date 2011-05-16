#if !defined _dohpunits_h
#define _dohpunits_h
#include "dohptype.h"
#include <petscviewer.h>

dEXTERN_C_BEGIN

typedef struct _p_dUnits *dUnits;
typedef struct _n_dUnit *dUnit;
extern PetscClassId dUNITS_CLASSID;
typedef enum {dUNITS_LENGTH,dUNITS_MASS,dUNITS_TIME,dUNITS_TEMPERATURE,dUNITS_MAX} dUnitsBaseType;
extern const char *const dUnitsBaseTypes[];
dErr dUnitsInitializePackage(const char path[]);
dErr dUnitsCreate(MPI_Comm comm,dUnits *);
dErr dUnitsSetFromOptions(dUnits);
dErr dUnitsSetBase(dUnits,dUnitsBaseType type,const char *longname,const char *shortname,dReal si,dReal scale,dUnit*);
dErr dUnitsGetBase(dUnits,dUnitsBaseType type,dUnit*);
dErr dUnitsCreateUnit(dUnits,const char *type,const char *longname,const char *shortname,dInt,const dReal[],dUnit*);
dErr dUnitsFindUnit(dUnits,const char *name,dUnit*);
dErr dUnitsView(dUnits,dViewer);
dErr dUnitsDestroy(dUnits*);
const char *dUnitQuantityName(dUnit);
const char *dUnitName(dUnit);
const char *dUnitShortName(dUnit);
const char *dUnitSIName(dUnit);
dScalar dUnitDimensionalize(dUnit,dScalar);
dScalar dUnitNonDimensionalize(dUnit,dScalar);
dScalar dUnitDimensionalizeSI(dUnit,dScalar);
dScalar dUnitNonDimensionalizeSI(dUnit,dScalar);

dEXTERN_C_END
#endif
