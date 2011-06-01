#include <dohpunits.h>
#include <dohp.h>

PetscClassId dUNITS_CLASSID;
static dBool dUnitsPackageInitialized;

struct _p_dUnits {
  PETSCHEADER(int);
  dInt nalloc;
  dUnit *list;
};
struct _n_dUnit {
  dUnits world;
  char *quantity;
  char *longname;               // Common name
  char *shortname;
  char *siname;
  dReal toCommon;               // Multiple of non-dimensional value to use common name
  dReal toSI;                   // Multiple of non-dimensional value to use SI name
  dReal expon[dUNITS_MAX];
};

const char *dUnitQuantityName(dUnit u)
{return u->quantity;}
const char *dUnitName(dUnit u)
{return u->longname;}
const char *dUnitShortName(dUnit u)
{return u->shortname;}
const char *dUnitSIName(dUnit u)
{return u->siname;}
dScalar dUnitDimensionalize(dUnit u,dScalar x)
{return x * u->toCommon;}
dScalar dUnitNonDimensionalize(dUnit u,dScalar x)
{return x / u->toCommon;}
dScalar dUnitDimensionalizeSI(dUnit u,dScalar x)
{return x * u->toSI;}
dScalar dUnitNonDimensionalizeSI(dUnit u,dScalar x)
{return x / u->toSI;}

const char *const dUnitsBaseTypes[] = {"LENGTH","MASS","TIME","TEMPERATURE","dUnitsBaseType","dUNITS_",0};
const char *const dUnitsBaseNamesSI[] = {"meter","kilogram","second","Kelvin"};
const char *const dUnitsBaseNamesShortSI[] = {"m","kg","s","K"};

static dErr dUnitsFinalizePackage(void)
{
  dFunctionBegin;
  dUnitsPackageInitialized = dFALSE;
  dFunctionReturn(0);
}
dErr dUnitsInitializePackage(const dUNUSED char path[])
{
  dErr err;

  dFunctionBegin;
  if (dUnitsPackageInitialized) dFunctionReturn(0);
  dUnitsPackageInitialized = dTRUE;
  err = PetscClassIdRegister("dUnits",&dUNITS_CLASSID);dCHK(err);
  err = PetscRegisterFinalize(dUnitsFinalizePackage);dCHK(err);
  dFunctionReturn(0);
}
dErr dUnitsCreate(MPI_Comm comm,dUnits *units)
{
  dUnits un;
  dErr err;

  dFunctionBegin;
#if !defined(PETSC_USE_DYNAMIC_LIBRARIES)
  err = dUnitsInitializePackage(PETSC_NULL);dCHK(err);
#endif
  err = PetscHeaderCreate(un,_p_dUnits,int,dUNITS_CLASSID,0,"dUnits",comm,dUnitsDestroy,dUnitsView);dCHK(err);
  un->nalloc = 128;
  err = dCallocA(un->nalloc,&un->list);dCHK(err);

  err = dUnitsSetBase(un,dUNITS_LENGTH,     "meter",   "m", 1.0,1.0,NULL);dCHK(err);
  err = dUnitsSetBase(un,dUNITS_MASS,       "kilogram","kg",1.0,1.0,NULL);dCHK(err);
  err = dUnitsSetBase(un,dUNITS_TIME,       "second",  "s", 1.0,1.0,NULL);dCHK(err);
  err = dUnitsSetBase(un,dUNITS_TEMPERATURE,"Kelvin",  "K", 1.0,1.0,NULL);dCHK(err);
  *units = un;
  dFunctionReturn(0);
}

dErr dUnitsSetFromOptions(dUnits un)
{
  dErr err;

  dFunctionBegin;
  dValidHeader(un,dUNITS_CLASSID,1);
  err = PetscOptionsBegin(((PetscObject)un)->comm,((PetscObject)un)->prefix,"Units manager","dUnits");dCHK(err);
  for (dUnitsBaseType btype = 0; btype < dUNITS_MAX; btype++) {
    char opt[256],help[256],uspec[256];
    dReal commonpersi = 1.0,scale = 1.0;
    dBool flg;
    err = PetscSNPrintf(opt,sizeof opt,"-units_%s",dUnitsBaseTypes[btype]);dCHK(err);
    err = PetscSNPrintf(uspec,sizeof uspec,"%s:%s:%f:%f",dUnitsBaseNamesSI[btype],dUnitsBaseNamesShortSI[btype],commonpersi,scale);
    err = PetscSNPrintf(help,sizeof help,"Common name:short name:one common unit of %s expressed in %s:common units per non-dimensionalized",dUnitsBaseTypes[btype],dUnitsBaseNamesSI[btype]);dCHK(err);
    err = PetscOptionsString(opt,help,"dUnitsSetBase",uspec,uspec,sizeof uspec,&flg);dCHK(err);
    if (flg) {
      char *longname,*shortname,*buf1,*buf2;
      longname = uspec;
      if (!(shortname = strchr(longname,':'))) dERROR(((dObject)un)->comm,PETSC_ERR_USER,"The field specification '%s' is ':' delimited",opt);
      *shortname++ = 0;
      if (!(buf1 = strchr(shortname,':'))) dERROR(((dObject)un)->comm,PETSC_ERR_USER,"The field specification for '%s' needs four arguments, but only two given",longname);
      *buf1++ = 0;
      if (!(buf2 = strchr(buf1,':'))) dERROR(((dObject)un)->comm,PETSC_ERR_USER,"The field specification for '%s' needs four arguments, but only three given",longname);
      *buf2++ = 0;
      if (sscanf(buf1,"%lf",&commonpersi) != 1) dERROR(((dObject)un)->comm,PETSC_ERR_USER,"Size of common unit '%s' could not be parsed from '%s'",longname,buf1);
      if (sscanf(buf2,"%lf",&scale) != 1) dERROR(((dObject)un)->comm,PETSC_ERR_USER,"Scale for common unit '%s' could not be parsed from '%s'",longname,buf2);
      err = dUnitsSetBase(un,btype,longname,shortname,commonpersi,scale,NULL);dCHK(err);
    }
  }
  err = PetscOptionsEnd();dCHK(err);
  dFunctionReturn(0);
}

dErr dUnitsSetBase(dUnits un,dUnitsBaseType type,const char *longname,const char *shortname,dReal si,dReal scale,dUnit *newunit)
{
  dUnit unit;
  dErr err;

  dFunctionBegin;
  if (!un->list[type]) {err = dCallocA(1,&un->list[type]);dCHK(err);}
  unit = un->list[type];
  err = dFree(unit->quantity);dCHK(err);
  err = dFree(unit->longname);dCHK(err);
  err = dFree(unit->shortname);dCHK(err);
  err = dFree(unit->siname);dCHK(err);
  unit->world = un;
  err = PetscStrallocpy(dUnitsBaseTypes[type],&unit->quantity);dCHK(err);
  err = PetscStrallocpy(longname,&unit->longname);dCHK(err);
  err = PetscStrallocpy(shortname,&unit->shortname);dCHK(err);
  err = PetscStrallocpy(dUnitsBaseNamesSI[type],&unit->siname);dCHK(err);
  unit->toSI = si * scale;
  unit->toCommon = scale;
  if (newunit) *newunit = unit;
  dFunctionReturn(0);
}
dErr dUnitsGetBase(dUnits un,dUnitsBaseType type,dUnit *unit)
{
  dFunctionBegin;
  *unit = un->list[type];
  dFunctionReturn(0);
}

static dErr dUnitsGetEmptyUnit_Private(dUnits un,dUnit *unit)
{
  dInt i;
  dErr err;

  dFunctionBegin;
  *unit = NULL;
  for (i=0; i<un->nalloc && un->list[i]; i++) ;
  if (i == un->nalloc) dERROR(((dObject)un)->comm,PETSC_ERR_SUP,"reallocation not implemented");
  err = dCallocA(1,&un->list[i]);dCHK(err);
  un->list[i]->world = un;
  *unit = un->list[i];
  dFunctionReturn(0);
}
typedef enum {dUNITS_NAME_LONG,dUNITS_NAME_SHORT,dUNITS_NAME_SI} dUnitsNameType;
static dErr dUnitsAssignName(dUnits un,const char *(*namer)(dUnit),const char *proposed,dInt n,const dReal expon[],char **assigned)
{
  dErr err;
  char buf[1024],*p = buf;
  dInt left = 1024;
  dFunctionBegin;
  if (proposed) {
    err = PetscStrallocpy(proposed,assigned);dCHK(err);
    dFunctionReturn(0);
  }
  buf[0] = 0;
  for (dInt i=0; i<n; i++) {
    const char *s;
    dUnit base;
    dInt len;
    if (expon[i] == 0) continue;
    err = dUnitsGetBase(un,i,&base);dCHK(err);
    s = namer(base);
    if (expon[i] == 1)
      len = snprintf(p,left,"%s ",s);
    else if (round(expon[i]) == expon[i])
      len = snprintf(p,left,"%s^%1.0f ",s,expon[i]);
    else
      len = snprintf(p,left,"%s^%f ",s,expon[i]);
    left -= len;
    p += len;
  }
  p[-1] = 0; // Kill trailing space
  err = PetscStrallocpy(buf,assigned);dCHK(err);
  dFunctionReturn(0);
}

// Logically collective
dErr dUnitsCreateUnit(dUnits un,const char *type,const char *longname,const char *shortname,dInt n,const dReal expon[],dUnit *newunit)
{
  dErr err;
  dUnit unit;

  dFunctionBegin;
  dValidHeader(un,dUNITS_CLASSID,1);
  if (n < 1 || n > dUNITS_MAX) dERROR(((dObject)un)->comm,PETSC_ERR_ARG_OUTOFRANGE,"The number of exponents %D must be positive, but no larger than %D",n,(dInt)dUNITS_MAX);
  dValidRealPointer(expon,5);
  dValidPointer(newunit,6);
  err = dUnitsGetEmptyUnit_Private(un,&unit);dCHK(err);
  err = PetscStrallocpy(type,&unit->quantity);dCHK(err);
  err = dUnitsAssignName(un,dUnitName,longname,n,expon,&unit->longname);dCHK(err);
  err = dUnitsAssignName(un,dUnitShortName,shortname,n,expon,&unit->shortname);dCHK(err);
  err = dUnitsAssignName(un,dUnitSIName,NULL,n,expon,&unit->siname);dCHK(err);
  unit->toSI = 1.0;
  unit->toCommon = 1.0;
  for (dInt i=0; i<n; i++) {
    dUnit base;
    err = dUnitsGetBase(un,i,&base);dCHK(err);
    unit->toCommon *= PetscPowScalar(dUnitDimensionalize(base,1.0),expon[i]);
    unit->toSI *= PetscPowScalar(dUnitDimensionalizeSI(base,1.0),expon[i]);
    unit->expon[i] = expon[i];
  }
  *newunit = unit;
  dFunctionReturn(0);
}
dErr dUnitsFindUnit(dUnits un,const char *name,dUnit *unit)
{
  dErr err;
  dFunctionBegin;
  *unit = NULL;
  for (dInt i=0; i<un->nalloc; i++) {
    dBool flg;
    dUnit t = un->list[i];
    err = PetscStrcmp(t->longname,name,&flg);dCHK(err);
    if (flg) {*unit = t; break;}
  }
  dFunctionReturn(0);
}
dErr dUnitsView(dUnits un,dViewer viewer)
{
  dBool iascii;
  dErr err;

  dFunctionBegin;
  dValidHeader(un,dUNITS_CLASSID,1);
  if (!viewer) {err = PetscViewerASCIIGetStdout(((dObject)un)->comm,&viewer);dCHK(err);}
  dValidHeader(viewer,PETSC_VIEWER_CLASSID,2);
  PetscCheckSameComm(un,1,viewer,2);

  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);dCHK(err);
  if (iascii) {
    err = PetscObjectPrintClassNamePrefixType((PetscObject)un,viewer,"Units Manager");dCHK(err);
    err = PetscViewerASCIIPushTab(viewer);dCHK(err);
    for (dInt i=0; i<un->nalloc && un->list[i]; i++) {
      dUnit u = un->list[i];
      err = PetscViewerASCIIPrintf(viewer,"%-12s: 1 internal unit = %10.4e %s (%s) = %10.4e %s\n",
                                   dUnitQuantityName(u),dUnitDimensionalize(u,1.0),dUnitName(u),dUnitShortName(u),dUnitDimensionalizeSI(u,1.0),dUnitSIName(u));dCHK(err);
      err = PetscViewerASCIIPrintf(viewer,"%-12s  1 %s = %10.4e %s\n","",dUnitShortName(u),dUnitDimensionalizeSI(u,dUnitNonDimensionalize(u,1.0)),dUnitSIName(u));dCHK(err);
    }
    err = PetscViewerASCIIPopTab(viewer);dCHK(err);
  } else dERROR(((dObject)un)->comm,PETSC_ERR_SUP,"Viewer type %s not supported",((PetscObject)viewer)->type_name);
  dFunctionReturn(0);
}
dErr dUnitsDestroy(dUnits *unp)
{
  dUnits un = *unp;
  dErr err;

  dFunctionBegin;
  if (!un) dFunctionReturn(0);
  PetscValidHeaderSpecific(un,dUNITS_CLASSID,1);
  if (--((PetscObject)un)->refct > 0) dFunctionReturn(0);
  for (dInt i=0; i<un->nalloc; i++) {
    dUnit u = un->list[i];
    if (u) {
      err = dFree(u->quantity);dCHK(err);
      err = dFree(u->longname);dCHK(err);
      err = dFree(u->shortname);dCHK(err);
      err = dFree(u->siname);dCHK(err);
    }
    err = dFree(un->list[i]);dCHK(err);
  }
  err = dFree(un->list);dCHK(err);
  err = PetscHeaderDestroy(unp);dCHK(err);
  dFunctionReturn(0);
}
dErr dOptionsRealUnits(const char opt[],const char text[],const char man[],dUnit unit,PetscReal defaultv,PetscReal *value,PetscBool *set)
{
  char text2[512];
  PetscReal default2,value2,one;
  dBool set2;
  dErr err;

  dFunctionBegin;
  if (unit) {
    default2 = dUnitDimensionalize(unit,defaultv);
    one = dUnitDimensionalize(unit,1);
    err = PetscSNPrintf(text2,sizeof text2,"%s [%G in %G %s]",text,defaultv,one,dUnitName(unit));dCHK(err);
    err = PetscOptionsReal(opt,text2,man,default2,&value2,&set2);dCHK(err);
    if (set2) *value = dUnitNonDimensionalize(unit,value2);
    if (set) *set = set2;
  } else {
    err = PetscSNPrintf(text2,sizeof text2,"%s [nondimensional]",text);dCHK(err);
    err = PetscOptionsReal(opt,text2,man,defaultv,value,set);dCHK(err);
  }
  dFunctionReturn(0);
}
