#if !defined _dohputil_h
#define _dohputil_h
#include "dohptype.h"

dEXTERN_C_BEGIN

extern dErr dRealTableView(dInt m,dInt n,const dReal mat[],dViewer viewer,const char *format,...);
extern dErr dIntTableView(dInt m,dInt n,const dInt mat[],dViewer viewer,const char *format,...);

dErr dNormsStart(dReal uerr[],dReal gerr[]);
dErr dNormsUpdate(dReal uerr[],dReal gerr[],dReal jw,dInt bs,const dScalar uu[],const dScalar u[],const dScalar duu[],const dScalar du[]);
dErr dNormsFinish(dReal uerr[],dReal gerr[]);

dEXTERN_C_END
#endif
