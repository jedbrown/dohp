#ifndef _DOHPSYS_H
#define _DOHPSYS_H

#include "dohptype.h"

dEXTERN_C_BEGIN

extern dErr dInitialize(int*,char***,const char*,const char*);
extern dErr dFinalize(void);

dEXTERN_C_END
#endif
