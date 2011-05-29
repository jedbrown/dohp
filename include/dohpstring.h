#ifndef _DOHPSTRING_H
#define _DOHPSTRING_H

#include <dohptype.h>

dEXTERN_C_BEGIN

extern dErr dStrcpyS(char dest[restrict],size_t n,const char src[restrict]);
extern dErr dFilePathSplit(const char *path,dInt *slash,dInt *dot);

dEXTERN_C_END

#endif
