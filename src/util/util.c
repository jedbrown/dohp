#include "dohptype.h"

/** This is similar to, but not the same as \c strcpy_s which is implemented in some C libraries
*
* It copies up to \c n-1 bytes from \c src into \c dest, then fills the rest of \c dest with \c 0.  It is an error if \c
* src does not fit.
*
* If \c dest has array type, a typical call would be
*
* err = dStrcpyS(dest,sizeof(dest),src);dCHK(err);
*
**/
dErr dStrcpyS(char dest[restrict],size_t n,const char src[restrict])
{
  char *restrict d = dest;
  const char *restrict s = src;

  dFunctionBegin;
  if (n && !dest) dERROR(1,"attempting to copy into a NULL string");
  while (--n && *d++ = *s++) {}
  if (!n) dERROR(1,"String truncated");
  err = dMemzero(d,n);dCHK(err);
  dFunctionReturn(0);
}

