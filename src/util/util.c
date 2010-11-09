#include <dohpstring.h>
#include <dohp.h>

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
  dErr err;

  dFunctionBegin;
  if (n && !dest) dERROR(PETSC_COMM_SELF,1,"attempting to copy into a NULL string");
  while (--n && (*d++ = *s++)) {}
  if (!n) dERROR(PETSC_COMM_SELF,1,"String truncated");
  err = dMemzero(d,n);dCHK(err);
  dFunctionReturn(0);
}

dErr dObjectGetComm(dObject obj,MPI_Comm *comm)
{ return PetscObjectGetComm(obj,comm); }

dErr dRealTableView(dInt m,dInt n,const dReal mat[],const char *name,dViewer viewer)
{
  dBool ascii;
  dErr err;

  dFunctionBegin;
  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  for (dInt i=0; i<m; i++) {
    if (name) {
      err = PetscViewerASCIIPrintf(viewer,"%10s[%2d][%2d:%2d] ",name,i,0,n);dCHK(err);
    }
    err = PetscViewerASCIIUseTabs(viewer,PETSC_NO);dCHK(err);
    for (dInt j=0; j<n; j++) {
      err = PetscViewerASCIIPrintf(viewer," % 9.5f",mat[i*n+j]);dCHK(err);
    }
    err = PetscViewerASCIIPrintf(viewer,"\n");dCHK(err);
    err = PetscViewerASCIIUseTabs(viewer,PETSC_YES);dCHK(err);
  }
  dFunctionReturn(0);
}
