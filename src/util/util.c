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

dErr dRealTableView(dInt m,dInt n,const dReal mat[],dViewer viewer,const char *format,...)
{
  va_list Argp;
  size_t fullLen;
  char name[4096];
  dBool ascii;
  dErr err;

  dFunctionBegin;
  va_start(Argp,format);
  err = PetscVSNPrintf(name,sizeof name,format,&fullLen,Argp);dCHK(err);
  va_end(Argp);
  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  if (!mat) {
    err = PetscViewerASCIIPrintf(viewer,"%10s: null %D*%D\n",name,m,n);dCHK(err);
    dFunctionReturn(0);
  }
  for (dInt i=0; i<m; i++) {
    err = PetscViewerASCIIPrintf(viewer,"%10s[%2d][%2d:%2d] ",name,i,0,n);dCHK(err);
    err = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);dCHK(err);
    for (dInt j=0; j<n; j++) {
      err = PetscViewerASCIIPrintf(viewer," % 9.5f",mat[i*n+j]);dCHK(err);
    }
    err = PetscViewerASCIIPrintf(viewer,"\n");dCHK(err);
    err = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dIntTableView(dInt m,dInt n,const dInt mat[],dViewer viewer,const char *format,...)
{
  va_list Argp;
  size_t fullLen;
  char name[4096];
  dBool ascii;
  dErr err;

  dFunctionBegin;
  va_start(Argp,format);
  err = PetscVSNPrintf(name,sizeof name,format,&fullLen,Argp);dCHK(err);
  va_end(Argp);
  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&ascii);dCHK(err);
  if (!ascii) dFunctionReturn(0);
  if (!mat) {
    err = PetscViewerASCIIPrintf(viewer,"%10s: null %D*%D\n",name,m,n);dCHK(err);
    dFunctionReturn(0);
  }
  for (dInt i=0; i<m; i++) {
    err = PetscViewerASCIIPrintf(viewer,"%10s[%2d][%2d:%2d] ",name,i,0,n);dCHK(err);
    err = PetscViewerASCIIUseTabs(viewer,PETSC_FALSE);dCHK(err);
    for (dInt j=0; j<n; j++) {
      err = PetscViewerASCIIPrintf(viewer," %3D",mat[i*n+j]);dCHK(err);
    }
    err = PetscViewerASCIIPrintf(viewer,"\n");dCHK(err);
    err = PetscViewerASCIIUseTabs(viewer,PETSC_TRUE);dCHK(err);
  }
  dFunctionReturn(0);
}

dErr dNormsStart(dReal uerr[],dReal gerr[])
{
  dErr err;

  dFunctionBegin;
  if (uerr) {err = dMemzero(uerr,3*sizeof(uerr[0]));dCHK(err);}
  if (gerr) {err = dMemzero(gerr,3*sizeof(gerr[0]));dCHK(err);}
  dFunctionReturn(0);
}

dErr dNormsUpdate(dReal uerr[],dReal gerr[],dReal jw,dInt bs,const dScalar uu[],const dScalar u[],const dScalar duu[],const dScalar du[])
{
  dFunctionBegin;
  for (dInt i=0; i<bs; i++) {
    if (uerr && uu && u) {
      dReal r = uu[i] - u[i];
      uerr[0] += dAbs(r)*jw;
      uerr[1] += dSqr(r)*jw;
      uerr[2] = dMax(uerr[2],dAbs(r));
    }
    if (gerr && duu && du) {
      dReal dr[3] = {duu[i*3+0] - du[i*3+0],
                     duu[i*3+1] - du[i*3+1],
                     duu[i*3+2] - du[i*3+2]};
      dReal gr2 = dSqr(dr[0]) + dSqr(dr[1]) + dSqr(dr[2]),grabs = dSqrt(gr2);
      gerr[0] += grabs*jw;
      gerr[1] += gr2*jw;
      gerr[2] = dMax(gerr[2],grabs);
    }
  }
  dFunctionReturn(0);
}

dErr dNormsFinish(dReal uerr[],dReal gerr[])
{
  dFunctionBegin;
  if (uerr) uerr[1] = dSqrt(uerr[1]);
  if (gerr) gerr[1] = dSqrt(gerr[1]);
  dFunctionReturn(0);
}

dErr dNormsAlgebraicScaled(dReal norms[3],Vec r)
{
  dErr err;
  dInt n,bs;

  dFunctionBegin;
  err = VecNorm(r,NORM_1_AND_2,norms);dCHK(err);
  err = VecNorm(r,NORM_INFINITY,&norms[2]);dCHK(err);
  err = VecGetSize(r,&n);dCHK(err);
  err = VecGetBlockSize(r,&bs);dCHK(err);
  n /= bs;
  if (!n) n = 1;                /* Handle empty subdomain case */
  /* Scale norms to mimic L^p instead of l^p */
  norms[0] /= n;
  norms[1] /= sqrt(1.*n);
  for (dInt i=0; i<3; i++) if (norms[i] < 1e-14) norms[i] = 1e-99; /* limit norms so that we don't see rounding error */
  dFunctionReturn(0);
}
