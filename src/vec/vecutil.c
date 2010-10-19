#include <dohpvec.h>
#include <dohp.h>

dErr VecCreateRedimensioned(Vec X,dInt bs,Vec *Y)
{
  dErr err;
  dInt n,xbs;

  dFunctionBegin;
  dValidHeader(X,VEC_CLASSID,1);
  if (bs < 1) dERROR(PETSC_ERR_ARG_OUTOFRANGE,"Block size must be at least 1, was %D",bs);
  dValidPointer(Y,3);

  err = VecGetLocalSize(X,&n);dCHK(err);
  err = VecGetBlockSize(X,&xbs);dCHK(err);
  err = VecCreate(((PetscObject)X)->comm,Y);dCHK(err);
  err = VecSetSizes(*Y,n/xbs*bs,PETSC_DETERMINE);dCHK(err);
  err = VecSetBlockSize(*Y,bs);dCHK(err);
  err = VecSetType(*Y,((PetscObject)X)->type_name);dCHK(err);
  dFunctionReturn(0);
}

dErr VecBlockView(Vec X,dViewer viewer)
{
  dErr err;
  dInt i,m,bs,rstart;
  dMPIInt rank;
  const dScalar *x;
  dBool ascii;

  dFunctionBegin;
  dValidHeader(X,VEC_CLASSID,1);
  if (!viewer) {err = PetscViewerASCIIGetStdout(((PetscObject)X)->comm,&viewer);dCHK(err);}
  dValidHeader(viewer,PETSC_VIEWER_CLASSID,2);
  err = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&ascii);dCHK(err);
  if (!ascii) dERROR(PETSC_ERR_SUP,"only ASCII");

  err = VecGetLocalSize(X,&m);dCHK(err);
  err = VecGetBlockSize(X,&bs);dCHK(err);
  err = VecGetOwnershipRange(X,&rstart,NULL);dCHK(err);
  err = MPI_Comm_rank(((PetscObject)X)->comm,&rank);dCHK(err);
  err = VecGetArrayRead(X,&x);dCHK(err);
  for (i=0; i<m; i+=bs) {
    switch (bs) {
      case 1:
      err = PetscViewerASCIISynchronizedPrintf(viewer,"[%d] %4D: %10G\n",rank,(rstart+i)/bs,PetscRealPart(x[i]));dCHK(err);
      break;
      case 2:
      err = PetscViewerASCIISynchronizedPrintf(viewer,"[%d] %4D: %10G %10G\n",rank,(rstart+i)/bs,PetscRealPart(x[i]),PetscRealPart(x[i+1]));dCHK(err);
      break;
      case 3:
      err = PetscViewerASCIISynchronizedPrintf(viewer,"[%d] %4D: %10G %10G %10G\n",rank,(rstart+i)/bs,PetscRealPart(x[i]),PetscRealPart(x[i+1]),PetscRealPart(x[i+2]));dCHK(err);
      break;
      default: dERROR(PETSC_ERR_SUP,"block size %D",bs);
    }
  }
  err = PetscViewerFlush(viewer);dCHK(err);
  err = VecRestoreArrayRead(X,&x);dCHK(err);
  dFunctionReturn(0);
}
