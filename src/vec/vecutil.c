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
