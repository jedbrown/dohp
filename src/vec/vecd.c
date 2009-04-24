#include "dohpvec.h"
#include "../src/vec/vec/impls/mpi/pvecimpl.h" /* To have access to Vec_MPI (.localrep) and VecCreate_MPI_Private */

/** Get the closed form of a Dohp vector
*
* @note Dohp vectors are basically just MPI vectors, the only difference is that instead of a local form, we have a
* closed form.  We subvert .localrep to mean the closed form.
*
**/
dErr VecDohpGetClosure(Vec v,Vec *c)
{
  Vec_MPI *vmpi;
  dTruth   isdohp;
  dErr     err;

  dFunctionBegin;
  err = PetscTypeCompare((dObject)v,VECDOHP,&isdohp);dCHK(err);
  if (!isdohp) dERROR(1,"Vector type %s does not have closure",((dObject)v)->type_name);
  vmpi = v->data;
  if (!vmpi->localrep) dERROR(1,"Vector has no closure");
  *c = vmpi->localrep;
  err = PetscObjectReference((dObject)*c);dCHK(err);
  dFunctionReturn(0);
}

dErr VecDohpRestoreClosure(Vec v,Vec *c)
{
  dErr err;

  dFunctionBegin;
  if (*c != ((Vec_MPI*)v->data)->localrep) dERROR(1,"attempting to restore incorrect closure");
  err = PetscObjectDereference((dObject)*c);dCHK(err);
  dFunctionReturn(0);
}

static dErr VecDuplicate_Dohp(Vec x,Vec *iny)
{
  Vec      y,xc,yc;
  Vec_MPI *ympi;
  dScalar *a;
  dErr     err;

  dFunctionBegin;
  dValidHeader(x,VEC_COOKIE,1);
  dValidPointer(y,2);
  *iny = 0;
  err = VecDohpGetClosure(x,&xc);dCHK(err);
  err = VecDuplicate(xc,&yc);dCHK(err);
  err = VecDohpRestoreClosure(x,&xc);dCHK(err);

  /* The rest is mostly the same as VecDuplicate_MPI, but we can't call that because it allocates memory.
  * Unfortunately, this is fragile if the VecMPI implementation changes.  I think this part of PETSc is quite stable and
  * I will be sufficiently involved to notice changes here. Famous last words. */
  err = VecCreate(((dObject)x)->comm,&y);dCHK(err);

  err = PetscMapDestroy(y->map);dCHK(err);
  y->map = x->map;
  y->map->refcnt++;

  err = VecGetArray(yc,&a);dCHK(err);
  err = VecCreate_MPI_Private(y,PETSC_FALSE,0,a);dCHK(err);
  err = VecRestoreArray(yc,&a);dCHK(err);
  ympi = y->data;
  err = dMemcpy(y->ops,x->ops,sizeof(struct _VecOps));dCHK(err);

  ympi->localrep = yc;             /* subverting .localrep to mean closed form */

  y->stash.donotstash   = x->stash.donotstash;
  y->stash.ignorenegidx = x->stash.ignorenegidx;

  err = PetscOListDuplicate(((dObject)x)->olist,&((dObject)y)->olist);dCHK(err);
  err = PetscFListDuplicate(((dObject)x)->qlist,&((dObject)y)->qlist);dCHK(err);
  if (x->mapping) {
    err = PetscObjectReference((dObject)x->mapping);dCHK(err);
    y->mapping = x->mapping;
  }
  if (x->bmapping) {
    err = PetscObjectReference((dObject)x->bmapping);dCHK(err);
    y->bmapping = x->bmapping;
  }
  y->map->bs   = x->map->bs;
  y->bstash.bs = x->bstash.bs;

  err = PetscObjectChangeTypeName((dObject)y,VECDOHP);dCHK(err);
  *iny = y;
  dFunctionReturn(0);
}

#if 0
static dErr VecDestroy_Dohp(Vec x)
{
  dErr err;

  dFunctionBegin;
  err = PetscObjectChangeTypeName((dObject)x,VECMPI);
#endif

dErr VecCreateDohp(MPI_Comm comm,dInt bs,dInt n,dInt nc,dInt nghosts,const dInt ghosts[],Vec *v)
{
  Vec_MPI *vmpi;
  dInt    *sghosts;
  Vec      vc,vg;
  dScalar *a;
  dErr     err;

  dFunctionBegin;
  dValidPointer(v,7);
  *v = 0;
  if (bs > 1) {
    err = dMallocA(nghosts,&sghosts);dCHK(err);
    for (dInt i=0; i<nghosts; i++) sghosts[i] = ghosts[i]*bs; /* Index ghosts by scalar offset instead of blocks */
  } else {
    sghosts = 0;
  }
  err = VecCreateGhostBlock(comm,bs,nc*bs,PETSC_DECIDE,nghosts,sghosts?sghosts:ghosts,&vc);dCHK(err);
  err = dFree(sghosts);dCHK(err);
  err = VecGetArray(vc,&a);dCHK(err);
  err = VecCreateMPIWithArray(comm,n*bs,PETSC_DECIDE,a,&vg);dCHK(err);
  err = VecRestoreArray(vc,&a);dCHK(err);
  err = VecSetBlockSize(vg,bs);dCHK(err);
  vmpi = vg->data;
  if (vmpi->localrep) dERROR(1,"Vector has localrep, expected no localrep");
  vmpi->localrep = vc;          /* subvert this field to mean closed rep */
  /* Since we subvect .localrep, VecDestroy_MPI will automatically destroy the closed form */
  vg->ops->duplicate = VecDuplicate_Dohp;
  //vg->ops->destroy   = VecDestroy_Dohp;
  /* It might be useful to set the (block) LocalToGlobal mapping here, but in the use case I have in mind, the user is
  * always working with the closed form anyway (in function evaluation).  The \e matrix does need a customized
  * LocalToGlobal mapping.
  */
  err = PetscObjectChangeTypeName((dObject)vg,VECDOHP);dCHK(err);
  *v = vg;
  dFunctionReturn(0);
}
