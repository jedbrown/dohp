#include <dohpvec.h>
#include <dohp.h>

#include "../src/vec/vec/impls/mpi/pvecimpl.h" /* To have access to Vec_MPI (.localrep) and VecCreate_MPI_Private */

static dErr VecStateSync_Private(Vec x,Vec y)
{
  dInt xstate,ystate;
  dErr err;

  dFunctionBegin;
  dValidHeader(x,VEC_CLASSID,1);
  dValidHeader(y,VEC_CLASSID,2);
  err = PetscObjectStateQuery((dObject)x,&xstate);dCHK(err);
  err = PetscObjectStateQuery((dObject)y,&ystate);dCHK(err);
  err = PetscObjectSetState((dObject)x,dMaxInt(xstate,ystate));dCHK(err);
  err = PetscObjectSetState((dObject)y,dMaxInt(xstate,ystate));dCHK(err);
  dFunctionReturn(0);
}


/** Get the closed form of a Dohp vector
*
* @note Dohp vectors are basically just MPI vectors, the only difference is that instead of a local form, we have a
* closed form.  We subvert .localrep to mean the closed form.
*
**/
dErr VecDohpGetClosure(Vec v,Vec *c)
{
  Vec_MPI *vmpi;
  dBool    isdohp;
  dErr     err;

  dFunctionBegin;
  dValidHeader(v,VEC_CLASSID,1);
  dValidPointer(c,2);
  err = PetscObjectTypeCompare((dObject)v,VECDOHP,&isdohp);dCHK(err);
  if (!isdohp) dERROR(PETSC_COMM_SELF,1,"Vector type %s does not have closure",((dObject)v)->type_name);
  vmpi = v->data;
  if (!vmpi->localrep) dERROR(PETSC_COMM_SELF,1,"Vector has no closure");
  *c = vmpi->localrep;
  err = VecStateSync_Private(v,*c);dCHK(err);
  err = PetscObjectReference((dObject)*c);dCHK(err);
  dFunctionReturn(0);
}

dErr VecDohpRestoreClosure(Vec v,Vec *c)
{
  dErr   err;
  dBool  isdohp;

  dFunctionBegin;
  dValidHeader(v,VEC_CLASSID,1);
  dValidPointer(c,2);
  err = PetscObjectTypeCompare((dObject)v,VECDOHP,&isdohp);dCHK(err);
  if (!isdohp) dERROR(PETSC_COMM_SELF,1,"Vector type %s does not have closure",((dObject)v)->type_name);
  if (*c != ((Vec_MPI*)v->data)->localrep) dERROR(PETSC_COMM_SELF,1,"attempting to restore incorrect closure");
  err = VecStateSync_Private(v,*c);dCHK(err);
  err = PetscObjectDereference((dObject)*c);dCHK(err);
  *c = NULL;
  dFunctionReturn(0);
}

dErr VecDohpZeroEntries(Vec v)
{
  dErr err;
  dBool  isdohp;
  Vec c;

  dFunctionBegin;
  dValidHeader(v,VEC_CLASSID,1);
  err = PetscObjectTypeCompare((dObject)v,VECDOHP,&isdohp);dCHK(err);
  if (!isdohp) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"Vector type %s",((dObject)v)->type_name);
  err = VecDohpGetClosure(v,&c);dCHK(err);
  err = VecZeroEntries(c);dCHK(err);
  err = VecDohpRestoreClosure(v,&c);dCHK(err);
  dFunctionReturn(0);
}

static dErr VecDuplicate_Dohp(Vec x,Vec *iny)
{
  Vec      y,xc,yc;
  Vec_MPI *ympi;
  dScalar *a;
  dErr     err;

  dFunctionBegin;
  dValidHeader(x,VEC_CLASSID,1);
  dValidPointer(iny,2);
  *iny = 0;
  err = VecDohpGetClosure(x,&xc);dCHK(err);
  err = VecDuplicate(xc,&yc);dCHK(err);
  err = VecDohpRestoreClosure(x,&xc);dCHK(err);

  /* The rest is mostly the same as VecDuplicate_MPI, but we can't call that because it allocates memory.
  * Unfortunately, this is fragile if the VecMPI implementation changes.  I think this part of PETSc is quite stable and
  * I will be sufficiently involved to notice changes here. Famous last words. */
  err = VecCreate(((dObject)x)->comm,&y);dCHK(err);

  err = PetscLayoutReference(x->map,&y->map);dCHK(err);

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
  y->map->bs   = x->map->bs;
  y->bstash.bs = x->bstash.bs;

  err = PetscObjectChangeTypeName((dObject)y,VECDOHP);dCHK(err);
  *iny = y;
  dFunctionReturn(0);
}

dErr VecCreateDohp(MPI_Comm comm,dInt bs,dInt n,dInt nc,dInt nghosts,const dInt ghosts[],Vec *v)
{
  Vec_MPI *vmpi;
  Vec      vc,vg;
  dScalar *a;
  dErr     err;

  dFunctionBegin;
  dValidPointer(v,7);
  *v = 0;
  err = VecCreateGhostBlock(comm,bs,nc*bs,PETSC_DECIDE,nghosts,ghosts,&vc);dCHK(err);
  err = VecGetArray(vc,&a);dCHK(err);
  err = VecCreateMPIWithArray(comm,bs,n*bs,PETSC_DECIDE,a,&vg);dCHK(err);
  err = VecRestoreArray(vc,&a);dCHK(err);
  err = VecSetBlockSize(vg,bs);dCHK(err);
  vmpi = vg->data;
  if (vmpi->localrep) dERROR(PETSC_COMM_SELF,1,"Vector has localrep, expected no localrep");
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

/** Create a cache for Dirichlet part of closure vector, and scatter from global closure to Dirichlet cache.

@arg[in] gvec Global vector
@arg[out] dcache New vector to hold the Dirichlet values
@arg[out] dscat Scatter from global closure to \a dcache

@note This could be local but it doesn't cost anything to make it global.
**/
dErr VecDohpCreateDirichletCache(Vec gvec,Vec *dcache,VecScatter *dscat)
{
  MPI_Comm comm;
  dErr     err;
  dBool    isdohp;
  IS       from;
  Vec      gc;
  dInt     n,nc,crstart;

  dFunctionBegin;
  dValidHeader(gvec,VEC_CLASSID,1);
  dValidPointer(dcache,2);
  dValidPointer(dscat,3);
  err = PetscObjectTypeCompare((PetscObject)gvec,VECDOHP,&isdohp);dCHK(err);
  if (!isdohp) dERROR(PETSC_COMM_SELF,PETSC_ERR_SUP,"Vec type %s",((PetscObject)gvec)->type_name);
  err = PetscObjectGetComm((PetscObject)gvec,&comm);dCHK(err);
  err = VecGetLocalSize(gvec,&n);dCHK(err);
  err = VecDohpGetClosure(gvec,&gc);dCHK(err);
  err = VecGetLocalSize(gc,&nc);dCHK(err);
  err = VecGetOwnershipRange(gc,&crstart,NULL);dCHK(err);
  err = VecCreateMPI(comm,nc-n,PETSC_DECIDE,dcache);dCHK(err);
  err = ISCreateStride(comm,nc-n,crstart+n,1,&from);dCHK(err);
  err = VecScatterCreate(gc,from,*dcache,NULL,dscat);dCHK(err);
  err = VecDohpRestoreClosure(gvec,&gc);dCHK(err);
  err = ISDestroy(&from);dCHK(err);
  /* \todo deal with rotations */
  dFunctionReturn(0);
}
