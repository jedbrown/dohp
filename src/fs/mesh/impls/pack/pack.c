#include "private/dmeshimpl.h"
#include <MBParallelConventions.h>
#include <stddef.h>

dErr dMeshCreate_Pack(dMesh mesh); /* The only exported function in this file */

/* This macro is in MBParallelComm.hpp and these should be the same */
#define MAX_SHARING_PROCS 64
#define MAX_NAME_LEN 128

struct SharedSet {
  dIInt ranks[MAX_SHARING_PROCS]; /* The canonical list of procs defines the set */
  dInt numranks;
  dMeshESH handle;
};

static int compareSharedSet(const void *avoid,const void *bvoid)
{
  const struct SharedSet *a = *(struct SharedSet *const*)avoid;
  const struct SharedSet *b = *(struct SharedSet *const*)bvoid;

  return memcmp(a->ranks,b->ranks,MAX_SHARING_PROCS*sizeof(dInt));
}

/**
* The \a buffer is sorted \a rank,\a set,\a tag with starts at \c kr[jr[ir[rank]+set]+tag].  The \a data is sorted \a
* tag,\a set,\a rank with starts at \c kt[jt[it[tag]+set]+rank].  \e Loading the \c Packer means populating \c data from
* the mesh interface, \e unloading the \c Packer means committing \c data to the mesh interface (with a reduction in
* each set if we are the owner).  \e Packing means moving from \c data order to \c buffer order, \e unpacking is the
* inverse.
*
* */
typedef struct {
  int myrank;
  dInt nr,ns,ne,nsr,nsre;       /* r=rank,s=set,e=entity, nsr=sum(snr)=sum(rns), nsre=sum(zipWith (*) snr sne) */
  dInt *snr,*rns,*sr;           /* snr[0..ns], rns[0..nr], sr[0..nsr],sr[is[s]+rr], rs[0..nsr],rs[ir[r]+ss] */
  dInt *ir,*is,*jd;             /* ir[0..nr], is[0..ns], ib[0..nr],ib[r], jd[0..nsr],jd[is[s]+rr] */
  dInt *rank,*irank;            /* rank[r] is the rank of proc with index r, irank[rank[r]]==r */

  dBool persistantcreated;
  MPI_Datatype *rtype;          /* Data spec for communication with each rank */
  PetscMPIInt *blen,*displ;     /* blen[ir[r]+ss] the block size and displacements for each set, used to create \a rtype */
  MPI_Request *sendreq,*recvreq;

  struct SharedSet *setarr;
  struct SharedSet **set;        /* Array of sets set[s].xxx must be initialized */
  dInt *estart,*sne;            /* Start and length of entity lists in \a ents */
  dMeshEH *ents;                /* All shared entities. */

  dMeshTag tag;
  dIInt ttype,tlen,tsize;       /* Types in ITAPS enum, length in values, size in bytes. */
  dInt maxtsize;                /* To decide whether the buffers need to be reallocated. */

  dInt datasize,minesize;
  void *data;                   /*   &data[jd[is[s]+r]*e*tsize] 0<=s<ns, 0<=r<snr[s], 0<=e<sne[s],  */
  void *mine;                   /*   &mine[e*tsize]             0<=e<ne */
} PackData;

/**
* \e Setup is required when the number or size of tags changes.
* 
*/
typedef struct {
  dInt nr,ns;
  char partition[dSTR_LEN];
  PackData *owned;
  PackData *unowned;
} Mesh_Pack;

/**
* There are two orderings, by rank (used for send/recieve) and by set (used during reduction operations).  We do not
* communicate unowned sets with any ranks except the owner.  Owned sets must be communicated both ways with all ranks
* that have copies.
* 
*/

typedef struct {
  PetscMPIInt numsets;
  char tagname[MAX_NAME_LEN];
} SharedHeader;

typedef struct {
  PetscMPIInt numranks;
  PetscMPIInt ranks[MAX_SHARING_PROCS];
} SharedSetHeader;

static const dInt iBase_SizeFromType[4] = {sizeof(int),sizeof(double),sizeof(void*),sizeof(char)};
static dBool CreatedMPITypes = false;
static MPI_Datatype MPITypeFromITAPS[4];
static MPI_Datatype dMPI_SHARED_HEADER,dMPI_SHARED_SET_HEADER,dMPI_ENT_HANDLE;

#if 0
static dErr dFindInt(dInt n,dInt a[],dInt v,dInt *i)
{
  dInt low,high;

  dFunctionBegin;
  low = 0; high = n;
  while (high-low > 1) {
    const dInt t = (low+high)/2;
    if (a[t] > v) high = t;
    else          low = t;
  }
  if (a[low] == v) *i = low;
  else dERROR(1,"Failed to find value %d in list [(%d,%d)..(%d,%d)].",v,0,a[0],n-1,a[n-1]);
  dFunctionReturn(0);
}
#endif

/**
* Creates all associations required for exchanging and reducing tags, but not anything specific to a particular tag
* 
* @param mesh 
* @param ns number of set handles
* @param sh array of set handles
* @param inpd new PackData, set up for the current entities
* 
* @return 
*/
static dErr MeshPackDataCreate(dMesh mesh,dInt ns,dMeshESH *sh,PackData **inpd)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  PackData *pd;
  iMesh_Instance mi = mesh->mi;
  dMeshTag pstatusTag,sprocTag,sprocsTag;
  PetscMPIInt rank,size;
  dInt *rcount,nsr,nr,ne,nsre;
  dErr err;

  dFunctionBegin;
  dValidPointer(inpd,4);
  *inpd = 0;
  err = dNew(PackData,&pd);dCHK(err);

  err = MPI_Comm_rank(comm,&rank);dCHK(err);
  err = MPI_Comm_size(comm,&size);dCHK(err);
  pd->myrank = rank;
  iMesh_getTagHandle(mi,PARALLEL_STATUS_TAG_NAME,&pstatusTag,&err,strlen(PARALLEL_STATUS_TAG_NAME));dICHK(mi,err);
  iMesh_getTagHandle(mi,PARALLEL_SHARED_PROC_TAG_NAME,&sprocTag,&err,strlen(PARALLEL_SHARED_PROC_TAG_NAME));dICHK(mi,err);
  iMesh_getTagHandle(mi,PARALLEL_SHARED_PROCS_TAG_NAME,&sprocsTag,&err,strlen(PARALLEL_SHARED_PROCS_TAG_NAME));dICHK(mi,err);

  err = dCalloc(size*sizeof(rcount[0]),&rcount);dCHK(err); /* number of sets each rank occurs in */

  /* Allocate everything that only depends on the number of sets */
  err = PetscMalloc6(ns,struct SharedSet,&pd->setarr, ns,struct SharedSet*,&pd->set, ns,dInt,&pd->is, ns,dInt,&pd->estart, ns,dInt,&pd->sne, ns,dInt,&pd->snr);dCHK(err);

  /* Create the canonical set representations so that we can sort them */
  for (dInt s=0; s<ns; s++) {
    char pstatus,*silly_status=&pstatus; /* silliness required due to the awkward iMesh interface */
    dIInt sproc,*sprocs,one_a=1,one_s,sprocs_a=sizeof(pd->setarr[0].ranks),sprocs_s;
    dInt r;

    sprocs = pd->setarr[s].ranks;
    iMesh_getEntSetData(mi,sh[s],pstatusTag,&silly_status,&one_a,&one_s,&err);dICHK(mi,err);
    if (!(pstatus & PSTATUS_SHARED)) dERROR(1,"shared packer contains an unshared set");

    /* make \c sprocs[] contain the list of copies, terminated by -1 */
    iMesh_getEntSetIntData(mi,sh[s],sprocTag,&sproc,&err);dICHK(mi,err);
    if (sproc == -1) {        /* There are multiple sharing procs */
      iMesh_getEntSetData(mi,sh[s],sprocsTag,(dIByte**)&sprocs,&sprocs_a,&sprocs_s,&err);dICHK(mi,err);
    } else {
      sprocs[0] = sproc; for (dInt i=1; i<MAX_SHARING_PROCS; i++) sprocs[i] = -1;
    }

    /* Find position 'r' to insert my rank */
    if (pstatus & PSTATUS_NOT_OWNED) {
      for (r=1; sprocs[r] < rank && sprocs[r] >= 0; r++) { /* Scan forward to my position */
        if (r == MAX_SHARING_PROCS-1) dERROR(1,"Too many sharing procs");
      }
    } else {                    /* I am the owner, insert at the beginning of the list */
      r = 0; 
    }

    {                           /* Insert my rank and determine the length in 'r' */
      dInt t0,t1;
      t0 = rank;
      while (t0 >= 0) {
        if (r == MAX_SHARING_PROCS) dERROR(1,"Too many sharing procs");
        t1 = sprocs[r];
        sprocs[r] = t0;
        t0 = t1;
        r++;
      }
    }
    pd->setarr[s].numranks = r;
    pd->setarr[s].handle = sh[s];
  }

  /* Sort the sets */
  for (dInt i=0; i<ns; i++) pd->set[i] = &pd->setarr[i];
  qsort(pd->set,ns,sizeof(pd->set[0]),compareSharedSet);

  nsr = 0;                      /* number of set-ranks seen */
  ne = 0;                       /* Number of entities seen (including duplicates if present in multiple sets) */
  nsre = 0;
  for (dInt s=0; s<ns; s++) {
    dIInt nents;
    dInt r;

    iMesh_getNumOfType(mi,pd->set[s]->handle,iBase_ALL_TYPES,&nents,&err);dICHK(mi,err);

    pd->sne[s] = nents;
    pd->estart[s] = ne;
    pd->is[s]  = nsr;           /* Offset of the first rank in set 's' */
    ne += nents;

    if (pd->set[s]->ranks[0] == rank) { /* We are the owner */
      r = pd->set[s]->numranks-1;       /* Number of remote ranks in set */
      pd->snr[s] = r;
      nsr += r; 
      nsre += r*nents;
      for (dInt j=0; j<r+1; j++) { /* count every rank */
        rcount[pd->set[s]->ranks[j]]++;
      }
    } else {                    /* We are not the owner */
      /* we only talk with the owner so we won't put the other sharing ranks into our communication structure, however
      * they are still part of the set specification */
      pd->snr[s] = 1;
      nsr += 1;
      nsre += nents;
      rcount[pd->set[s]->ranks[0]]++; /* count the owner */
      rcount[rank]++;                 /* count me, just for consistency check */
    }
  }

  if (rcount[rank] != ns) dERROR(1,"[%d] I should be in every set rcount=%d, ns=%d",rank,rcount[rank],ns);

  /* Count the number of distinct neighbor ranks */
  rcount[rank] = 0;             /* don't count me */
  nr = 0;
  for (dInt i=0; i<size; i++) {
    if (rcount[i]) nr++;
  }

  pd->ns = ns;
  pd->nsr = nsr;
  pd->nr = nr;
  pd->ne = ne;
  pd->nsre = nsre;

  /* Allocate for the rank mappings */
  err = PetscMalloc4(nr,dInt,&pd->ir,nr,dInt,&pd->rns,nr,dInt,&pd->rank,size,dInt,&pd->irank);dCHK(err);
  {
    dInt ir = 0,r = 0;
    for (dInt i=0; i<size; i++) {
      if (rcount[i]) {
        pd->rank[ir] = i;               /* initialize rank_index -> rank */
        pd->ir[ir] = r; r += rcount[i]; /* start of each rank's sets */
        pd->rns[ir] = rcount[i];        /* number of sets for this rank */
        pd->irank[i] = ir;
        ir++;
      } else {
        pd->irank[i] = -1;
      }
    }
    if (r != pd->nsr) dERROR(1,"counts do not match");
  }

  /* Set up everything which does not depend on tag size */
  err = PetscMalloc5(nsr,MPI_Datatype,&pd->rtype,nsr,PetscMPIInt,&pd->blen,nsr,PetscMPIInt,&pd->displ,nr,MPI_Request,&pd->recvreq,nr,MPI_Request,&pd->sendreq);dCHK(err);
  err = PetscMalloc3(nsr,dInt,&pd->sr,nsr,dInt,&pd->jd,ne,dMeshEH,&pd->ents);dCHK(err);
  /* We do not initialize the \e values in \a jb and \a jd since these are dependent on the tag size */
  err = PetscMemzero(rcount,size*sizeof(rcount[0]));dCHK(err); /* Count the number of times we hit each rank */
  for (dInt s=0; s<ns; s++) {
    dInt i,r,rr,arrsize,arralloc,start,end;
    dMeshEH *arr;
    if (pd->set[s]->ranks[0] == rank) { /* I am the owner, we need all the subsets */
      start = 1; end = 1+pd->snr[s];
    } else {                    /* Just look at the owner */
      if (pd->snr[s] != 1) dERROR(1,"We are not the owner, but we still have multiple ranks in this set");
      start = 0; end = 1;
    }
    for (i=start,r=0; i != end; i++,r++) {
      rr = pd->set[s]->ranks[i];  /* the \e real rank */
      pd->sr[pd->is[s]+r] = pd->irank[rr];
      if (!(0 <= pd->irank[rr] && pd->irank[rr] < pd->nr)) dERROR(1,"[%d] invalid real rank %d",rank,rr);
      if (!(0 <= pd->ir[pd->irank[rr]] && pd->ir[pd->irank[rr]] < pd->nsr)) dERROR(1,"pd->ir invalid");
    }
    arr = &pd->ents[pd->estart[s]];
    arralloc = pd->sne[s];
    iMesh_getEntities(mi,pd->set[s]->handle,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,&arr,&arralloc,&arrsize,&err);dICHK(mi,err);
    if (arrsize != pd->sne[s]) dERROR(1,"incorrect sizes");
  } 
  err = dFree(rcount);dCHK(err);
  *inpd = pd;
  dFunctionReturn(0);
}

static dErr MeshPackDataPrepare(dMesh mesh,PackData *pd,dMeshTag tag)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  iMesh_Instance mi = mesh->mi;
  dIInt tsize,tlen,ttype;
  dInt pos,rr,*rcount;
  dErr err;

  dFunctionBegin;
  dValidHeader(mesh,dMESH_COOKIE,1);
  dValidPointer(pd,2);
  if (pd->tag == tag) dFunctionReturn(0);  /* If the tag has not changed, there is nothing to do */
  iMesh_getTagSizeValues(mi,tag,&tlen,&err);dICHK(mi,err);
  iMesh_getTagType(mi,tag,&ttype,&err);dICHK(mi,err);
  tsize = tlen * iBase_SizeFromType[ttype];
  if (pd->maxtsize < tsize) {  /* Reallocate the buffers */
    err = PetscFree2(pd->data,pd->mine);dCHK(err); /* the first time around, these will be NULL which is safe with PetscFree2 */
    pd->maxtsize   = tsize;
    pd->datasize   = pd->nsre * tsize;
    pd->minesize   = pd->ne * tsize;
    err = PetscMalloc2(pd->datasize,char,&pd->data,pd->minesize,char,&pd->mine);dCHK(err);
  }
  pd->tag = tag;
  pd->tsize = tsize;
  pd->tlen = tlen;
  pd->ttype = ttype;

  if (pd->persistantcreated) {
    for (dInt r=0; r<pd->nr; r++) {
      err = MPI_Type_free(&pd->rtype[r]);dCHK(err);
      err = MPI_Request_free(&pd->sendreq[r]);dCHK(err);
      err = MPI_Request_free(&pd->recvreq[r]);dCHK(err);
    }
  }

  /* Initialize the \a jd array, it is in set-rank order (the start of each rank within a set) */
  err = dCalloc(pd->nr*sizeof(rcount[0]),&rcount);dCHK(err);
  pos = 0;
  for (dInt s=0; s<pd->ns; s++) {
    for (dInt r=0; r<pd->snr[s]; r++) {
      pd->jd[pd->is[s]+r] = pos; /* start in the data array */
      rr = pd->sr[pd->is[s]+r];  /* The rank index */
      {
        int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        printf("[%d] s=%d/%d r=%d/%d rr=%d/%d(%d) sne=%d/%d\n",rank,s,pd->ns,r,pd->snr[s],rr,pd->nr,pd->rank[rr],pd->sne[s],pd->ne); fflush(stdout);
      }
      pd->displ[pd->ir[rr]+rcount[rr]] = pos;
      pd->blen[pd->ir[rr]+rcount[rr]] = pd->sne[s] * tsize;
      rcount[rr]++;
      pos += pd->sne[s] * tsize;
    }
  }

  for (dInt r=0; r<pd->nr; r++) {
    if (rcount[r] != pd->rns[r]) dERROR(1,"rank-set counts do not agree, rank %d(%d) %d %d",r,pd->rank[r],pd->rns[r],rcount[r]);
  }

  for (dInt r=0; r<pd->nr; r++) {
    err = MPI_Type_indexed(rcount[r],pd->blen+pd->ir[r],pd->displ+pd->ir[r],MPI_BYTE,&pd->rtype[r]);dCHK(err);
    err = MPI_Type_commit(&pd->rtype[r]);dCHK(err);
    err = MPI_Send_init(pd->data,1,pd->rtype[r],pd->rank[r],((dObject)mesh)->tag,comm,&pd->sendreq[r]);dCHK(err);
    err = MPI_Recv_init(pd->data,1,pd->rtype[r],pd->rank[r],((dObject)mesh)->tag,comm,&pd->recvreq[r]);dCHK(err);
  }
  pd->persistantcreated = true;
  err = dFree(rcount);dCHK(err);
  dFunctionReturn(0);
}

/** 
* Load the tag values from the mesh into the \a data array.
*
* @pre dMeshPackDataPrepare() has been called with the current tag

* @param mesh 
* @param pd 
* 
* @return 
*/
static dErr MeshPackDataLoad(dMesh mesh,PackData *pd)
{
  iMesh_Instance mi = mesh->mi;
  dErr err;

  dFunctionBegin;
  /* Load \a mine from the mesh */
  {
    dIInt used;
    iMesh_getArrData(mi,pd->ents,pd->ne,pd->tag,(dIByte**)&pd->mine,&pd->minesize,&used,&err);dICHK(mi,err);
  }

  /* Copy \a mine into \a data */
  for (dInt s=0; s<pd->ns; s++) {
    for (dInt r=0; r<pd->snr[s]; r++) {
      err = dMemcpy((char*)pd->data+pd->jd[pd->is[s]+r],(char*)pd->mine+pd->estart[s]*pd->tsize,pd->sne[s]*pd->tsize);dCHK(err);
    }
  }
  dFunctionReturn(0);
}
  
static dErr MeshPackDataUnload(dMesh mesh,PackData *pd)
{
  iMesh_Instance mi = mesh->mi;
  dErr err;

  dFunctionBegin;
#if defined(DOHP_USE_DEBUG)
  for (int i=0; i<pd->minesize/(int)sizeof(int); i++) ((int*)pd->mine)[i] = -1; /* should all be overwritten */
#endif
  /* Put \a data into \a mine, copy over existing data if there are multiple sets */
  for (dInt s=0; s<pd->ns; s++) {
    for (dInt r=0; r<pd->snr[s]; r++) {
      err = dMemcpy((char*)pd->mine+pd->estart[s]*pd->tsize,(char*)pd->data+pd->jd[pd->is[s]+r],pd->sne[s]*pd->tsize);dCHK(err);
    }
  }
  iMesh_setArrData(mi,pd->ents,pd->ne,pd->tag,pd->mine,pd->ne*pd->tsize,&err);dICHK(mi,err);
  dFunctionReturn(0);
}

static dErr PackDataDestroy(PackData *pd)
{
  dErr err;

  dFunctionBegin;
  //err = PetscFree6(pd->set,pd->is,pd->estart,pd->sne,pd->snr,pd->spacksize);dCHK(err);
  err = PetscFree4(pd->ir,pd->rns,pd->rank,pd->irank);dCHK(err);
  //err = PetscFree3(pd->sr,pd->rs,pd->jd);dCHK(err);
  //err = PetscFree5(pd->ib,pd->ibsize,pd->ibused,pd->mpireq,pd->ents);dCHK(err);
  err = PetscFree2(pd->data,pd->mine);dCHK(err);
  dFunctionReturn(0);
}

static dErr CreateMPITypes(void)
{
  MPI_Datatype typeSH[] = {MPI_INT,MPI_CHAR};
  MPI_Aint dispSH[] = {offsetof(SharedHeader,numsets),offsetof(SharedHeader,tagname)};
  int sizeSH[] = {1,MAX_NAME_LEN};
  MPI_Datatype typeSSH[] = {MPI_INT,MPI_INT};
  MPI_Aint dispSSH[] = {offsetof(SharedSetHeader,numranks),offsetof(SharedSetHeader,ranks)};
  int sizeSSH[] = {1,MAX_SHARING_PROCS};
  dErr err;

  dFunctionBegin;
  if (CreatedMPITypes) dFunctionReturn(0);
  err = MPI_Type_create_struct(2,sizeSH,dispSH,typeSH,&dMPI_SHARED_HEADER);dCHK(err);
  err = MPI_Type_commit(&dMPI_SHARED_HEADER);dCHK(err);
  err = MPI_Type_create_struct(2,sizeSSH,dispSSH,typeSSH,&dMPI_SHARED_SET_HEADER);dCHK(err);
  err = MPI_Type_commit(&dMPI_SHARED_SET_HEADER);dCHK(err);
  err = MPI_Type_dup(MPI_UNSIGNED_LONG,&dMPI_ENT_HANDLE);dCHK(err);
  MPITypeFromITAPS[0] = MPI_INT;
  MPITypeFromITAPS[1] = MPI_DOUBLE;
  MPITypeFromITAPS[2] = dMPI_ENT_HANDLE;
  MPITypeFromITAPS[3] = MPI_BYTE;
  CreatedMPITypes = true;
  dFunctionReturn(0);
}


static dErr PackView(dMesh mesh,Mesh_Pack *pack)
{
  iMesh_Instance mi = mesh->mi;
  int err;
  int rank,nents;
  int osi = pack->owned->ns,usi = pack->unowned->ns;
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscSynchronizedPrintf(MPI_COMM_WORLD,"[%d] owned sets = %d   unowned sets = %d\n",rank,osi,usi);
  for (int i=0; i<osi; i++) {
    int ntype[4];
    iMesh_getNumOfType(mi,pack->owned->set[i]->handle,iBase_ALL_TYPES,&nents,&err);
    for (int j=0; j<4; j++) {
      iMesh_getNumOfType(mi,pack->owned->set[i]->handle,j,&ntype[j],&err);
    }
    PetscSynchronizedPrintf(MPI_COMM_WORLD,"[%d] owned[%d] ents=%d(%d,%d,%d,%d)\n",rank,i,nents,ntype[0],ntype[1],ntype[2],ntype[3]);
  }
  for (int i=0; i<usi; i++) {
    int ntype[4];
    iMesh_getNumOfType(mi,pack->unowned->set[i]->handle,iBase_ALL_TYPES,&nents,&err);
    for (int j=0; j<4; j++) {
      iMesh_getNumOfType(mi,pack->unowned->set[i]->handle,j,&ntype[j],&err);
    }
    PetscSynchronizedPrintf(MPI_COMM_WORLD,"[%d] unowned[%d] ents=%d(%d,%d,%d,%d) owner=%d members=[%d %d %d %d %d ...]\n",
                            rank,i,nents,ntype[0],ntype[1],ntype[2],ntype[3],
                            pack->unowned->rank[pack->unowned->sr[i]],
                            pack->unowned->set[i]->ranks[0],
                            pack->unowned->set[i]->ranks[1],
                            pack->unowned->set[i]->ranks[2],
                            pack->unowned->set[i]->ranks[3],
                            pack->unowned->set[i]->ranks[4]
                            );
    {
      dMeshEH *ents;
      int ea=0,es,*data,da=0,ds;
      iMesh_getEntities(mi,pack->unowned->set[i]->handle,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,&ents,&ea,&es,&err);dICHK(mi,err);
      iMesh_getIntArrData(mi,ents,es,pack->unowned->tag,&data,&da,&ds,&err);dICHK(mi,err);
      if (es != ds) dERROR(1,"iMesh inconsistency");
      for (int j=0; j<ds; j++) {
        PetscSynchronizedPrintf(MPI_COMM_WORLD,"[%d]   owner[%d] = %d\n",rank,j,data[j]);
      }
      free(ents); free(data);
    }
  }
  PetscSynchronizedFlush(MPI_COMM_WORLD);
  return 0;
}

/**
* Sends the requested tag from owner to shared.  It sends one message per neighbor.  The tags are assumed to be on all
* interface entities, at least on the owner.  The tag (with same name) must exist on the non-owning proc, but it may not
* yet be attached to the interface entities.
*
* @param mesh the mesh, should have neighboring entities set up
* @param tag tag handle
*
* @return
*/
static dErr dMeshTagBcast_Pack(dMesh mesh,dMeshTag tag)
{
  Mesh_Pack *pack = mesh->data;
  dBool haveowned,haveunowned;
  dErr err;

  dFunctionBegin;
  MPI_Barrier(MPI_COMM_WORLD);
  err = PetscPrintf(PETSC_COMM_WORLD,"---------------- beginning prepare\n");dCHK(err);
  MPI_Barrier(MPI_COMM_WORLD);
  err = MeshPackDataPrepare(mesh,pack->owned,tag);dCHK(err);
  MPI_Barrier(MPI_COMM_WORLD);
  err = PetscPrintf(PETSC_COMM_WORLD,"---------------- done preparing owned\n");dCHK(err);
  MPI_Barrier(MPI_COMM_WORLD);
  err = MeshPackDataPrepare(mesh,pack->unowned,tag);dCHK(err);
  MPI_Barrier(MPI_COMM_WORLD);
  err = PetscPrintf(PETSC_COMM_WORLD,"---------------- done preparing unowned\n");dCHK(err);
  MPI_Barrier(MPI_COMM_WORLD);

  haveowned = !!pack->owned->nr;
  haveunowned = !!pack->unowned->nr;
  err = MeshPackDataLoad(mesh,pack->owned);dCHK(err);
  if (haveunowned) { err = MPI_Startall(pack->unowned->nr,pack->unowned->recvreq);dCHK(err); }
  if (haveowned) { err = MPI_Startall(pack->owned->nr,pack->owned->sendreq);dCHK(err); }
  err = MPI_Waitall(pack->unowned->nr,pack->unowned->recvreq,MPI_STATUSES_IGNORE);dCHK(err); /* Wait on the receives */
  err = MeshPackDataUnload(mesh,pack->unowned);dCHK(err);
  err = MPI_Waitall(pack->owned->nr,pack->owned->sendreq,MPI_STATUSES_IGNORE);dCHK(err); /* Make sure the sends all completed */

  MPI_Barrier(MPI_COMM_WORLD);
  err = PackView(mesh,pack);dCHK(err);
  MPI_Barrier(MPI_COMM_WORLD);
  dFunctionReturn(0);
}

static dErr dMeshLoad_Pack(dMesh mesh)
{
  iMesh_Instance mi = mesh->mi;
  Mesh_Pack *pack = mesh->data;
  char options[dSTR_LEN];
  dMeshTag pstatusTag;
  dIInt nallset;
  dMeshESH *allset;
  size_t fnamelen;
  dErr err;

  dFunctionBegin;
  err = PetscStrlen(mesh->infile,&fnamelen);dCHK(err);
  err = PetscSNPrintf(options,sizeof(options),"PARALLEL=BCAST_DELETE;PARTITION=%s;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;CPUTIME;%s",pack->partition,mesh->inoptions);dCHK(err);
  iMesh_load(mi,0,mesh->infile,options,&err,(int)fnamelen,(int)sizeof(options));dICHK(mi,err);

  /* Get the interface sets */
  iMesh_getTagHandle(mi,PARALLEL_STATUS_TAG_NAME,&pstatusTag,&err,strlen(PARALLEL_STATUS_TAG_NAME));dICHK(mi,err);
  iMesh_getNumEntSets(mi,mesh->root,1,&nallset,&err);dICHK(mi,err);
  {
    dInt alloc = nallset;
    err = dMalloc(nallset*sizeof(dMeshESH),&allset);dCHK(err);
    iMesh_getEntSets(mi,mesh->root,1,&allset,&alloc,&nallset,&err);dICHK(mi,err);
  }

  /* Partition based on ownership, create owned and unowned PackData */
  {
    dMeshESH *oset,*uset;       /* sets in each partition */
    dInt osi=0,usi=0;           /* number of sets */
    err = PetscMalloc2(nallset,dMeshESH,&oset,nallset,dMeshESH,&uset);

    for (dInt i=0; i<nallset; i++) {
      char pstatus,*stupid = &pstatus;
      int one_a=1,one_s;

      iMesh_getEntSetData(mi,allset[i],pstatusTag,&stupid,&one_a,&one_s,&err);dICHK(mi,err);
      if (!(pstatus & PSTATUS_SHARED)) continue;
      if (pstatus & PSTATUS_NOT_OWNED) {
        uset[usi++] = allset[i];
      } else {
        oset[osi++] = allset[i];
      }
    }
    err = MeshPackDataCreate(mesh,osi,oset,&pack->owned);dCHK(err);
    err = MeshPackDataCreate(mesh,usi,uset,&pack->unowned);dCHK(err);
    err = PetscFree2(oset,uset);dCHK(err);
  }
  err = dFree(allset);dCHK(err);
  { int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    printf("[%d] %s finished, %p %p\n",rank,__func__,(void*)pack->owned,(void*)pack->unowned);
  }
  dFunctionReturn(0);
}

static dErr dMeshDestroy_Pack(dMesh mesh)
{
  Mesh_Pack *pack = mesh->data;
  dErr err;

  dFunctionBegin;
  err = PackDataDestroy(pack->owned);dCHK(err);
  err = PackDataDestroy(pack->unowned);dCHK(err);
  err = PetscFree(pack->owned);dCHK(err);
  err = PetscFree(pack->unowned);dCHK(err);
  err = dFree(mesh->data);dCHK(err); /* Like this so the pointer is zeroed */
  dFunctionReturn(0);
}

static dErr dMeshSetFromOptions_Pack(dMesh mesh)
{
  static const char defaultPartition[] = "PARALLEL_PARTITION";
  Mesh_Pack *pack = mesh->data;
  dBool flg;
  dErr err;

  dFunctionBegin;
  if (!pack->partition[0]) {err = PetscStrcpy(pack->partition,defaultPartition);dCHK(err);}
  err = PetscOptionsString("-dmesh_partition","Name of partition tag","dMeshSetInFile",defaultPartition,pack->partition,sizeof(pack->partition),&flg);dCHK(err);
  dFunctionReturn(0);
}

dErr dMeshCreate_Pack(dMesh mesh)
{
  Mesh_Pack *pack;
  dErr err;

  dFunctionBegin;
  err = CreateMPITypes();dCHK(err);
  err = dNewLog(mesh,Mesh_Pack,&pack);dCHK(err);
  mesh->data = (void*)pack;
  err = PetscNew(PackData,&pack->owned);dCHK(err);
  err = PetscNew(PackData,&pack->unowned);dCHK(err);

  iMesh_newMesh("PARALLEL",&mesh->mi,&err,(int)strlen("PARALLEL"));dICHK(mesh->mi,err);

  mesh->ops->view           = 0;
  mesh->ops->destroy        = dMeshDestroy_Pack;
  mesh->ops->tagbcast       = dMeshTagBcast_Pack;
  mesh->ops->setfromoptions = dMeshSetFromOptions_Pack;
  mesh->ops->load           = dMeshLoad_Pack;
  dFunctionReturn(0);
}
