#include "private/dmeshimpl.h"
#include "uthash.h"
#include <MBParallelConventions.h>

dErr dMeshCreate_Pack(dMesh mesh); /* The only exported function in this file */

/* This macro is in MBParallelComm.hpp and these should be the same */
#define MAX_SHARING_PROCS 64
#define MAX_NAME_LEN 128

struct SharedSet {
  dInt ranks[MAX_SHARING_PROCS]; /* The canonical list of procs defines the set */
  dInt setnumber;
  UT_hash_handle hh;
};

/**
* The \a buffer is sorted \a rank,\a set,\a tag with starts at \c kr[jr[ir[rank]+set]+tag].  The \a data is sorted \a
* tag,\a set,\a rank with starts at \c kt[jt[it[tag]+set]+rank].  \e Loading the \c Packer means populating \c data from
* the mesh interface, \e unloading the \c Packer means committing \c data to the mesh interface (with a reduction in
* each set if we are the owner).  \e Packing means moving from \c data order to \c buffer order, \e unpacking is the
* inverse.
*
* */
typedef struct {
  dInt nr,ns,ne,nsr,nsre;       /* r=rank,s=set,e=entity, nsr=sum(snr)=sum(rns), nsre=sum(zipWith (*) snr sne) */
  dInt *snr,*rns,*sr,*rs;       /* snr[0..ns], rns[0..nr], sr[0..nsr],sr[is[s]+rr], rs[0..nsr],rs[ir[r]+ss] */
  dInt *ir,*is,*jd;             /* ir[0..nr], is[0..ns], ib[0..nr],ib[r], jd[0..nsr],jd[is[s]+rr] */
  dInt *ib,*ibsize,*ibused;     /* Buffer index, size available, size used */
  dInt *rank,*irank;            /* rank[r] is the rank of proc with index r, irank[rank[r]]==r */

  PetscMPIInt mpitag;
  MPI_Request *mpireq;
  struct SharedSet *set;        /* Array of sets set[s].xxx must be initialized */
  struct SharedSet *sethash;    /* Used to lookup the set a message from a proc corresponds to */
  dMeshESH *sh;                 /* Set handles, indexed by set number */
  dInt *estart,*sne;            /* Start and length of entity lists in \a ents */
  dInt *spacksize;
  dMeshEH *ents;                /* All shared entities. */

  dMeshTag tag;
  dIInt ttype,tlen,tsize;       /* Types in ITAPS enum, length in values, size in bytes. */
  dInt maxtsize;                /* To decide whether the buffers need to be reallocated. */

  dInt buffersize,datasize,minesize;
  void *buffer;                 /* &buffer[ib[r]] 0<=r<nr */
  void *data;                   /*   &data[jd[is[s]+r]*e*tsize] 0<=s<ns, 0<=r<snr[s], 0<=e<sne[s],  */
  void *mine;                   /*   &mine[e*tsize]             0<=e<ne */
} PackData;

/**
* \e Setup is required when the number or size of tags changes.
* 
*/
typedef struct {
  dInt nr,ns;
  char *partition;
  PackData *owned;
  PackData *unowned;
} Mesh_Pack;

static dErr MeshPackDataSetUp(dMesh mesh,PackData *pd);

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
  else dERROR(1,"Not found.");
  dFunctionReturn(0);
}

/**
* Creates all the other associations required for exchanging and reducing tags.
* 
* @pre The set handles array (.sh) with size (.ns) has been initialized
* 
* @param mesh 
* @param pd 
* @param owner 
* 
* @return 
*/
static dErr MeshPackDataSetUp(dMesh mesh,PackData *pd)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  iMesh_Instance mi = mesh->mi;
  const dInt ns = pd->ns;
  dMeshTag pstatusTag,sprocTag,sprocsTag;
  dInt rank,size,*rcount,nsr,nr,ne,nsre;
  dErr err;

  dFunctionBegin;
  err = MPI_Comm_rank(comm,&rank);dCHK(err);
  err = MPI_Comm_size(comm,&size);dCHK(err);
  iMesh_getTagHandle(mi,PARALLEL_STATUS_TAG_NAME,&pstatusTag,&err,strlen(PARALLEL_STATUS_TAG_NAME));dICHK(mi,err);
  iMesh_getTagHandle(mi,PARALLEL_SHARED_PROC_TAG_NAME,&sprocTag,&err,strlen(PARALLEL_SHARED_PROC_TAG_NAME));dICHK(mi,err);
  iMesh_getTagHandle(mi,PARALLEL_SHARED_PROCS_TAG_NAME,&sprocsTag,&err,strlen(PARALLEL_SHARED_PROCS_TAG_NAME));dICHK(mi,err);

  err = dCalloc(size*sizeof(rcount[0]),&rcount);dCHK(err); /* number of sets each rank occurs in */

  /* Allocate everything that only depends on the number of sets */
  err = PetscMalloc6(ns,struct SharedSet,&pd->set, ns,dInt,&pd->is, ns,dInt,&pd->estart, ns,dInt,&pd->sne, ns,dInt,&pd->snr, ns,dInt,&pd->spacksize);dCHK(err);

  nsr = 0;                      /* number of set-ranks seen */
  ne = 0;                       /* Number of entities seen (including duplicates if present in multiple sets) */
  nsre = 0;
  for (dInt s=0; s<ns; s++) {
    char pstatus,*silly_status=&pstatus;
    dIInt sproc,sprocs[MAX_SHARING_PROCS],*silly_sprocs=sprocs,one_a=1,one_s,sprocs_a=sizeof(sprocs),sprocs_s,nents;
    dInt r;

    iMesh_getEntSetData(mi,pd->sh[s],pstatusTag,&silly_status,&one_a,&one_s,&err);dICHK(mi,err);
    if (!(pstatus & PSTATUS_SHARED)) dERROR(1,"shared packer contains an unshared set");
      
    /* make \c sprocs[] contain the list of copies, terminated by -1 */
    iMesh_getEntSetIntData(mi,pd->sh[s],sprocTag,&sproc,&err);dICHK(mi,err);
    if (sproc == -1) {        /* There are multiple sharing procs */
      iMesh_getEntSetData(mi,pd->sh[s],sprocsTag,(dIByte**)&silly_sprocs,&sprocs_a,&sprocs_s,&err);dICHK(mi,err);
    } else {
      sprocs[0] = sproc; sprocs[1] = -1;
    }

    /* Find my position in the list */
    if (pstatus & PSTATUS_NOT_OWNED) {                        /* Insert my rank in the list */
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
    iMesh_getNumOfType(mi,pd->sh[s],iBase_ALL_TYPES,&nents,&err);dICHK(mi,err);

    pd->sne[s] = nents;
    pd->estart[s] = ne;
    pd->is[s]  = nsr;           /* Offset of the first rank in set 's' */
    pd->snr[s] = r-1;           /* Number of remote ranks in set */
    nsr += r-1;                 /* Count all but me */
    ne += nents;
    nsre += (r-1)*nents;

    pd->set[s].setnumber = s;
    err = dMemcpy(pd->set[s].ranks,sprocs,MAX_SHARING_PROCS*sizeof(sprocs[0]));dCHK(err);
    //HASH_ADD(hh,pd->sethash,ranks,sizeof(pd->set[0].ranks),&pd->set[s]);

    for (dInt j=0; j<r; j++) {
      rcount[sprocs[j]]++; /* Update set count for each rank in set */
    }
  }

  /* I should be in every set */
  if (rcount[rank] != ns) dERROR(1,"rank counts do not match");

  /* Count the number of distinct neighbor ranks */
  rcount[rank] = 0;             /* don't count me */
  nr = 0;
  for (dInt i=0; i<size; i++) {
    if (rcount[i]) nr++;
  }

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
  err = dFree(rcount);dCHK(err);

  /* Set up everything which does not depend on tag size */
  err = PetscMalloc3(nsr,dInt,&pd->sr,nsr,dInt,&pd->rs,nsr,dInt,&pd->jd);dCHK(err);
  err = PetscMalloc5(nr,dInt,&pd->ib,nr,dInt,&pd->ibsize,nr,dInt,&pd->ibused,nr,MPI_Request,&pd->mpireq,ne,dMeshEH,&pd->ents);dCHK(err);
  /* We do not initialize the \e values in \a jb and \a jd since these are dependent on the tag size */
  for (dInt s=0; s<ns; s++) {
    dInt i,r,rr,arrsize,arralloc;
    dMeshEH *arr;
    for (i=0,r=0; i<pd->snr[s]+1; i++,r++) {
      rr = pd->set[s].ranks[i];  /* the \e real rank */
      if (rr == rank) continue; /* skip me */
      pd->sr[pd->is[s]+r] = pd->irank[rr];
      pd->rs[pd->irank[rr]] = s;
    }
    arr = &pd->ents[pd->estart[s]];
    arralloc = pd->sne[s];
    iMesh_getEntities(mi,pd->sh[s],iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,&arr,&arralloc,&arrsize,&err);dICHK(mi,err);
    if (arrsize != pd->sne[s]) dERROR(1,"incorrect sizes");
  } 
  dFunctionReturn(0);
}

static dErr MeshPackDataPrepare(dMesh mesh,PackData *pd,const dMeshTag tag)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  iMesh_Instance mi = mesh->mi;
  dIInt tsize,tlen,ttype;
  PetscMPIInt inc,sheadsize,ssetheadsize;
  dInt pos;
  dErr err;

  dFunctionBegin;
  /* If the tag has changed, recompute all offsets (this is really only necessary if \a tsize has changed) */
  if (pd->tag == tag) dFunctionReturn(0);
  iMesh_getTagSizeValues(mi,tag,&tlen,&err);dICHK(mi,err);
  iMesh_getTagType(mi,tag,&ttype,&err);dICHK(mi,err);
  tsize = tlen * iBase_SizeFromType[ttype];
  if (pd->maxtsize < tsize) {  /* Reallocate the buffers */
    err = PetscFree2(pd->buffer,pd->data);dCHK(err);
    /* Arbitrary packing estimate */
    pd->maxtsize   = tsize;
    pd->buffersize = pd->nr*(dInt)sizeof(SharedHeader) + pd->nsr*(dInt)sizeof(SharedSetHeader) + pd->nsre*tsize + pd->nsr*64;
    pd->datasize   = pd->nsre * tsize;
    pd->minesize   = pd->ne * tsize;
    err = PetscMalloc3(pd->buffersize,char,&pd->buffer,pd->datasize,char,&pd->data,pd->minesize,char,&pd->mine);dCHK(err);
  }
  pd->tag = tag;
  pd->tsize = tsize;
  pd->tlen = tlen;
  pd->ttype = ttype;

  /* Initialize the \a jd array, it is in set-rank-order (the start of each rank within a set) */
  pos = 0;
  for (dInt s=0; s<pd->ns; s++) {
    for (dInt r=0; r<pd->snr[s]; r++) {
      pd->jd[pd->is[s]+r] = pos;
      pos += pd->sne[s] * tsize;
    }
    err = MPI_Pack_size(pd->sne[s]*tsize,MPI_BYTE,comm,&inc);dCHK(err);
    pd->spacksize[s] = inc;
  }

  err = MPI_Pack_size(1,dMPI_SHARED_HEADER,comm,&sheadsize);dCHK(err);
  err = MPI_Pack_size(1,dMPI_SHARED_SET_HEADER,comm,&ssetheadsize);dCHK(err);

  /* Initialize the \a ib array, this is in rank-order (start of the rank) */
  pos = 0;
  for (dInt r=0; r<pd->nr; r++) {
    pd->ib[r] = pos;
    pos += sheadsize + pd->rns[r]*ssetheadsize;
    for (dInt s=0; s<pd->rns[r]; s++) {
      pos += pd->spacksize[pd->rs[pd->ir[r]+s]];
    }
    pd->ibsize[r] = pos - pd->ib[r];
  }
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
      err = dMemcpy((char*)pd->data+pd->jd[pd->is[s]+r],(char*)pd->mine+pd->estart[s],pd->sne[s]*pd->tsize);dCHK(err);
    }
  }
  dFunctionReturn(0);
}

static dErr MeshPackDataPack(dMesh mesh,PackData *pd)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  SharedHeader head;
  SharedSetHeader shead;
  void *outbuf;
  PetscMPIInt position,outsize;
  dErr err;

  dFunctionBegin;
  iMesh_getTagName(mesh->mi,pd->tag,head.tagname,&err,sizeof(head.tagname));dICHK(mesh->mi,err);
  for (dInt r=0; r<pd->nr; r++) {
    outbuf = (char*)pd->buffer + pd->ib[r];
    position = 0;
    outsize = pd->ibsize[r];

    head.numsets = pd->rns[r];    /* The message header */
    err = MPI_Pack(&head,1,dMPI_SHARED_HEADER,outbuf,outsize,&position,comm);dCHK(err);
    for (dInt s=0; s<pd->rns[r]; s++) {
      dInt ss,rr;
      ss = pd->rs[pd->ir[r]+s]; /* The \e real set index */
      shead.numranks = pd->snr[ss];
      err = dMemcpy(shead.ranks,pd->set[ss].ranks,MAX_SHARING_PROCS);dCHK(err);
      err = MPI_Pack(&shead,1,dMPI_SHARED_SET_HEADER,outbuf,outsize,&position,comm);dCHK(err);

      err = dFindInt(shead.numranks,shead.ranks,r,&rr);dCHK(err); /* Find my set-rank index */
      err = MPI_Pack((char*)pd->data+pd->jd[pd->is[ss]+rr],pd->sne[ss]*pd->tsize,MPI_BYTE,outbuf,outsize,&position,comm);dCHK(err);
    }
    pd->ibused[r] = position;
  }
  dFunctionReturn(0);
}

static dErr MeshPackDataUnpack(dMesh mesh,PackData *pd)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  SharedHeader head;
  SharedSetHeader shead;
  char tagname[MAX_NAME_LEN];
  void *inbuf;
  PetscMPIInt position,insize;
  PetscTruth match;
  dErr err;

  dFunctionBegin;
  iMesh_getTagName(mesh->mi,pd->tag,tagname,&err,sizeof(tagname));dICHK(mesh->mi,err);
  for (dInt r=0; r<pd->nr; r++) {
    inbuf = (char*)pd->buffer+pd->ib[r];
    insize = pd->ibsize[r];
    position = 0;

    err = MPI_Unpack(inbuf,insize,&position,&head,1,dMPI_SHARED_HEADER,comm);dCHK(err);
    if (head.numsets != pd->rns[r]) dERROR(1,"Incorrect number of sets in header.");
    err = PetscMemcmp(head.tagname,tagname,sizeof(tagname),&match);dCHK(err);
    if (!match) dERROR(1,"Tag names do not match.");
    for (dInt s=0; s<pd->rns[r]; s++) {
      dInt ss,rr;
      ss = pd->rs[pd->ir[r]+s]; /* The \e real set index */

      err = MPI_Unpack(inbuf,insize,&position,&shead,1,dMPI_SHARED_SET_HEADER,comm);dCHK(err);
      if (shead.numranks != pd->snr[ss]) dERROR(1,"Number of ranks in set does not agree.");
      err = PetscMemcmp(shead.ranks,pd->set[ss].ranks,shead.numranks*sizeof(shead.ranks[0]),&match);dCHK(err);
      if (!match) dERROR(1,"Ranks in set do not agree.");

      err = dFindInt(shead.numranks,shead.ranks,r,&rr);dCHK(err); /* Find my set-rank index */
      err = MPI_Unpack(inbuf,insize,&position,(char*)pd->data+pd->jd[pd->is[ss]+rr],pd->sne[ss]*pd->tsize,MPI_BYTE,comm);dCHK(err);
    }
  }
  dFunctionReturn(0);
}
  
static dErr MeshPackDataUnload(dMesh mesh,PackData *pd)
{
  iMesh_Instance mi = mesh->mi;
  dErr err;

  dFunctionBegin;
  /* Put \a data into \a mine, copy over existing data if there are multiple sets */
  for (dInt s=0; s<pd->ns; s++) {
    for (dInt r=0; r<pd->snr[s]; r++) {
      err = dMemcpy((char*)pd->mine+pd->estart[s],(char*)pd->data+pd->jd[pd->is[s]+r],pd->sne[s]*pd->tsize);dCHK(err);
    }
  }
  iMesh_setArrData(mi,pd->ents,pd->ne,pd->tag,pd->mine,pd->ne*pd->tsize,&err);dICHK(mi,err);
  dFunctionReturn(0);
}

static dErr MeshPackDataPostSends(dMesh mesh,PackData *pd)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  PetscMPIInt mpitag = ((dObject)mesh)->tag;
  dErr err;

  dFunctionBegin;
  for (dInt r=0; r<pd->nr; r++) {
    err = MPI_Isend((char*)pd->buffer+pd->ib[r],pd->ibused[r],MPI_PACKED,pd->rank[r],mpitag,comm,pd->mpireq+r);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr MeshPackDataPostRecvs(dMesh mesh,PackData *pd)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  PetscMPIInt mpitag = ((dObject)mesh)->tag;
  dErr err;

  dFunctionBegin;
  for (dInt r=0; r<pd->nr; r++) {
    err = MPI_Irecv((char*)pd->buffer+pd->ib[r],pd->ibsize[r],MPI_PACKED,pd->rank[r],mpitag,comm,pd->mpireq+r);dCHK(err);
  }
  dFunctionReturn(0);
}

static dErr PackDataDestroy(PackData *pd)
{
  dErr err;

  dFunctionBegin;
  err = PetscFree6(pd->set,pd->is,pd->estart,pd->sne,pd->snr,pd->spacksize);dCHK(err);
  err = PetscFree4(pd->ir,pd->rns,pd->rank,pd->irank);dCHK(err);
  err = PetscFree3(pd->sr,pd->rs,pd->jd);dCHK(err);
  err = PetscFree5(pd->ib,pd->ibsize,pd->ibused,pd->mpireq,pd->ents);dCHK(err);
  err = PetscFree3(pd->buffer,pd->data,pd->mine);dCHK(err);
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
  dErr err;

  dFunctionBegin;
  err = MeshPackDataPrepare(mesh,pack->owned,tag);dCHK(err);
  err = MeshPackDataPrepare(mesh,pack->unowned,tag);dCHK(err);
  err = MeshPackDataLoad(mesh,pack->owned);dCHK(err);
  err = MeshPackDataPack(mesh,pack->owned);dCHK(err);
  err = MeshPackDataPostSends(mesh,pack->owned);dCHK(err);
  err = MeshPackDataPostRecvs(mesh,pack->unowned);dCHK(err);
  err = MPI_Waitall(pack->unowned->nr,pack->unowned->mpireq,MPI_STATUSES_IGNORE);dCHK(err); /* Wait on the receives */
  err = MeshPackDataUnpack(mesh,pack->unowned);dCHK(err);
  err = MeshPackDataUnload(mesh,pack->unowned);dCHK(err);
  err = MPI_Waitall(pack->owned->nr,pack->owned->mpireq,MPI_STATUSES_IGNORE);dCHK(err); /* Make sure the sends all completed */
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

  /* Partition based on ownership */
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

    {
      PetscSynchronizedPrintf(MPI_COMM_WORLD,"owned sets = %d   unowned sets = %d\n",osi,usi);
      PetscSynchronizedFlush(MPI_COMM_WORLD);
    }

    pack->owned->ns = osi;
    pack->unowned->ns = usi;
    err = dMalloc(osi*sizeof(dMeshESH),&pack->owned->sh);dCHK(err);
    err = dMalloc(usi*sizeof(dMeshESH),&pack->unowned->sh);dCHK(err);
    err = dMemcpy(pack->owned->sh,oset,osi*sizeof(oset[0]));dCHK(err);
    err = dMemcpy(pack->unowned->sh,uset,usi*sizeof(uset[0]));dCHK(err);
    err = PetscFree2(oset,uset);dCHK(err);
  }
  err = dFree(allset);dCHK(err);

  err = MeshPackDataSetUp(mesh,pack->owned);dCHK(err);
  err = MeshPackDataSetUp(mesh,pack->unowned);dCHK(err);
  dFunctionReturn(0);
}

static dErr dMeshDestroy_Pack(dMesh mesh)
{
  Mesh_Pack *pack = mesh->data;
  dErr err;

  dFunctionBegin;
  err = PackDataDestroy(pack->owned);dCHK(err);
  err = PackDataDestroy(pack->unowned);dCHK(err);
  err = PetscFree2(pack->owned,pack->unowned);dCHK(err);
  err = dFree(mesh->data);dCHK(err);
  dFunctionReturn(0);
}

static dErr dMeshSetFromOptions_Pack(dMesh mesh)
{
  Mesh_Pack *pack = mesh->data;
  char str[dSTR_LEN];
  dBool flg;
  dErr err;

  dFunctionBegin;
  err = PetscOptionsString("-dmesh_partition","Name of partition tag","dMeshSetInFile",pack->partition,str,sizeof(str),&flg);dCHK(err);
  if (flg) {
    err = PetscStrfree(pack->partition);dCHK(err);
    err = PetscStrallocpy(str,&pack->partition);dCHK(err);
  }
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
  err = PetscMalloc2(1,PackData,&pack->owned,1,PackData,&pack->unowned);dCHK(err);
  err = PetscStrallocpy("PARALLEL_PARTITION",&pack->partition);dCHK(err);

  iMesh_newMesh("PARALLEL",&mesh->mi,&err,(int)strlen("PARALLEL"));dICHK(mesh->mi,err);

  mesh->ops->view           = 0;
  mesh->ops->destroy        = dMeshDestroy_Pack;
  mesh->ops->tagbcast       = dMeshTagBcast_Pack;
  mesh->ops->setfromoptions = dMeshSetFromOptions_Pack;
  mesh->ops->load           = dMeshLoad_Pack;
  dFunctionReturn(0);
}
