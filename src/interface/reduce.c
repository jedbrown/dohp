#include "private/dohpimpl.h"
#include <MBParallelConventions.h>

/* This macro is in MBParallelComm.hpp and these should be the same */
#define MAX_SHARING_PROCS 64
#define MAX_NAME_LEN 128

struct dMeshSharedSet {
  dMeshESH handle;
  dInt owner,numcopies;
  dInt copyranks[MAX_SHARING_PROCS];
}

struct dMeshSharedHeader {
  PetscMPIInt numtags,numents,datasize,owner,numcopies;
  PetscMPIInt copyranks[MAX_SHARING_PROCS];
};

struct dMeshSharedTagHeader {
  PetscMPIInt type,numvalues;
  char name[MAX_NAME_LEN];
};

static dBool CreatedMPITypes = false;
static MPI_Datatype MPITypeFromITAPS[4];
static MPI_Datatype dMPI_SHARED_HEADER,dMPI_SHARED_TAG_HEADER,dMPI_ENT_HANDLE;

dErr CreateMPITypes()
{
  const MPI_Datatype typeSH[] = {MPI_INT,MPI_INT};
  const MPI_Aint dispSH[] = {offsetof(struct dMeshSharedHeader,numtags),offsetof(struct dMeshSharedHeader,copyranks)};
  const int sizeSH[] = {5,MAX_SHARING_PROCS};
  const MPI_Datatype typeSTH[] = {MPI_INT,MPI_CHAR};
  const MPI_Aint dispSTH[] = {offsetof(struct dMeshSharedHeader,type),offsetof(struct dMeshSharedHeader,name)};
  const int sizeSTH[] = {2,MAX_NAME_LEN};
  dErr err;

  dFunctionBegin;
  if (CreatedMPITypes) dFunctionReturn(0);
  err = MPI_Type_create_struct(2,sizeSH,dispSH,typeSH,&dMPI_SHARED_HEADER);dCHK(err);
  err = MPI_Type_create_struct(2,sizeSTH,dispSTH,typeSTH,&dMPI_SHARED_TAG_HEADER);dCHK(err);
  err = MPI_Type_dup(MPI_UNSIGNED_LONG,dMPI_ENT_HANDLE);dCHK(err);
  MPITypeFromITAPS[0] = MPI_INT;
  MPITypeFromITAPS[1] = MPI_DOUBLE;
  MPITypeFromITAPS[2] = dMPI_ENT_HANDLE;
  MPITypeFromITAPS[3] = MPI_BYTE;
  CreatedMPITypes = true;
  dFunctionReturn(0);
}

/** 
* Determine an upper bound on the amount of space required to pack the tags
* 
* @param mesh the mesh object
* @param sset the shared set to pack, used to determine the number of entities
* @param ntags number of tags to be packed
* @param tag tag handles to pack
* @param bytes upper bound on the number of bytes required
* 
* @return 
*/
static dErr packTagsSize(dMesh mesh,const struct dMeshSharedSet *sset,dInt ntags,const dMeshTag tag[],int *bytes)
{
  iMesh_Instance mi = mesh->mi;
  int type,size;
  dInt needed,inc;
  dErr err;

  dFunctionBegin;
  err = MPI_Pack_size(1,dMPI_SHARED_HEADER,comm,&inc);dCHK(err);
  iMesh_getNumOfType(mi,sset->handle,iBase_ALL_TYPES,&count,&err);dICHK(mi,err);
  needed = inc;
  for (int i=0; i<ntags; i++) {
    iMesh_getTagType(mi,tag[i],&type,&err);dICHK(mi,err);
    iMesh_getTagSizeValues(mi,tag[i],&size,&err);dICHK(mi,err);
    err = MPI_Pack_size(1,dMPI_SHARED_TAG_HEADER,comm,&inc);dCHK(err);
    needed += inc;
    err = MPI_Pack_size(count*size,MPITypeFromITAPS[type],comm,&inc);dCHK(err);
    needed += inc;
  }
  *bytes = needed;
  dFunctionReturn(0);
}

/** 
* Pack the values associated with \a tag on the entity set in \a sset into \a outbuf.
* 
* @param mesh the mesh object
* @param sset struct holding a shared set (normally interface values)
* @param ntags number of tags to pack
* @param tag tag handles to pack
* @param outbuf buffer to pack into
* @param outsize size of buffer in bytes
* @param position position in the buffer (in bytes)
* 
* @return err
*/
static dErr packHeader(dMesh mesh,const struct dMeshSharedSet *sset,dInt ntags,const dMeshTag tag[],void *outbuf,int outsize,int *position)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  iMesh_Instance mi = mesh->mi;
  struct dMeshSharedTagHeader tagheader;
  struct dMeshSharedHeader header;
  MeshListEH ents=MLZ;
  void *values;
  int maxtagsize,numcopies,datasize,valuesAlloc,valuesSize;
  dErr err;

  dFunctionBegin;
  iMesh_getEntities(mi,sset->handle,iBase_ALL_TYPES,iMesh_ALL_TOPOLOGIES,MLREF(ents),&err);dICHK(mi,err);

  /* find out how much space we need for the largest tag values */
  for (int i=0; i<ntags; i++) {
    int bytes;
    iMesh_getTagSizeBytes(mi,tag[i],&bytes,&err);dICHK(mi,err);
    if (maxtagsize < bytes) maxtagsize = bytes;
    datasize += bytes*ents.s;
  }
  valuesAlloc = maxtagsize*ents.s;
  err = dMalloc(valuesAlloc,&values);dCHK(err);

  /* Determine how many copies we have. */
  numcopies = 0;
  for (int i=0; i<MAX_SHARING_PROCS; i++) {
    if (sharing[i] < 0) break;
    numcopies++;
  }
         
  /* Set up the message header. */
  header.numtags   = (PetscMPIInt)ntags;
  header.numents   = (PetscMPIInt)ents.s;
  header.datasize  = (PetscMPIInt)datasize;
  header.owner     = (PetscMPIInt)sset->owner;
  header.numcopies = (PetscMPIInt)sset->numcopies;
  for (int i=0; i<MAX_SHARING_PROCS; i++) {
    header.copyranks[i] = (PetscMPIInt)sset->copyranks[i];
  }

  /* Pack the header */
  err = MPI_Pack(&header,1,dMPI_SHARED_HEADER,outbuf,outsize,position,comm);dCHK(err);

  /* Pack each tag and values */
  for (int i=0; i<ntags; i++) {
    iMesh_getTagType(mi,tag[i],&type,&err);dICHK(mi,err);
    iMesh_getTagSizeValues(mi,tag[i],&size,&err);dICHK(mi,err);
    iMesh_getArrData(mi,ents.v,ents.s,tag[i],(char**)&values,&valuesAlloc,&valuesSize,&err);dICHK(mi,err);

    /* set up the tag header */
    tagheader.type      = type;
    tagheader.numvalues = size;
    iMesh_getTagName(mi,tag[i],tagheader.name,&err,MAX_NAME_LEN);dICHK(mi,err);

    /* pack buffer */
    err = MPI_Pack(&tagheader,1,dMPI_SHARED_TAG_HEADER,outbuf,outsize,position,comm);dCHK(err);
    err = MPI_Pack(values,valuesSize,MPITypeFromITAPS[type],outbuf,outsize,position,comm);dCHK(err);
  }
  err = dFree(values);dCHK(err);
  MeshListFree(ents);
  dFunctionReturn(0);
}

/** 
* This is to check on the necessary data size if for some reason you don't already know how much room is needed to store
* the tag values.
* 
* @param inbuf buffer with the
* @param insize 
* @param datasize 
* 
* @return 
*/
static dErr peekHeaderSize(dMesh mesh,void *inbuf,int insize,const int *position,dInt *datasize)
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  struct dMeshSharedHeader header;
  int dummypos = *position;
  dErr err;

  dFunctionBegin;
  err = MPI_Unpack(inbuf,insize,&dummypos,&header,1,dMPI_SHARED_HEADER,comm);dCHK(err);
  *datasize = header.datasize;
  dFunctionReturn(0);
}

/** 
* Unpack the tag data.  Does not actually set the tag values, but confirms that they match (at least in debug mode).
* 
* @param mesh mesh object
* @param sset shared set
* @param numtags number of tags
* @param tag array of tags
* @param inbuf buffer to unpack
* @param insize size of buffer in bytes
* @param position current position in buffer
* @param outbuf buffer to unpack into
* @param outsize size in bytes
* @param outposition position in bytes
* @param tagdata array of offsets, length \a numtags
* 
* @return err
*/
static dErr unpackHeader(dMesh mesh,const struct dMeshSharedSet *sset,dInt numtags,const dMeshTag tag[],void *inbuf,int insize,int *position,void *outbuf,int outsize,int *outposition,void *tagdata[])
{
  MPI_Comm comm = ((dObject)mesh)->comm;
  iMesh_Instance mi = mesh->mi;
  struct dMeshSharedTagHeader tagheader;
  struct dMeshSharedHeader header;
  void *values;
  dInt numents;
  PetscMPIInt typesize;
  MPI_Datatype mpitype;
  dErr err;

  dFunctionBegin;
  iMesh_getNumOfType(mi,sset->handle,iBase_ALL_TYPES,&numents,&err);dICHK(mi,err);

  err = MPI_Unpack(inbuf,insize,position,&header,1,dMPI_SHARED_HEADER,comm);dCHK(err);

  /* Sanity checks */
  if (sset->owner != header.owner) dERROR(1,"Owner rank %d does not match unpacked %d",sset->owner,header.owner);
  if (sset->numcopies != header.numcopies) dERROR(1,"Number of copies %d does not match unpacked %d",sset->numcopies,header.numcopies);
  for (int i=0; i<numcopies; i++) {
    if (sset->copyranks[i] != header.copyranks[i]) dERROR(1,"Copy %d rank %d does not match unpacked %d",i,sset->copyranks[i],header.copyranks[i]);
  }
  if (numents != header.numents) dERROR(1,"Number of entities %d does not match unpacked %d",numents,header.numents);
  if (header.datasize > outsize - outposition) dERROR(1,"Not enough space in buffer");
  if (numtags != header.numtags) dERROR(1,"Number of tags %d does not match unpacked %d",numtags,header.numtags);

  /* Unpack each set of tag values into the buffer */
  for (int i=0; i<numtags; i++) {
    err = MPI_Unpack(inbuf,insize,position,&tagheader,1,dMPI_SHARED_TAG_HEADER,comm);dCHK(err);

    /* Sanity checks for this tag header */
    {
      int type,numvalues;
      char name[MAX_NAME_LEN];
      iMesh_getTagType(mi,tag[i],&type,&err);dICHK(mi,err);
      iMesh_getTagSizeValues(mi,tag[i],&numvalues,&err);dICHK(mi,err);
      if (type != tagheader.type) dERROR(1,"Tag type %d does not match unpacked %d",type,tagheader.type);
      if (numvalues != tagheader.numvalues) dERROR(1,"Number of values %d does not match unpacked %d",numvalues,tagheader.numvalues);
      iMesh_getTagName(mi,tag[i],name,&err,MAX_NAME_LEN);dICHK(mi,err);
      if (!strcmp(name,tagheader.name)) dERROR(1,"Tag name %s do not match unpacked %s",name,tagheader.name);
    }

    /* Unpack into the buffer and set the tagdata pointer to point at the start. */
    mpitype = MPITypeFromITAPS[tagheader.type];
    err = MPI_Type_size(mpitype,&typesize);
    tagdata[i] = outbuf + outposition;
    err = MPI_Unpack(inbuf,insize,position,tagdata[i],numvalues*numents,mpitype,comm);dCHK(err);
    outposition += numvalues*numents*typesize;
  }

  dFunctionReturn(0);
}
