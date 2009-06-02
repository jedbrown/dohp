#include <stdio.h>
#include <MBCore.hpp>
#include <MBRange.hpp>

#define CHK(err) do {                                                   \
    if (err != MB_SUCCESS) {                                            \
      fprintf(stderr,"Error %d in %s at %s:%d\n",err,__FUNCTION__,__FILE__,__LINE__); \
      return err;                                                       \
    }                                                                   \
  } while (0)

MBErrorCode infer_entities(MBInterface *iface,int from_dim,int to_dim)
{
  MBRange from,to;
  MBErrorCode err;
  err = iface->get_entities_by_dimension(0,from_dim,from);CHK(err);
  err = iface->get_adjacencies(from,to_dim,true,to);CHK(err);
  return MB_SUCCESS;
}

int main(int argc,char *argv[])
{
  char *infile,*outfile;
  MBInterface* iface = new MBCore;
  MBErrorCode err;

  if (argc != 3) return 1;
  infile = argv[1];
  outfile = argv[2];
  err = iface->load_mesh(infile);CHK(err);
  err = infer_entities(iface,3,2);CHK(err);
  err = infer_entities(iface,3,1);CHK(err);
  err = infer_entities(iface,2,1);CHK(err);

  iface->write_mesh(outfile);
  delete iface;
  return 0;
}
