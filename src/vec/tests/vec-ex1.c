static const char help[] = "Test vector closure\n\n";

#include <dohpvec.h>
#include <dohpsys.h>
#include <dohp.h>

int main(int argc,char *argv[])
{
  MPI_Comm     comm;
  dMPIInt      rank,size;
  PetscViewer  viewer;
  Vec          x,xc,xl,y,yc,z;
  IS           isg;
  VecScatter   scatter;
  dScalar     *a;
  dInt         nghost = 2,ghosts[2],n=3,bs=2,xn,xbs;
  dErr         err;

  err = dInitialize(&argc,&argv,0,help);dCHK(err);
  comm = PETSC_COMM_WORLD;
  err = MPI_Comm_size(comm,&size);dCHK(err);
  err = MPI_Comm_rank(comm,&rank);dCHK(err);
  if (size != 2) dERROR(PETSC_COMM_SELF,1,"This example must be run with 2 processes");
  err = PetscViewerASCIIGetStdout(comm,&viewer);dCHK(err);

  ghosts[0] = n*((size+rank-1)%size)+1; /* second block of left neighbor, periodically */
  ghosts[1] = n*((size+rank+1)%size);   /* first block of right neighbor, periodically */
  err = PetscSynchronizedPrintf(comm,"[%d] ghosts %D %D\n",rank,ghosts[0],ghosts[1]);dCHK(err);
  err = PetscSynchronizedFlush(comm);dCHK(err);

  err = VecCreateDohp(comm,bs,n-1,n,nghost,ghosts,&x);dCHK(err);
  err = VecGetLocalSize(x,&xn);dCHK(err);
  if (xn != (n-1)*bs) dERROR(PETSC_COMM_SELF,1,"local size %d, expected %d",xn,(n-1)*bs);
  err = VecGetBlockSize(x,&xbs);dCHK(err);
  if (xbs != bs) dERROR(PETSC_COMM_SELF,1,"block size %d, expected %d",xbs,bs);

  err = PetscPrintf(comm,"Empty vector\n");
  err = VecView(x,viewer);dCHK(err);

  err = VecGetArray(x,&a);dCHK(err);
  for (dInt i=0; i<xn; i++) a[i] = 1000+rank*100+i;
  err = VecRestoreArray(x,&a);dCHK(err);
  err = PetscPrintf(comm,"Global unclosed form\n");dCHK(err);
  err = VecView(x,viewer);dCHK(err);
  err = VecDohpGetClosure(x,&xc);dCHK(err);
  err = PetscPrintf(comm,"Global closed form\n");dCHK(err);
  err = VecView(xc,viewer);dCHK(err);
  err = VecGhostGetLocalForm(xc,&xl);dCHK(err);

  /* Create a global vector that will hold the local values on each process, and a scatter from that to the local form.
  * This makes viewing deterministic, where as using PETSC_VIEWER_STDOUT_SELF requires manual synchronization. */
  err = VecCreateMPI(comm,(n+nghost)*bs,PETSC_DECIDE,&z);dCHK(err);
  err = ISCreateStride(comm,(n+nghost)*bs,(n+nghost)*bs*rank,1,&isg);dCHK(err);
  err = VecScatterCreate(z,isg,xl,NULL,&scatter);dCHK(err);
  err = ISDestroy(&isg);dCHK(err);
  err = VecScatterBegin(scatter,xl,z,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd(scatter,xl,z,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = PetscPrintf(comm,"Before VecGhostUpdateBegin/End\n");dCHK(err);
  err = VecView(z,viewer);dCHK(err);
  err = VecGhostUpdateBegin(xc,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecGhostUpdateEnd(xc,INSERT_VALUES,SCATTER_FORWARD);dCHK(err);
  err = VecScatterBegin(scatter,xl,z,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = VecScatterEnd(scatter,xl,z,INSERT_VALUES,SCATTER_REVERSE);dCHK(err);
  err = PetscPrintf(comm,"After VecGhostUpdateBegin/End\n");dCHK(err);
  err = VecView(z,viewer);dCHK(err);
  err = VecScatterDestroy(&scatter);dCHK(err);
  err = VecDestroy(&z);dCHK(err);
  err = VecGhostRestoreLocalForm(xc,&xl);dCHK(err);

  err = VecDuplicate(x,&y);dCHK(err);
  err = VecCopy(x,y);dCHK(err);
  err = PetscPrintf(comm,"Global unclosed form of y\n");dCHK(err);
  err = VecView(y,viewer);dCHK(err);
  err = VecDohpGetClosure(y,&yc);dCHK(err);
  err = PetscPrintf(comm,"Global closed form of y\n");dCHK(err);
  err = VecView(yc,viewer);dCHK(err);
  err = VecZeroEntries(yc);dCHK(err);
  err = VecShift(xc,1000);dCHK(err);
  err = VecCopy(xc,yc);dCHK(err);
  err = PetscPrintf(comm,"Global closed form of y, after copy of shifted closure\n");dCHK(err);
  err = VecView(yc,viewer);dCHK(err);

  err = VecDohpRestoreClosure(x,&xc);dCHK(err);
  err = VecDohpRestoreClosure(y,&yc);dCHK(err);
  err = VecDestroy(&x);dCHK(err);
  err = VecDestroy(&y);dCHK(err);
  err = dFinalize();dCHK(err);
  return 0;
}
