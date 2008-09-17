#ifndef SPECFUNC
#define SPECFUNC

#include <math.h>
#include <iostream>

#ifdef MPI_ON
#include "mpi.h"
#endif

extern int DIM;
extern int ThisThread;
extern int TotalThreadsNumber;

long int nint(double);

void rnd_seed();
double rnd();

/* gamma and error functions used from g++ math library
double erf(double);
double gam(double);
*/

void PrintErrorLog(const char*);
void PrintErrorLog(long int,const char*,const char*);

void StampSignature(char*);
void exit(long int,const char*,const char* =NULL);
void PrintLineMark(long int,char*,char* =NULL);

template<class T>
void PrintLineMark(long int nline ,char* fname ,T code) {
  long thread;
  char *buffer=new char[sizeof(T)*TotalThreadsNumber];

#ifdef MPI_ON
  MPI_Gather((char*)&code,sizeof(T),MPI_CHAR,buffer,sizeof(T),MPI_CHAR,0,MPI_COMM_WORLD);
#else
  char *ptr;
  int i;  

  for (ptr=(char*)&code,i=0;i<sizeof(T);i++,ptr++) buffer[i]=*ptr;  
#endif

  if (ThisThread==0) {
    std::cout << "linemark: line=" << nline << ", file=" << fname, " << code="; 
    for (thread=0;thread<TotalThreadsNumber;thread++) std::cout << "  " << ((T*)buffer)[thread]; 

    std::cout << std::endl;
  }

  delete [] buffer;
}



template<class TMesh>
bool GetGradient(double* gradQ,double cellQ,double* Q,long int ncell,TMesh &grid) {
  int counter,idim,pface;
  long int neib;
  double dx2,x[4][3],x0[3],df[4];
  double A[3][3],aa[3][3],af[3],detaa;

  switch (DIM) {
  case 1:
    grid.cell[ncell].GetCellCenter(x0);
    counter=0,gradQ[0]=0.0,dx2=0.0;

    for (pface=0;pface<DIM+1;pface++) if ((neib=grid.cell[ncell].neighbour_cellno[pface])>=0) {
      counter++;

      grid.cell[neib].GetCellCenter(x[pface]);
      gradQ[0]+=(Q[pface]-cellQ)*(x[pface][0]-x0[0]);
      dx2=pow(x[pface][0]-x0[0],2);
    }

    gradQ[0]/=(counter!=0) ? dx2 : 1.0; 

    if (counter!=DIM+1) return false;
    break;
  case 2:
    grid.cell[ncell].GetCellCenter(x0);
    counter=0,gradQ[0]=0.0,gradQ[1]=0.0;

    for (pface=0;pface<DIM+1;pface++) if ((neib=grid.cell[ncell].neighbour_cellno[pface])>=0) {
      counter++;

      grid.cell[neib].GetCellCenter(x[pface]);
      A[pface][0]=x[pface][0]-x0[0],A[pface][1]=x[pface][1]-x0[1];
      df[pface]=Q[pface]-cellQ;
    }

    if (counter!=DIM+1) return false;

    aa[0][0]=pow(A[0][0],2)+pow(A[1][0],2)+pow(A[2][0],2);
    aa[0][1]=A[0][0]*A[0][1]+A[1][0]*A[1][1]+A[2][0]*A[2][1];
    aa[1][0]=aa[0][1];
    aa[1][1]=pow(A[0][1],2)+pow(A[1][1],2)+pow(A[2][1],2);

    af[0]=A[0][0]*df[0]+A[1][0]*df[1]+A[2][0]*df[2];
    af[1]=A[0][1]*df[0]+A[1][1]*df[1]+A[2][1]*df[2];

    detaa=aa[0][0]*aa[1][1]-aa[1][0]*aa[0][1];

    gradQ[0]=(af[0]*aa[1][1]-af[1]*aa[0][1])/detaa;
    gradQ[1]=(aa[0][0]*af[1]-aa[1][0]*af[0])/detaa;

    break;
  default:
    printf("Error: GetGradient. DIM=%i is not implemented\n",DIM);
    exit(__LINE__,__FILE__);
  }

  return true;
}

//=========================================================
//calculation of CRC-32
class CRC32 {
private:
  unsigned long crc_accum,crc_table[256];

  //generate the table of CRC remainders for all possible bytes 
  void generare_crc_table() { 
    register int i, j;  
    register unsigned long crc_accum;

    for (i=0;i<256;i++) { 
      crc_accum=((unsigned long)i<<24);
      for (j=0;j<8;j++) crc_accum=(crc_accum&0x80000000L) ? (crc_accum<<1)^0x04c11db7L : crc_accum<<1;
      crc_table[i]=crc_accum; 
    }
  };

public: 

  CRC32 () {
    crc_accum=0;
    generare_crc_table();
  };

  //update the CRC on the data block one byte at a time
  template <class T> 
  void add(T* buffer, long int size) {
    char *data_blk_ptr=(char*)buffer;
    long int data_blk_size=size*sizeof(T); 
    register long int i,j;

    for (j=0;j<data_blk_size;j++) { 
      i=((int)(crc_accum>>24)^ *data_blk_ptr++)&0xff;
      crc_accum=(crc_accum<<8)^crc_table[i]; 
    }
  } 

  template <class T>
  void add(T t) {
    add(&t,1);
  } 

  void clear() {
    crc_accum=0;
  };

  unsigned long checksum() { 
    return crc_accum;
  }; 

  void PrintChecksum(long int nline,char* fname) {
    char message[1000];
    
    sprintf(message," line=%ld, file=%s",nline,fname);
    PrintChecksum(message);
  };

  void PrintChecksum(char* message=NULL) {
    unsigned long int *buffer=new unsigned long int[TotalThreadsNumber];
    long int thread;

    buffer[0]=checksum();

#ifdef MPI_ON
    MPI_Gather(buffer,1,MPI_UNSIGNED_LONG,buffer,1,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);
#endif

    if (ThisThread==0) {
      CRC32 cumulativeSignature;
      cumulativeSignature.add(buffer,TotalThreadsNumber);

      if (message!=NULL) printf("CRC32 checksum, cumulativeSignature=0x%lx, message=%s:\n",cumulativeSignature.checksum(),message);
      else printf("CRC32 checksum, cumulativeSignature=0x%lx:\n",cumulativeSignature.checksum());

      for (thread=0;thread<TotalThreadsNumber;thread++) printf("thread=%ld, sum=0x%lx\n",thread,buffer[thread]);
    }

    delete [] buffer;
  };

};

 
  

#endif
   
