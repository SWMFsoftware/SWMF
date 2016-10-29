//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
#ifndef SPECFUNC
#define SPECFUNC

#include "mpi.h"

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "rnd.h"
#include "constants.h"


#define _STDOUT_ERRORLOG_MODE__ON_   0
#define _STDOUT_ERRORLOG_MODE__OFF_  1
#define _STDOUT_ERRORLOG_MODE_ _STDOUT_ERRORLOG_MODE__ON_


//extern int DIM;
extern int ThisThread;
extern int TotalThreadsNumber;

long int nint(double);



//void rnd_seed(int seed=-1);
//double rnd();




/////////////////////
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
void PrintLineMark(long int nline ,char* fname ,T code,const char *msg=NULL) {
  long thread;
  char *buffer=new char[sizeof(T)*TotalThreadsNumber];

#ifdef MPI_ON
  MPI_Gather((char*)&code,sizeof(T),MPI_CHAR,buffer,sizeof(T),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
#else
  char *ptr;
  int i;  

  for (ptr=(char*)&code,i=0;i<sizeof(T);i++,ptr++) buffer[i]=*ptr;  
#endif

  if (ThisThread==0) {
    std::cout << "linemark: line=" << nline << ", file=" << fname;
    if (msg!=NULL) std::cout << ", msg=" << msg;
    std::cout << ": code=";

    for (thread=0;thread<TotalThreadsNumber;thread++) std::cout << "  " << ((T*)buffer)[thread]; 

    std::cout << std::endl;
  }

  delete [] buffer;
}


#ifdef DIM
template<class TMesh>
bool GetGradient(double* gradQ,double cellQ,double* Q,long int ncell,TMesh &grid) {
  int counter,pface;
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
    printf("$PREFIX:Error: GetGradient. DIM=%i is not implemented\n",DIM);
    exit(__LINE__,__FILE__);
  }

  return true;
}
#endif

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

  void PrintChecksum(long int nline,const char* fname) {
    char message[1000];
    
    sprintf(message," line=%ld, file=%s",nline,fname);
    PrintChecksum(message);
  };

  void PrintChecksum(char* message=NULL) {
    unsigned long int *buffer=new unsigned long int[TotalThreadsNumber];
    long int thread;

    buffer[0]=checksum();

#ifdef MPI_ON
    unsigned long int bufferRecv[TotalThreadsNumber];

    MPI_Gather(buffer,1,MPI_UNSIGNED_LONG,bufferRecv,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    memcpy(buffer,bufferRecv,TotalThreadsNumber*sizeof(unsigned long int));
#endif

    if (ThisThread==0) {
      CRC32 cumulativeSignature;
      cumulativeSignature.add(buffer,TotalThreadsNumber);

      if (message!=NULL) printf("$PREFIX:CRC32 checksum, cumulativeSignature=0x%lx, message=%s:\n",cumulativeSignature.checksum(),message);
      else printf("$PREFIX:CRC32 checksum, cumulativeSignature=0x%lx:\n",cumulativeSignature.checksum());

      for (thread=0;thread<TotalThreadsNumber;thread++) printf("$PREFIX:thread=%ld, sum=0x%lx\n",thread,buffer[thread]);
    }

    delete [] buffer;
  };

};


//=========================================================
//Vector Rotations
namespace VectorRotation {
  inline void Along_Z_direction(double *l,double angle) {
    double temp[2],cosAngle,sinAngle;

    cosAngle=cos(angle);
    sinAngle=sin(angle);

    temp[0]=cosAngle*l[0]-sinAngle*l[1];
    temp[1]=sinAngle*l[0]+cosAngle*l[1];

    memcpy(l,temp,2*sizeof(double));
  }
}


//=========================================================
//Vector Operations
namespace Vector3D {
  inline void CrossProduct(double *res,double *a,double *b) {
    res[0]=a[1]*b[2]-a[2]*b[1];
    res[1]=a[2]*b[0]-a[0]*b[2];
    res[2]=a[0]*b[1]-a[1]*b[0];
  }

  inline double Length(double *x) {
    return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  }

  inline double DotProduct(double *a,double *b) {
    int i;
    double res=0.0;

    for (i=0;i<3;i++) res+=a[i]*b[i];
    return res;
  }

  inline void Orthogonalize(double *PrimaryVector,double *OrthogonalVector) {
     double l2=0.0,c=0.0;
     int i;

     //get the dot products of the vectors 
     for (i=0;i<3;i++) {
        l2+=PrimaryVector[i]*PrimaryVector[i];
        c+=PrimaryVector[i]*OrthogonalVector[i];
     }

     //get the orthogonal component of 'OrthogonalVector'
     c/=l2;

     for (i=0;i<3;i++) {
       OrthogonalVector[i]-=c*PrimaryVector[i];
      }
   }

  inline double ParallelComponentLength(double *Vector,double *Axis) {
    int i;
    double l=0.0,c=0.0;

    for (i=0;i<3;i++) {
      l+=pow(Axis[i],2);
      c+=Vector[i]*Axis[i];
    }

    return c/sqrt(l);
  }
 

  //determine an orthogonal frame of rederence: z is input; e1 and e2 and orthogonal to z and form a right handed frame of reference
  inline void GetNormFrame(double *e0,double *e1,double *z) {
    double l,e[3];
    int idim;

    l=Length(z);
    for (idim=0;idim<3;idim++) e[idim]=z[idim]/l;

    //get e0
    if (fabs(e[0])>1.0E-5) e0[0]=-e[1],e0[1]=e[0],e0[2]=0.0;
    else e0[0]=0.0,e0[1]=e[2],e0[2]=-e[1];

    l=Length(e0);
    for (idim=0;idim<3;idim++) e0[idim]/=l;

    //get e1: e1=z \times e0
    CrossProduct(e1,z,e0);
  }

  inline void Normalize(double *x) {
    double l=Length(x);

    for (int idim=0;idim<3;idim++) x[idim]/=l;
  }

  //distribute the vector direction
  namespace Distribution {
    //uniform distribution of the
    inline void Uniform(double *a) {
      for (int i=0;i<3;i++) a[i]=sqrt(-log(rnd()))*cos(2.0*Pi*rnd());
      Vector3D::Normalize(a);
    }
  }
}
  

//=========================================================
//Relativistic functions
namespace Relativistic {
  inline double Speed2E(double Speed,double mass) {
    return mass*pow(SpeedOfLight,2)*(1.0/sqrt(1.0-pow(Speed/SpeedOfLight,2))-1.0);
  }

  inline double E2Speed(double E,double mass) {
    double mc2;

    mc2=mass*SpeedOfLight*SpeedOfLight;
    return SpeedOfLight*sqrt(E*(E+2.0*mc2))/(E+mc2);
  }
}

#endif
   
