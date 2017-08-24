//$Id$
//array of bitwise flags

#include <stdio.h>
#include <math.h>
#include <cmath>

#include "specfunc.h"

#ifndef _DO_NOT_LOAD_GLOBAL_H_
#include "global.h"
#endif

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#include <omp.h>
#endif //_PIC_COMPILATION_MODE_ == _PIC_COMPILATION_MODE__HYBRID_

class cBitArray {
public:
  int size; //size of the bit array
  
private:
  int *data; //the ptr to the integeter array

public:
  cBitArray() {
    size = 0;
    data =NULL;
  }

  cBitArray(int num) {
    if (num>0) {
      size = num;
      int nint = (int) ceil((double)num/(sizeof(int)*8));
      data = new int [nint];
      for (int i=0; i<nint; i++ ) data[i]=0;
    }
    else exit(__LINE__,__FILE__,"Error: parameter 'num' must be more that zero");
  }

  ~cBitArray() {  // destructor
    if (data!=NULL) delete [] data;
  }

  bool SetTrue(int k) {  // set the k-th element in the bit array to true
    if (k<0 || k>=size) return false;
    
    int i = k/(8*sizeof(int));
    int j = k%(8*sizeof(int));
    
    unsigned int flag = 1;
    flag = flag << j;
    int & tmp=data[i];

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    #pragma omp atomic
    #endif
    tmp |=  flag;

    return true;
  }

  bool SetFalse(int k) { //set the k-th element in the bit array to false
    if (k<0 || k>=size) return false;
    
    int i = k/(8*sizeof(int));
    int j = k%(8*sizeof(int));
    
    unsigned int flag = 1;
    flag = flag << j;
    flag = ~flag;
    int & tmp = data[i];

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    #pragma omp atomic
    #endif
    tmp &=  flag;

    return true;
  }
  
  // get the value (true or false) of the k-th element of the bit array
  bool Get(int k) {
    if ((k<0) || (k>=size)) {
      char msg[200];

      sprintf(msg,"Error: %d is out of the range of the bit array", k);
      exit(__LINE__,__FILE__,msg);
    }

    int i = k/(8*sizeof(int));
    int j = k%(8*sizeof(int));
    
    unsigned int flag = 1;
    flag = flag << j;
    
    return ((data[i] & flag)==0) ? false : true;
  }
};
