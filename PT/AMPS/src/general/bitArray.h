//$Id$
//array of bitwise flags

#include <stdio.h>
#include <math.h>
#include <cmath>

class cBitArray {

public:
  int size;
  
private:
  int * data;

public:
  cBitArray(){
    size = 0;
    data =NULL;
  }

  cBitArray(int num){
    if (num>0){
      size = num;
      int nint = (int) ceil((double)num/(sizeof(int)*8));
      data = new int [nint];
      for (int i=0; i<nint; i++ ) data[i]=0;
    }
  }

  ~cBitArray(){
    if (data!=NULL) delete [] data;
  }

  bool SetTrue(int k){
    if (k<0 || k>=size) return false;
    
    int i = k/(8*sizeof(int));
    int j = k%(8*sizeof(int));
    
    unsigned int flag = 1;
    flag = flag << j;
    int & tmp=data[i];
#pragma omp atomic
    tmp |=  flag;
    return true;
  }

  bool SetFalse(int k){
    if (k<0 || k>=size) return false;
    
    int i = k/(8*sizeof(int));
    int j = k%(8*sizeof(int));
    
    unsigned int flag = 1;
    flag = flag << j;
    flag = ~flag;
    int & tmp = data[i];
#pragma omp atomic    
    tmp &=  flag;
    return true;
  }
 
  bool Get(int k){
    if (k<0 || k>=size) return false;
    
    int i = k/(8*sizeof(int));
    int j = k%(8*sizeof(int));
    
    unsigned int flag = 1;
    flag = flag << j;
    
    if (data[i] & flag){
      return true;
    }else{
      return false;
    }
  }
  

};
