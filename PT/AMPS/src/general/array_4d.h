//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===================================================
//$Id$
//===================================================

#ifndef ARRAY_4D
#define ARRAY_4D

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "specfunc.h"

template <class T>
class array_4d {
protected:
  T* data;
  long int size_dim0,size_dim1,size_dim2,size_dim3;
  long int ndim0_ndim1,ndim0_ndim1_ndim2;

public:

  array_4d() {
    data=NULL;
    size_dim0=0;
    size_dim1=0;
    size_dim2=0;
    size_dim3=0;
  };

//===================================================
 ~array_4d() { 
    if (data!=NULL) delete [] data;
    data=NULL; 
  };

//===================================================
  array_4d(long int n0,long int n1,long int n2,long int n3) {
    if ((n0<=0)||(n1<=0)||(n2<=0)||(n3<=0)) {
      printf("Error: allocation of array_4d object\n");
      printf("with negative number of elemens\n");
      exit(__LINE__,__FILE__);
    } 

   data=new T[n0*n1*n2*n3];
   size_dim0=n0;
   size_dim1=n1;
   size_dim2=n2;
   size_dim3=n3;

   ndim0_ndim1=n0*n1; 
   ndim0_ndim1_ndim2=n0*n1*n2;
  };

//===================================================
  void init(long int n0,long int n1,long int n2,long int n3) {
    if ((n0<=0)||(n1<=0)||(n2<=0)||(n3<=0)) {
      printf("Error: allocation of array_4d object\n");
      printf("with negative number of elemens\n");
      exit(__LINE__,__FILE__);
    }
   
    if (size_dim0!=0) {
      printf("Error: initialization of allocated of array_4d object\n");
      exit(__LINE__,__FILE__);
    }

   data=new T[n0*n1*n2*n3];
   size_dim0=n0;
   size_dim1=n1;
   size_dim2=n2;
   size_dim3=n3;

   ndim0_ndim1=n0*n1;
   ndim0_ndim1_ndim2=n0*n1*n2;
  };
  
//===================================================
  T   operator () (long int i0,long int i1,long int i2,long int i3) const { 
    return data[i0+size_dim0*i1+ndim0_ndim1*i2+ndim0_ndim1_ndim2*i3]; 
  };

//===================================================
  T & operator () (long int i0,long int i1,long int i2,long int i3) { 
    return data[i0+size_dim0*i1+ndim0_ndim1*i2+ndim0_ndim1_ndim2*i3]; 
  };

//===================================================
  array_4d<T>& operator = (const array_4d<T>& v) {
    long int i,imax;
  
    imax=size_dim0*size_dim1*size_dim2*size_dim3;
    for (i=0;i<imax;i++) data[i]=v.data[i];
 
    return *this;
  };

//===================================================
  array_4d<T>& operator = (T f) {
    long int i,imax;

    imax=size_dim0*size_dim1*size_dim2*size_dim3;
    for (i=0;i<imax;i++) data[i]=f;
 
    return *this;
  };

//===================================================
  friend array_4d<T> operator + (const array_4d<T> &v1,const array_4d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)
    ||(v1.size_dim3!=v2.size_dim3)) {
      printf("Error: add two array_4d<T> of different length.\n");
      exit(__LINE__,__FILE__);
    }

    long int i,imax;
    array_4d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2,v1.size_dim3);

    imax=v1.size_dim0*v1.size_dim1*v1.size_dim2*v1.size_dim3;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]+v2.data[i]; 
   
    return v3;
  };

//===================================================
  friend array_4d<T> operator - (const array_4d<T> &v1,const array_4d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)
    ||(v1.size_dim3!=v2.size_dim3)) {
      printf("Error: add two array_4d<T> of different length.\n");
      exit(__LINE__,__FILE__);
    }

    long int i,imax; 
    array_4d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2,v1.size_dim1);

    imax=v1.size_dim0*v1.size_dim1*v1.size_dim2*v1.size_dim3;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]-v2.data[i];

    return v3;
  };

//===================================================
  friend array_4d<T> operator * (const array_4d<T> &v1, const T t) {
    long int i,imax;
    array_4d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2,v1.size_dim3);

    imax=v1.size_dim0*v1.size_dim1*v1.size_dim2*v1.size_dim3;
    for (i=0;i<imax;i++) v3.data[i]=t*v1.data[i];

    return v3;
  };

//===================================================
  friend array_4d<T> operator / (const array_4d<T> &v1, const T t) {
    if (t == 0) {
      printf("Error: divide vector by 0.\n");
      exit(__LINE__,__FILE__);
    }
    long int i,imax;
    array_4d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2,v1.size_dim3);

    imax=v1.size_dim0*v1.size_dim1*v1.size_dim2*v1.size_dim3;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]/t;

    return v3;
  };

//===================================================
  friend array_4d<T>& operator += (array_4d<T> &v1,const array_4d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2) 
    ||(v1.size_dim3!=v2.size_dim3)) {
      printf("Error: add two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    long int i,imax;

    imax=v1.size_dim0*v1.size_dim1*v1.size_dim2*v1.size_dim3;
    for (i=0;i<imax;i++) v1.data[i]+=v2.data[i];

    return v1;
  };

//===================================================
  friend array_4d<T>& operator -= (array_4d<T> &v1,const array_4d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)
    ||(v1.size_dim3!=v2.size_dim3))  {
      printf("Error: add two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    long int i,imax;

    imax=v1.size_dim0*v1.size_dim1*v1.size_dim2*v1.size_dim3;
    for (i=0;i<imax;i++) v1.data[i]-=v2.data[i];

    return v1;
  };

//===================================================
  friend array_4d<T>& operator *= (array_4d<T> &v1,const T t) {
    long int i,imax;

    imax=v1.size_dim0*v1.size_dim1*v1.size_dim2*v1.size_dim3;
    for (i=0;i<imax;i++) v1.data[i]*=t;

    return v1;
  };

//===================================================
  friend array_4d<T>& operator /= (array_4d<T> &v1,const T t) {
    if (t == 0) {
      printf("Error: divide array_4d<T> by 0.\n");
      exit(__LINE__,__FILE__);
    }

    long int i,imax;

    imax=v1.size_dim0*v1.size_dim1*v1.size_dim2*v1.size_dim3;
    for (i=0;i<imax;i++) v1.data[i]/=t;

    return v1;
  };

//===================================================
};

#endif
