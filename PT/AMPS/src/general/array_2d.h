//===================================================
//$Id$
//===================================================

#ifndef ARRAY_2D
#define ARRAY_2D

#include <math.h>
#include <stdlib.h>
#include <iostream.h>

using namespace std;

template <class T>
class array_2d {
private:
  T* data;
  long int size_dim0,size_dim1;

public:

  array_2d() {
    data=NULL;
    size_dim0=0;
    size_dim1=0;
  };

//===================================================
 ~array_2d() { 
   if (data!=NULL) delete [] data;
   data=NULL;
  };

//===================================================
  array_2d(long int n0,long int n1) {
    if ((n0<=0)||(n1<=0)) {
      cerr << "Error: allocation of array_2d object " << endl;
      cerr << "with negative number of elemens" << endl;
      exit(0);
    } 

   data=new T[n0*n1];
   size_dim0=n0;
   size_dim1=n1;
  };

//===================================================
  void init(long int n0,long int n1) {
    if ((n0<=0)||(n1<=0)) {
      cerr << "Error: allocation of array_2d object " << endl;
      cerr << "with negative number of elemens" << endl;
      exit(0);
    }
   
    if (size_dim0!=0) {
      cerr << "Error: initialization of allocated of array_2d object " << endl;
      exit(0);
    }

   data=new T[n0*n1];
   size_dim0=n0;
   size_dim1=n1;
  };
  
//===================================================
  T   operator () (long int i0,long int i1) const { 
    return data[i0+size_dim0*i1]; 
  };

//===================================================
  T & operator () (long int i0,long int i1) { 
    return data[i0+size_dim0*i1]; 
  };

//===================================================
  array_2d<T>& operator = (const array_2d<T>& v) {
    long int i,imax;
  
    imax=size_dim0*size_dim1;
    for (i=0;i<imax;i++) data[i]=v.data[i];
 
    return *this;
  };

//===================================================
  array_2d<T>& operator = (T f) {
    long int i,imax;

    imax=size_dim0*size_dim1;
    for (i=0;i<imax;i++) data[i]=f;
 
    return *this;
  };

//===================================================
  friend array_2d<T> operator + (const array_2d<T> &v1,const array_2d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)) {
      cerr << "Error: add two array_2d<T> of different length.\n";
      exit(0);
    }

    long int i,imax;
    array_2d<T> v3(v1.size_dim0,v1.size_dim1);

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]+v2.data[i]; 
   
    return v3;
  };

//===================================================
  friend array_2d<T> operator - (const array_2d<T> &v1,const array_2d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)) {
      cerr << "Error: add two array_2d<T> of different length.\n";
      exit(0);
    }

    long int i,imax; 
    array_2d<T> v3(v1.size_dim0,v1.size_dim1);

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]-v2.data[i];

    return v3;
  };

//===================================================
  friend array_2d<T> operator * (const array_2d<T> &v1, const T t) {
    long int i,imax;
    array_2d<T> v3(v1.size_dim0,v1.size_dim1);

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v3.data[i]=t*v1.data[i];

    return v3;
  };

//===================================================
  friend array_2d<T> operator / (const array_2d<T> &v1, const T t) {
    if (t == 0) {
      cerr << "Error: divide vector by 0.\n";
      exit(0);
    }
    long int i,imax;
    array_2d<T> v3(v1.size_dim0,v1.size_dim1);

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]/t;

    return v3;
  };

//===================================================
  friend array_2d<T>& operator += (array_2d<T> &v1,const array_2d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)) { 
      cerr << "Error: add two vectors of different length.\n";
      exit(0);
    }

    long int i,imax;

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v1.data[i]+=v2.data[i];

    return v1;
  };

//===================================================
  friend array_2d<T>& operator -= (array_2d<T> &v1,const array_2d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)) {
      cerr << "Error: add two vectors of different length.\n";
      exit(0);
    }

    long int i,imax;

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v1.data[i]-=v2.data[i];

    return v1;
  };

//===================================================
  friend array_2d<T>& operator *= (array_2d<T> &v1,const T t) {
    long int i,imax;

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v1.data[i]*=t;

    return v1;
  };

//===================================================
  friend array_2d<T>& operator /= (array_2d<T> &v1,const T t) {
    if (t == 0) {
      cerr << "Error: divide array_2d<T> by 0.\n";
      exit(0);
    }

    long int i,imax;

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v1.data[i]/=t;

    return v1;
  };

//===================================================
  friend ostream &operator << (ostream &o, const array_2d<T> &v) {
    long int offset;
    for(long int i0=0;i0<v.size_dim0;i0++)
    for(long int i1=0;i1<v.size_dim1;i1++) {
      offset=i0+v.size_dim0*i1;
      o <<v.data[offset]<<" ("<<i0<<","<<i1<<")"<<endl;
    }

    return(o);
  };

};

#endif
