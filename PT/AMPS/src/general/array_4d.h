//===================================================
//$Id$
//===================================================

#ifndef ARRAY_4D
#define ARRAY_4D

#include <math.h>
#include <stdlib.h>
#include <iostream.h>

using namespace std;

template <class T>
class array_4d {
private:
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
      cerr << "Error: allocation of array_4d object " << endl;
      cerr << "with negative number of elemens" << endl;
      exit(0);
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
      cerr << "Error: allocation of array_4d object " << endl;
      cerr << "with negative number of elemens" << endl;
      exit(0);
    }
   
    if (size_dim0!=0) {
      cerr << "Error: initialization of allocated of array_4d object " << endl;
      exit(0);
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
      cerr << "Error: add two array_4d<T> of different length.\n";
      exit(0);
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
      cerr << "Error: add two array_4d<T> of different length.\n";
      exit(0);
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
      cerr << "Error: divide vector by 0.\n";
      exit(0);
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
      cerr << "Error: add two vectors of different length.\n";
      exit(0);
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
      cerr << "Error: add two vectors of different length.\n";
      exit(0);
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
      cerr << "Error: divide array_4d<T> by 0.\n";
      exit(0);
    }

    long int i,imax;

    imax=v1.size_dim0*v1.size_dim1*v1.size_dim2*v1.size_dim3;
    for (i=0;i<imax;i++) v1.data[i]/=t;

    return v1;
  };

//===================================================
  friend ostream &operator << (ostream &o, const array_4d<T> &v) {
    long int offset;
    for(long int i0=0;i0<v.size_dim0;i0++)
    for(long int i1=0;i1<v.size_dim1;i1++)
    for(long int i2=0;i2<v.size_dim2;i2++) 
    for(long int i3=0;i3<v.size_dim3;i3++) {
      offset=i0+v.size_dim0*i1+v.ndim0_ndim1*i2+v.ndim0_ndim1_ndim2*i3;
      o <<v.data[offset]<<" ("<<i0<<","<<i1<<","<<i2<<","<<i3<<")"<<endl;
    }

    return(o);
  };

};

#endif
