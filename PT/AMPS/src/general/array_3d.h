//===================================================
//$Id$
//===================================================

#ifndef ARRAY_3D
#define ARRAY_3D

#include <math.h>
#include <stdlib.h>
#include <iostream.h>

using namespace std;

template <class T>
class array_3d {
private:
  T* data;
  long int size_dim0,size_dim1,size_dim2;
  long int ndim0_ndim1;

public:

  array_3d() {
    data=NULL;
    size_dim0=0;
    size_dim1=0;
    size_dim2=0;
  };

//===================================================
 ~array_3d() {
   if (data!=NULL) delete [] data;
   data=NULL;
 };

//===================================================
  array_3d(long int n0,long int n1,long int n2) {
    if ((n0<=0)||(n1<=0)||(n2<=0)) {
      cerr << "Error: allocation of array_3d object " << endl;
      cerr << "with negative number of elemens" << endl;
      exit(0);
    } 

   data=new T[n0*n1*n2];
   size_dim0=n0;
   size_dim1=n1;
   size_dim2=n2;

   ndim0_ndim1=n0*n1; 
  };

//===================================================
  long int size() const {
    return size_dim0*size_dim1*size_dim2;
  }

//===================================================
  void init(long int n0,long int n1,long int n2) {
    if ((n0<=0)||(n1<=0)||(n2<=0)) {
      cerr << "Error: allocation of array_3d object " << endl;
      cerr << "with negative number of elemens" << endl;
      exit(0);
    }
   
    if (size_dim0!=0) {
      cerr << "Error: initialization of allocated of array_3d object " << endl;
      exit(0);
    }

   data=new T[n0*n1*n2];
   size_dim0=n0;
   size_dim1=n1;
   size_dim2=n2;

   ndim0_ndim1=n0*n1;
  };
  
//===================================================
  T   operator () (long int i0,long int i1,long int i2) const { 
    return data[i0+size_dim0*i1+ndim0_ndim1*i2]; 
  };

//===================================================
  T & operator () (long int i0,long int i1,long int i2) { 
    return data[i0+size_dim0*i1+ndim0_ndim1*i2]; 
  };

//===================================================
  array_3d<T>& operator = (const array_3d<T>& v) {
    long int offset;

    for(long int i0=0;i0<v.size_dim0;i0++) 
    for(long int i1=0;i1<v.size_dim1;i1++)
    for(long int i2=0;i2<v.size_dim2;i2++) {
      offset=i0+v.size_dim0*i1+v.ndim0_ndim1*i2;
      data[offset]=v.data[offset];
    }

    return *this;
  };

//===================================================
  array_3d<T>& operator = (T f) {
    long int offset;
 
    for(long int i0=0;i0<size_dim0;i0++)
    for(long int i1=0;i1<size_dim1;i1++)
    for(long int i2=0;i2<size_dim2;i2++) {
      offset=i0+size_dim0*i1+ndim0_ndim1*i2;
      data[offset]=f;
    }

    return *this;
  };

//===================================================
  friend array_3d<T> operator + (const array_3d<T> &v1,const array_3d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)) {
      cerr << "Error: add two array_3d<T> of different length.\n";
      exit(0);
    }

    long int offset;
    array_3d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2);
    for(long int i0=0;i0<v1.size_dim0;i0++)
    for(long int i1=0;i1<v1.size_dim1;i1++)
    for(long int i2=0;i2<v1.size_dim2;i2++) {
      offset=i0+v1.size_dim0*i1+v1.ndim0_ndim1*i2;
      v3.data[offset]=v1.data[offset]+v2.data[offset];
    }

    return v3;
  };

//===================================================
  friend array_3d<T> operator - (const array_3d<T> &v1,const array_3d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)) {
      cerr << "Error: add two array_3d<T> of different length.\n";
      exit(0);
    }

    long int offset;
    array_3d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2);
    for(long int i0=0;i0<v1.size_dim0;i0++)
    for(long int i1=0;i1<v1.size_dim1;i1++)
    for(long int i2=0;i2<v1.size_dim2;i2++) {
      offset=i0+v1.size_dim0*i1+v1.ndim0_ndim1*i2;
      v3.data[offset]=v1.data[offset]-v2.data[offset];
    }

    return v3;
  };

//===================================================
  friend array_3d<T> operator * (const array_3d<T> &v1, const T t) {
    long int offset;
    array_3d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2);
    for(long int i0=0;i0<v1.size_dim0;i0++)
    for(long int i1=0;i1<v1.size_dim1;i1++)
    for(long int i2=0;i2<v1.size_dim2;i2++) {
      offset=i0+v1.size_dim0*i1+v1.ndim0_ndim1*i2;
      v3.data[offset]=t*v1.data[offset];
    }
    
    return v3;
  };

//===================================================
  friend array_3d<T> operator / (const array_3d<T> &v1, const T t) {
    if (t == 0) {
      cerr << "Error: divide vector by 0.\n";
      exit(0);
    }
    long int offset;
    array_3d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2);
    for(long int i0=0;i0<v1.size_dim0;i0++)
    for(long int i1=0;i1<v1.size_dim1;i1++)
    for(long int i2=0;i2<v1.size_dim2;i2++) {
      offset=i0+v1.size_dim0*i1+v1.ndim0_ndim1*i2;
      v3.data[offset]=v1.data[offset]/t;
    }

    return v3;
  };

//===================================================
  friend array_3d<T>& operator += (array_3d<T> &v1,const array_3d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)) {
      cerr << "Error: add two vectors of different length.\n";
      exit(0);
    }
    long int offset;
    for(long int i0=0;i0<v1.size_dim0;i0++)
    for(long int i1=0;i1<v1.size_dim1;i1++)
    for(long int i2=0;i2<v1.size_dim2;i2++) {
      offset=i0+v1.size_dim0*i1+v1.ndim0_ndim1*i2;
      v1.data[offset]+=v2.data[offset];
    }

    return v1;
  };

//===================================================
  friend array_3d<T>& operator -= (array_3d<T> &v1,const array_3d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)) {
      cerr << "Error: add two vectors of different length.\n";
      exit(0);
    }
    long int offset;
    for(long int i0=0;i0<v1.size_dim0;i0++)
    for(long int i1=0;i1<v1.size_dim1;i1++)
    for(long int i2=0;i2<v1.size_dim2;i2++) {
      offset=i0+v1.size_dim0*i1+v1.ndim0_ndim1*i2;
      v1.data[offset]-=v2.data[offset];
    }

    return v1;
  };

//===================================================
  friend array_3d<T>& operator *= (array_3d<T> &v1,const T t) {
    long int offset;
    for(long int i0=0;i0<v1.size_dim0;i0++)
    for(long int i1=0;i1<v1.size_dim1;i1++)
    for(long int i2=0;i2<v1.size_dim2;i2++) {
      offset=i0+v1.size_dim0*i1+v1.ndim0_ndim1*i2;
      v1.data[offset]*=t;
    }
    return v1;
  };

//===================================================
  friend array_3d<T>& operator /= (array_3d<T> &v1,const T t) {
    if (t == 0) {
      cerr << "Error: divide array_3d<T> by 0.\n";
      exit(0);
    }

    long int offset;
    for(long int i0=0;i0<v1.size_dim0;i0++)
    for(long int i1=0;i1<v1.size_dim1;i1++)
    for(long int i2=0;i2<v1.size_dim2;i2++) {
      offset=i0+v1.size_dim0*i1+v1.ndim0_ndim1*i2;
      v1.data[offset]/=t;
    }

    return v1;
  };

//===================================================
  friend ostream &operator << (ostream &o, const array_3d<T> &v) {
    long int offset;
    for(long int i0=0;i0<v.size_dim0;i0++)
    for(long int i1=0;i1<v.size_dim1;i1++)
    for(long int i2=0;i2<v.size_dim2;i2++) {
      offset=i0+v.size_dim0*i1+v.ndim0_ndim1*i2;
      o <<v.data[offset]<<" ("<<i0<<","<<i1<<","<<i2<<")"<<endl;
    }

    return(o);
  };

};

#endif
