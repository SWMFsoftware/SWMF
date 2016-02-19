/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite  
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// input_array.h
// A sample user-defined data type for illustrating ConfigFile
// Operators << and >> are defined to allow writing to and reading from files
// Richard J. Wagner 24 May 2004
// Modified P. Henri 8 June 2011
// corrected by Markidis

#ifndef INPUT_ARRAY_H
#define INPUT_ARRAY_H

#include <iostream>
#include <iosfwd>

struct array_int {
  int a, b, c, d, e, f,g,h,i,j,k,l;

  array_int() {}
  array_int(int u1, int u2, int u3, int u4, int u5, int u6, int u7, int u8, int u9, int u10, int u11, int u12):
	  	      a(u1), b(u2), c(u3), d(u4), e(u5), f(u6), g(u7), h(u8), i(u9), j(u10), k(u11), l(u12){}
  array_int(const array_int & orig):a(orig.a), b(orig.b), c(orig.c), d(orig.d), e(orig.e), f(orig.f), g(orig.g), h(orig.h), i(orig.i), j(orig.j), k(orig.k), l(orig.l){}

  array_int & operator=(const array_int & orig) {
    a = orig.a;
    b = orig.b;
    c = orig.c;
    d = orig.d;
    e = orig.e;
    f = orig.f;
    g = orig.g;
    h = orig.h;
    i = orig.i;
    j = orig.j;
    k = orig.k;
    l = orig.l;
    return *this;
  }
};


inline std::ostream & operator<<(std::ostream & os, const array_int & t) {

  os << t.a << " " << t.b << " " << t.c << " " << t.d << " " << t.e << " " << t.f<< " " << t.g<< " " << t.h<< " " << t.i<< " " << t.j<< " " << t.k<< " " << t.l;
  return os;
}


inline std::istream & operator>>(std::istream & is, array_int & t) {

  is >> t.a >> t.b >> t.c >> t.d >> t.e >> t.f >> t.g >> t.h >> t.i >> t.j >> t.k >> t.l;
  return is;
}


struct array_double {
  double a, b, c, d, e, f,g,h,i,j,k,l;

  array_double() {}

  array_double(int u1, int u2, int u3, int u4, int u5, int u6, int u7, int u8, int u9, int u10, int u11, int u12):
	               a(u1), b(u2), c(u3), d(u4), e(u5), f(u6), g(u7), h(u8), i(u9), j(u10), k(u11), l(u12){}

  array_double(const array_double & orig):a(orig.a), b(orig.b), c(orig.c), d(orig.d), e(orig.e), f(orig.f), g(orig.g), h(orig.h), i(orig.i), j(orig.j), k(orig.k), l(orig.l){}

  array_double & operator=(const array_double & orig) {
    a = orig.a;
    b = orig.b;
    c = orig.c;
    d = orig.d;
    e = orig.e;
    f = orig.f;
    g = orig.g;
    h = orig.h;
    i = orig.i;
    j = orig.j;
    k = orig.k;
    l = orig.l;
    return *this;
  }
};


inline std::ostream & operator<<(std::ostream & os, const array_double & t) {

  os << t.a << " " << t.b << " " << t.c << " " << t.d << " " << t.e << " " << t.f << " " << t.g<< " " << t.h<< " " << t.i<< " " << t.j<< " " << t.k<< " " << t.l;
  return os;
}


inline std::istream & operator>>(std::istream & is, array_double & t) {
  is >> t.a >> t.b >> t.c >> t.d >> t.e >> t.f >> t.g >> t.h >> t.i >> t.j >> t.k >> t.l;
  return is;
}


struct array_bool {
  bool a, b, c, d, e, f,g,h,i,j,k,l;

  array_bool() {}
  array_bool(int u1, int u2, int u3, int u4, int u5, int u6, int u7, int u8, int u9, int u10, int u11, int u12):
	  	  	   a(u1), b(u2), c(u3), d(u4), e(u5), f(u6), g(u7), h(u8), i(u9), j(u10), k(u11), l(u12){}
  array_bool(const array_bool & orig):a(orig.a), b(orig.b), c(orig.c), d(orig.d), e(orig.e), f(orig.f),g(orig.g), h(orig.h), i(orig.i), j(orig.j), k(orig.k), l(orig.l){}

  array_bool & operator=(const array_bool & orig) {
    a = orig.a;
    b = orig.b;
    c = orig.c;
    d = orig.d;
    e = orig.e;
    f = orig.f;
    g = orig.g;
	h = orig.h;
	i = orig.i;
	j = orig.j;
	k = orig.k;
	l = orig.l;
    return *this;
  }
};


inline std::ostream & operator<<(std::ostream & os, const array_bool & t) {
  os << t.a << " " << t.b << " " << t.c << " " << t.d << " " << t.e << " " << t.f << " " << t.g<< " " << t.h<< " " << t.i<< " " << t.j<< " " << t.k<< " " << t.l;
  return os;
}


inline std::istream & operator>>(std::istream & is, array_bool & t) {
  is >> t.a >> t.b >> t.c >> t.d >> t.e >> t.f >> t.g >> t.h >> t.i >> t.j >> t.k >> t.l;
  return is;
}

#endif
