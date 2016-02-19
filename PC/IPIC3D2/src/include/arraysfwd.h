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

/* forward declaration for array classes */
#ifndef arraysfwd_h
#define arraysfwd_h
#include "ipicdefs.h" // for pfloat

// Unfortunately C++ does not allow typedefs for templates:
//
// template <class type> typedef type* array_fetch1;
// template <class type> typedef type** array_fetch2;
// template <class type> typedef type*** array_fetch2;
//
// So to make operator[] return an animal that I can assign to
// in a way that will work independent of whether the underlying
// array access is by means of a chained pointer or a calculated
// 1-dimensional index, I resort to declaration macros.
// 
#if defined(FLAT_ARRAYS) || defined(CHECK_BOUNDS)

  //#define const_arr_ref3(type) iPic3D::const_array_ref3<type>
  //#define const_arr_ref2(type) iPic3D::const_array_ref2<type>
  //#define const_arr_ref1(type) iPic3D::const_array_ref1<type>
  #define const_arr_get3(type) iPic3D::const_array_get3<type>
  #define const_arr_get2(type) iPic3D::const_array_get2<type>
  #define const_arr_get1(type) iPic3D::const_array_get1<type>
  #define arr_fetch3(type) iPic3D::array_fetch3<type>
  #define arr_fetch2(type) iPic3D::array_fetch2<type>
  #define arr_fetch1(type) iPic3D::array_fetch1<type>
  
#else // FLAT_ARRAYS

  #define const_arr_get3(type) type const*const*const*
  #define const_arr_get2(type) type const*const*
  #define const_arr_get1(type) type const*
  #define arr_fetch3(type) type***
  #define arr_fetch2(type) type**
  #define arr_fetch1(type) type*

#endif

#if defined(FLAT_ARRAYS) || defined(CHECK_BOUNDS)
  #define convert_to_arr3(arg) (arg.fetch_arr3())
#else
  #define convert_to_arr3(arg) (arg)
#endif
//template <class type>
//inline type*** fetch_arr3(array_fetch3(type)&in)
//{
//  #if defined(FLAT_ARRAYS) || defined(CHECK_BOUNDS)
//  return in.fetch_arr3();
//  #else
//  return in; // the argument is what we want, so return it
//  #endif
//}

namespace iPic3D
{
  template <class T>
  class const_array_ref3;
  template <class T>
  class const_array_ref4;
  template <class T>
  class array_ref1;
  template <class T>
  class array_ref2;
  template <class T>
  class array_ref3;
  template <class T>
  class array_ref4;
  template <class T>
  class array1;
  template <class T>
  class array2;
  template <class T>
  class array3;
  template <class T>
  class array4;
  template <class T>
  class array_fetch1;
  template <class T>
  class const_array_get1;
  template <class T>
  class const_array_get2;
  template <class T>
  class const_array_get3;
  template <class T>
  class array_fetch2;
  template <class T>
  class array_fetch3;
}

// These aliases are defined for the following flexibilization purposes:
// - to avoid filling the code with template brackets
//   (i.e., to minimize explicitly template-dependent code).
// - so that they can be redefined according to the user's
//   preferred array implementation.
//
typedef iPic3D::array_ref1<int> arr1_int;
typedef iPic3D::array_ref2<int> arr2_int;
typedef iPic3D::array_ref3<int> arr3_int;
typedef iPic3D::array_ref4<int> arr4_int;
//
typedef iPic3D::const_array_ref3<void*> const_arr3_ptr;
typedef iPic3D::array_ref3<void*> arr3_ptr;
//
typedef iPic3D::const_array_ref3<double> const_arr3_double;
typedef iPic3D::const_array_ref4<double> const_arr4_double;
typedef iPic3D::const_array_ref4<pfloat> const_arr4_pfloat;
typedef iPic3D::array_ref1<double> arr1_double;
typedef iPic3D::array_ref2<double> arr2_double;
typedef iPic3D::array_ref3<double> arr3_double;
typedef iPic3D::array_ref4<double> arr4_double;
typedef iPic3D::array1<int> array1_int;
typedef iPic3D::array2<int> array2_int;
typedef iPic3D::array3<int> array3_int;
typedef iPic3D::array4<int> array4_int;
typedef iPic3D::array1<double> array1_double;
typedef iPic3D::array2<double> array2_double;
typedef iPic3D::array3<double> array3_double;
typedef iPic3D::array4<double> array4_double;
typedef iPic3D::array4<pfloat> array4_pfloat;
// This directive should be consistent with the directives in Alloc.h
//#if defined(FLAT_ARRAYS) || defined(CHECK_BOUNDS)
//typedef iPic3D::array_fetch1<double> arr1_double_fetch;
//typedef iPic3D::array_fetch2<double> arr2_double_fetch;
//typedef iPic3D::array_fetch3<double> arr3_double_fetch;
//#else
//typedef double*    arr1_double_fetch;
//typedef double**   arr2_double_fetch;
//typedef double***  arr3_double_fetch;
//#endif
typedef arr_fetch1(double) arr1_double_fetch;
typedef arr_fetch2(double) arr2_double_fetch;
typedef arr_fetch3(double) arr3_double_fetch;

typedef const_arr_get1(double) arr1_double_get;
typedef const_arr_get2(double) arr2_double_get;
typedef const_arr_get3(double) arr3_double_get;
typedef const_arr_get1(pfloat) arr1_pfloat_get;

#endif
