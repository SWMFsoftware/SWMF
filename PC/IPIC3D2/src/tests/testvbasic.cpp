
// I started with Reza Rahman's
// "Intel Xeon Phi Coprocessor Architecture and Tools:
//   The Guide for Application Developers".
// I referred to:
// * shuffle intrinsics documentation:
//     http://software.intel.com/en-us/node/460872
// * vector mask intrinsics:
//     http://goo.gl/SLFmT7
// * intel xeon phi instruction set architecture reference manual:
//    http://software.intel.com/sites/default/files/forum/278102/327364001en.pdf
// * for double precision the mask bits are the bottom 8 bits:
//    http://stackoverflow.com/questions/18742863/manipulating-masks-for-doubles-on-xeon-phi
// * intrinsics for type casting:
//    http://software.intel.com/en-us/node/485381
// * masked swizzle parameters:
//    http://goo.gl/IXNDdP
// * SIMD MIC classes, e.g. F64vec8, M512 (__m512):
//    http://goo.gl/vMKWo5
// * SIMD Xeon classes (seem to be better documented):
//    http://vhdplus.sourceforge.net/doc/class_f64vec4.html

#define MICVEC_DEFINE_OUTPUT_OPERATORS
#include <iostream>
#include <micvec.h>
//#include <immintrin.h> // for __mm512_int2mask(0xf);
#include <stdio.h>
#include <stdint.h> // for int64_t
#define NO_MPI
#include "../utility/debug.cpp"
#include "Alloc.h"
#include "mic_basics.h"

#define printexpr(var) std::cout << "line " << __LINE__ << ": " \
  << #var << " = " << var << std::endl;

#define print_F64vec8(in) \
  { \
  double* var = (double*) &in; \
  printf("line %d: %s = (", __LINE__, #in); \
  for(int i=0; i<7; i++) printf("%g, ", var[i]); \
  printf("%g)\n", var[7]); \
  }

#define print_int_arr(var) \
  { \
  printf("line %d: %s = (", __LINE__, #var); \
  for(int i=0; i<15; i++) printf("%d, ", var[i]); \
  printf("%d)\n", var[15]); \
  }

using namespace std;

void testI8()
{
  _MM_PERM_ENUM p32;
  __declspec(align(64)) I64vec8 inputData(7,6,5,4,3,2,1,0);
  __declspec(align(64)) I64vec8 outputData;

  int64_t* idata = (int64_t*) &inputData;
  for(int i=0;i<8;i++) printf("idata[%d]=%d\n",i,idata[i]);

  std::cout << "input = " << inputData << endl;
  printexpr(inputData);

  // swizzle input data and print
  //std::cout << "swizzle data for pattern 'cdab' \n" << inputData.cdab() << endl;
  printexpr(inputData.cdab()); // use for 2x2 matrix transpose
  printexpr(inputData.badc()); // use for 4x4 matrix transpose of 2x2 blocks
  printexpr(inputData.dacb()); // use for cross product
  printexpr(inputData.aaaa());
  printexpr(inputData.bbbb());
  printexpr(inputData.cccc());
  printexpr(inputData.dddd());
}

void testF8()
{
  _MM_PERM_ENUM p32;
  __declspec(align(64)) F64vec8 data1(1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.8);
  __declspec(align(64)) F64vec8 data2(2.7,2.6,2.5,2.4,2.3,2.2,2.1,2.8);
  __declspec(align(64)) F64vec8 outp1;
  __declspec(align(64)) F64vec8 outp2;
  __declspec(align(64)) F64vec8 outputData;

  double* idata = (double*) &data1;
  for(int i=0;i<8;i++) printf("idata[%d]=%g\n",i,idata[i]);

  std::cout << "input = " << data1 << endl;

  // swizzle input data and print
  //std::cout << "swizzle data for pattern 'cdab' \n" << inputData.cdab() << endl;
  printexpr(data1.cdab());
  printexpr(data1.badc()); // use for 2x2 transpose
  printexpr(data1.dacb());
  printexpr(data1.aaaa());
  printexpr(data1.bbbb());
  printexpr(data1.cccc());
  printexpr(data1.dddd());
}

void print_data(double data[8][8], int line)
{
  printf("=== at line %d: data is: ===\n", line);
  for(int i=0; i<8; i++)
  {
    printf("  ");
    for(int j=0; j<7; j++)
    {
      printf("%g, ", data[i][j]);
    }
    printf("%g)\n", data[i][7]);
  }
  printf("=== end of data ===\n");
}

void test_transpose_8x8_double()
{
  __declspec(align(64)) double data[8][8];
  // initialize with indexing data
  for(int i=0;i<8;i++)
  for(int j=0;j<8;j++)
    data[i][j]=(1+i)+.1*(j+1);

  print_data(data,__LINE__);
  transpose_8x8_double(data);
  print_data(data,__LINE__);
}

// copy applying mask
inline void masked_copy(__m512d& dst, __m512d src, __mmask8 mask)
{
    dst = _mm512_mask_mov_pd(dst,mask,src);
}

int test_copy_methods()
{
  // const causes compiler error
  /*const*/ F64vec8 u(8,7,6,5,4,3,2,1);
  F32vec16 v = _mm512_cvtpd_pslo(u); 
  F32vec16 vinv = F32vec16(1.)/v;
  printexpr(vinv);
  F64vec8 vinvd = _mm512_cvtpslo_pd(vinv);
  printexpr(vinvd);
  printexpr(F64vec8(1.)/u);
  // store the part of the stream before the first 64-bit
  // byte-aligned address preceding 
  F64vec8 arr = F64vec8(0.);
  //double arr[8]ALLOC_ALIGNED;
  //for(int i=0;i<8;i++) arr[i] = 0.;
  print_F64vec8(u);
  print_F64vec8(arr);
  copy012to012(arr, u);
  print_F64vec8(arr);
  copy456to456(arr, u);
  print_F64vec8(arr);
  copy456to012(arr, u);
  print_F64vec8(arr);
  copy012to456(arr, u);
  print_F64vec8(arr);
  //_mm512_packstorehi_epi32(&arr[4],*(__m512i*)&u);
  //print_F64vec8(arr);
  // _mm512_packstorelo_epi32(&arr[5],*(__m512i*)&u);

  const F64vec8 u1(8.1,7.1,6.1,5.1,4.1,3.1,2.1,1.1);
  const F64vec8 u2(8.2,7.2,6.2,5.2,4.2,3.2,2.2,1.2);
  // see mask_values in mic_basics.h
  const __mmask8 rmask_11101110 = _mm512_int2mask(0x77); //mask=01110111
  const F64vec8 out = _mm512_mask_mov_pd(u1,rmask_11101110,u2);
  print_F64vec8(out);
}

// calculate the cross product of the corresponding
// physical vectors in the two halves of the arguments:
//
// for(int i=0;i<8;i+=4)
// {
//   w[i+0] = u[i+1]*v[i+2] - u[i+2]*v[i+1];
//   w[i+1] = u[i+2]*v[i+0] - u[i+0]*v[i+2];
//   w[i+2] = u[i+0]*v[i+1] - u[i+1]*v[i+0];
// }
void test_cross_product()
{
  // calculate u cross v
  F64vec8 u(1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1);
  F64vec8 v(2.8,2.7,2.6,2.5,2.4,2.3,2.2,2.1);
  F64vec8 w = cross_product(u,v);
  printexpr(w);
  for(int i=0;i<8;i+=4)
  {
    w[i+0] = u[i+1]*v[i+2] - u[i+2]*v[i+1];
    w[i+1] = u[i+2]*v[i+0] - u[i+0]*v[i+2];
    w[i+2] = u[i+0]*v[i+1] - u[i+1]*v[i+0];
  }
  printexpr(w);
}

void test_integer_double_conversions()
{
  const F64vec8 u(8.8,7.7,-6.6,-5.5,-4.4,3.3,2.2,1.1);
  printexpr(u);
  // Rounding control values; these can be one of the following:
  // * _MM_ROUND_MODE_NEAREST - round to nearest (even)
  // * _MM_ROUND_MODE_DOWN - round toward negative infinity
  // * _MM_ROUND_MODE_UP - round toward positive infinity
  // * _MM_ROUND_MODE_TOWARD_ZERO - round toward zero
  
  // theoretically this returns an unsigned int, but it seems
  // to work for a signed int as well.  The official function
  // to call is _mm512_cvt_roundpd_epi32lo, as documented at
  // http://software.intel.com/en-us/node/461034, but it seems
  // that this function does not actually exist.
  const MyI32vec16 v = _mm512_cvtfxpnt_roundpd_epi32lo(u, _MM_ROUND_MODE_DOWN);
  const MyI32vec16 u_nearest = _mm512_cvtfxpnt_roundpd_epi32lo(u, _MM_ROUND_MODE_NEAREST);
  // why is this guy messed up?
  F64vec8 uround = _mm512_cvtepu32lo_pd(v);
  printexpr(uround);
  uround = floor(u);
  printexpr(uround);
  uround = max(F64vec8(0.),uround);
  printexpr(uround);
  int* varr = (int*)&v;
  print_int_arr(varr);
  int* u_nearest_arr = (int*)&u_nearest;
  print_int_arr(u_nearest_arr);
}

inline MyI32vec16 init_I32vec16(int x, int y, int z, int t)
{ return MyI32vec16(0,0,0,0,0,0,0,0,t,z,y,x,t,z,y,x); }
inline MyI32vec16 init_I32vec16(int x, int y, int z)
{ return MyI32vec16(0,0,0,0,0,0,0,0,0,z,y,x,0,z,y,x); }

void test_init()
{
  I32vec16 a4 = init_I32vec16(1,2,3,4);
  I32vec16 a3 = init_I32vec16(1,2,3);
  a4 = a4 + I32vec16(1);
  int* a4_arr = (int*)&a4;
  int* a3_arr = (int*)&a3;
  print_int_arr(a4_arr);
  print_int_arr(a3_arr);
  //I32vec16 a(1,2);
  //int* a_arr = (int*)&a;
  //print_int_arr(a_arr);
}

void test_maxmin()
{
  F64vec8 u(-8,7,-6,-5,4,3,2,1);
  F64vec8 v(9,-8,-6,3,4,3,2,1);
  F64vec8 M = _mm512_max_pd(u,v);
  F64vec8 gM = _mm512_gmax_pd(u,v);
  print_F64vec8(M);
  print_F64vec8(gM);
}

void test_compare()
{
  I32vec16 u(0,0,0,0,0,0,0,0,8,7,-6,5,4,3,2,1);
  I32vec16 v(0,0,0,0,0,0,0,0,8,7,-6,5,-4,3,2,1);
  bool eq = !test_ne(u,v);
  dprint(eq);
  //I32vec16 M = _mm512_max_pd(u,v);
  //I32vec16 m = _mm512_min_pd(u,v);
  // returns a mask
  __mmask16 mask_ne = _mm512_cmp_epi32_mask(u,v,_MM_CMPINT_NE);
  __mmask16 mask_lt = _mm512_cmp_epi32_mask(u,v,_MM_CMPINT_LT);
  __mmask16 mask_ge = _mm512_cmp_epi32_mask(u,v,_MM_CMPINT_GE);
  __mmask16 mask_eq = _mm512_cmp_epi32_mask(u,v,_MM_CMPINT_EQ);
  int int_ne = _mm512_mask2int(mask_ne);
  int int_lt = _mm512_mask2int(mask_lt);
  int int_ge = _mm512_mask2int(mask_ge);
  int int_eq = _mm512_mask2int(mask_eq);
  dprint(int_ne);
  dprint(int_lt);
  dprint(int_ge);
  dprint(int_eq);
  // print_F64vec8(_mm512_max_pd(u,v));
  // print_F64vec8(_mm512_gmax_pd(u,v));
}

int main()
{
  //testI8();
  //testF8();
  test_transpose_8x8_double();
  //test_copy_methods();
  //test_cross_product();
  //test_integer_double_conversions();
  //test_init();
  test_compare();
}
