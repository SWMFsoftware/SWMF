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

#ifndef mic_basics_h
#define mic_basics_h

// documentation on intrinsics:
//
// * intel preprocessor macros:
//     http://goo.gl/CMS4WY
// * SIMD C++ classes:
//     http://software.intel.com/en-us/node/462716
// * AVX-512 intrinsics (for KNL):
//     http://software.intel.com/en-us/node/485150

#if defined(__MIC__)
  #include <micvec.h>
  #define printexpr(var) std::cout << "line " << __LINE__ << ": " \
    << #var << " = " << var << std::endl;

  // === section: extending intrinsics classes ===

  #define MICVEC_DEFINE_OUTPUT_OPERATORS
  //#include <iosfwd>
  #include <iostream>
  // for some bizarre reason this, along with arithmetic and subscripting
  // operators, is not defined for this class, so I do it here.
  inline std::ostream& operator<<(std::ostream& os, const I32vec16 obj)
  {
    int* ref = (int*)&obj;
    os << "(";
    for(int i=0;i<15;i++)
      os << ref[i] << ",";
    os << ref[15] << ")";
    return os;
  }
  
  // unfortunately the built-in I32vec16 is missing
  // a lot of the functionality in e.g. F64vec8,
  // and unfortunately C++ does not provide a good
  // and simple mechanism to extend the interface of a class,
  // so a lot of boilerplate is necessary.
  class MyI32vec16: public I32vec16
  {
   public:
    //using I32vec16::operator+;
    int& operator[](int i)
    {
      int* memory = (int*)this;
      return memory[i];
    }
    // unfortunately we cannot inherit constructors.
    MyI32vec16(int x15,int x14,int x13,int x12,int x11,int x10,int x9,int x8,
               int x7, int x6, int x5, int x4, int x3, int x2, int x1,int x0):
      I32vec16(x15,x14,x13,x12,x11,x10,x9,x8,x7,x6,x5,x4,x3,x2,x1,x0)
    {}
    MyI32vec16(): I32vec16()
    {}
    MyI32vec16(__m512i in): I32vec16(in)
    {}
    MyI32vec16(int in): I32vec16(in)
    {}
    MyI32vec16 operator*(MyI32vec16 in)
    {
      return _mm512_mulhi_epi32(*this,in);
    }
    MyI32vec16 operator+(MyI32vec16 in)
    {
      return _mm512_add_epi32(*this, in);
    }
    MyI32vec16 operator-(MyI32vec16 in)
    {
      return _mm512_sub_epi32(*this, in);
    }
  };


  // === section: extending intrinsics methods ===
  //
  // In these subsections I follow the same organization as used in
  //   http://software.intel.com/en-us/node/460888

  // === subsection: initialization (load) intrinsics ===
  inline F64vec8 make_F64vec8(double x, double y, double z, double t)
  { return F64vec8(t,z,y,x,t,z,y,x); }
  inline F64vec8 make_F64vec8(double x, double y, double z)
  { return F64vec8(0,z,y,x,0,z,y,x); }
  inline I32vec16 make_I32vec16(int x, int y, int z, int t)
  { return I32vec16(0,0,0,0,0,0,0,0,t,z,y,x,t,z,y,x); }
  inline I32vec16 make_I32vec16(int x, int y, int z)
  { return I32vec16(0,0,0,0,0,0,0,0,0,z,y,x,0,z,y,x); }
  
  // === subsection: arithmetic intrinsics ===

  // add v3 to the product of v1 and v2 and return the result
  inline I32vec16 fused_multiply_add(I32vec16 v1, I32vec16 v2, I32vec16 v3)
  {
    return _mm512_fmadd_epi32(v1,v2,v3);
  }
  inline I32vec16 multiply(I32vec16 v1, I32vec16 v2)
  {
    return _mm512_mullo_epi32(v1,v2);
  }

  inline F64vec8 cross_product(F64vec8 u, F64vec8 v)
  {
    F64vec8 us = u.dacb();
    F64vec8 vs = v.dacb();
    return us*vs.dacb() - us.dacb()*vs;
  }

  // === subsection: conversion intrinsics ===

  inline I32vec16 round_down(F64vec8 vec)
  { return _mm512_cvtfxpnt_roundpd_epi32lo(vec, _MM_ROUND_MODE_DOWN); }
  inline I32vec16 round_to_nearest(F64vec8 vec)
  { return _mm512_cvtfxpnt_roundpd_epi32lo(vec, _MM_ROUND_MODE_NEAREST); }

  // === subsection: compare intrinsics ===

  // set bits of lanes that are not equal
  inline int test_ne(I32vec16 u, I32vec16 v)
  {
    __mmask16 mask_ne = _mm512_cmp_epi32_mask(u,v,_MM_CMPINT_NE);
    return _mm512_mask2int(mask_ne);
  }
  // return true if all bits are eq, else return false
  //inline bool test_eq(I32vec16 u, I32vec16 v)
  //{
  //  return !test_ne(u,v);
  //}

  // set bits of elements where u<v
  inline int test_lt(F64vec8 u, F64vec8 v)
  {
    __mmask8 mask_lt = _mm512_cmp_pd_mask(u,v,_MM_CMPINT_LT);
    return _mm512_mask2int(mask_lt);
  }
  
  // set bits of elements where u>v
  inline int test_gt(F64vec8 u, F64vec8 v)
  {
    __mmask8 mask_gt = _mm512_cmp_pd_mask(u,v,_MM_CMPINT_GT);
    return _mm512_mask2int(mask_gt);
  }

  // === subsection: min/max operations ===

  inline F64vec8 maximum(F64vec8 u, F64vec8 v)
  {
    return _mm512_max_pd(u,v);
    //return _mm512_gmax_pd(u,v);
  }
  inline F64vec8 minimum(F64vec8 u, F64vec8 v)
  {
    return _mm512_min_pd(u,v);
    //return _mm512_gmin_pd(u,v);
  }

  // === subsection: memory and initialization ===

  // mask_values
  //
  // rmask  hex   mask  dec   oct
  // 0000 : 0x0 = 0000 =  0 =   0
  // 1000 : 0x1 = 0001 =  1 =   1
  // 0100 : 0x2 = 0010 =  2 =   2
  // 1100 : 0x3 = 0011 =  3 =   3
  // 0010 : 0x4 = 0100 =  4 =   4
  // 1010 : 0x5 = 0101 =  5 =   5
  // 0110 : 0x6 = 0110 =  6 =   6
  // 1110 : 0x7 = 0111 =  7 =   7
  // 0001 : 0x8 = 1000 =  8 =  01
  // 1001 : 0x9 = 1001 =  9 =  02
  // 0101 : 0xa = 1010 = 10 =  03
  // 1101 : 0xb = 1011 = 11 =  04
  // 0011 : 0xc = 1100 = 12 =  05
  // 1011 : 0xd = 1101 = 13 =  06
  // 0111 : 0xe = 1110 = 14 =  07
  // 1111 : 0xf = 1111 = 15 = 010
  //
  // rmask_11001100 (reversed mask) would equal mask_00110011
  // rmask bits go from lowest to highest, but mask bits
  // go from highest to lowest.

  // return concatenation of low half of lo and low half of hi
  inline __m512 cat_low_halves(__m512 lo, __m512 hi)
  {
    const __mmask16 rmask_00001111 = _mm512_int2mask(0xff00); //mask=11110000
    // low half uses lo, and high half uses low half of hi
    return _mm512_mask_permute4f128_ps(lo, rmask_00001111, hi,_MM_PERM_BADC);
  }
  
  // return concatenation of hgh half of lo and hgh half of hi
  inline __m512 cat_hgh_halves(__m512 lo, __m512 hi)
  {
    const __mmask16 rmask_11110000 = _mm512_int2mask(0x00ff); //mask=00001111
    // high half uses hi, and low half uses high half of lo
    return _mm512_mask_permute4f128_ps(hi, rmask_11110000, lo,_MM_PERM_BADC);
  }
  
  inline F64vec8 cat_low_halves(void* lo, void* hi)
  { return (F64vec8) cat_low_halves(*(__m512*)lo, *(__m512*)hi); }
  inline F64vec8 cat_hgh_halves(void* lo, void* hi)
  { return (F64vec8) cat_hgh_halves(*(__m512*)lo, *(__m512*)hi); }
  //
  inline F64vec8 cat_low_halves(F64vec8 lo, F64vec8 hi)
  { return (F64vec8) cat_low_halves(*(__m512*)(&lo), *(__m512*)(&hi)); }
  inline F64vec8 cat_hgh_halves(F64vec8 lo, F64vec8 hi)
  { return (F64vec8) cat_hgh_halves(*(__m512*)(&lo), *(__m512*)(&hi)); }
  
  // copy low vector of src to low vector of dst
  inline void copy012to012(F64vec8& dst, F64vec8 src)
  {
    const __mmask8 rmask_11100000 = _mm512_int2mask(0x07); // mask=00000111
    // masked copy
    dst = _mm512_mask_mov_pd(dst,rmask_11100000,src);
  }
  // copy hgh vector of src to hgh vector of dst
  inline void copy456to456(F64vec8& dst, F64vec8 src)
  {
    const __mmask8 rmask_00001110 = _mm512_int2mask(0x70); // mask=01110000
    dst = _mm512_mask_mov_pd(dst,rmask_00001110,src);
  }
  // copy low vector of src to hgh vector of dst
  inline void copy012to456(F64vec8& dst, F64vec8 src)
  {
    // mask=0000000000111111 because 32-bit elements
    const __mmask16 rmask_1111110000000000 = _mm512_int2mask(0x03f);
    // write the first three elements of src beginning at dst[4]
    _mm512_mask_packstorelo_epi32(&dst[4],rmask_1111110000000000,*(__m512i*)&src);
  }
  // copy hgh vector of src to low vector of dst
  inline void copy456to012(F64vec8& dst, F64vec8 src)
  {
    // mask=0011111111111111 because 32-bit elements
    const __mmask16 rmask_1111111111111100 = _mm512_int2mask(0x3fff);
    // write all but the highest element to the left of dst[4].
    _mm512_mask_packstorehi_epi32(&dst[4],rmask_1111111111111100,*(__m512i*)&src);
  }
  
  // copy vectors of src to vectors of dst
  inline void copy012and456(F64vec8& dst, F64vec8 src)
  {
    const __mmask8 rmask_11101110 = _mm512_int2mask(0x77); //mask=01110111
    dst = _mm512_mask_mov_pd(dst,rmask_11101110,src);
  }
  
  // broadcast s0 into mask=0 slots, s1 into mask=1 slots
  inline F64vec8 broadcast_mask_blend(double s0, double s1, __mmask8 mask)
  {
    return _mm512_mask_blend_pd(mask, F64vec8(s0), F64vec8(s1));
    //
    // broadcast s0 into vector
    //const F64vec8 t1 = _mm512_extload_pd(&s0,
    //  _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
    //
    // broadcast s1 into same vector with mask
    //return (F64vec8) _mm512_mask_extload_pd(t1, mask, &s1,
    //  _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
  }

  // === section: methods for fast transpose ===

  // transpose each 2x2 block in the two rows of a 2x8 matrix
  // used to transpose the 2x2 blocks of the 8x8 matrix
  //
  // This transpose changes data in the form
  //   [a0 a1 a2 a3 a4 a5 a6 a7]
  //   [b0 b1 b2 b3 b4 b5 b6 b7]
  // into the form
  //   [a0 b0 a2 b2 a4 b4 a6 b6]
  //   [a1 b1 a3 b3 a5 b5 a7 b7]
  //
  inline void trans2x2(F64vec8& out1, F64vec8& out2, F64vec8 in1, F64vec8 in2)
  {
    // out1: copy odds from evens of data2, i.e.:
    //
    //   static const bool odd[] = {0,1,0,1,0,1,0,1};
    //   for(int i=0; i<8; i++) out1[i] = odd[i] ? in2[i-1] : in1[i];
    //
    const __mmask16 rmask_01010101 = _mm512_int2mask(0xaa); // mask=10101010
    // replace unmasked in1 with swizzled in2
    out1 = F64vec8(_mm512_mask_swizzle_pd(__m512d(in1),
      rmask_01010101, __m512d(in2),_MM_SWIZ_REG_CDAB));
    //cout << "swizzle for pattern 'cdab' with rmask=01010101\n"
    //  << out1 << endl;
  
    // out2: copy evens from odds of in1, i.e.:
    //
    //   static const bool even[] = {1,0,1,0,1,0,1,0};
    //   for(int i=0;i<8;i++) out2[i] = even[i]? in1[i+1] : out2[i];
    //
    const __mmask16 rmask_10101010 = _mm512_int2mask(0x55); // mask=01010101
    // replace unmasked in2 with swizzled in1
    out2 = F64vec8(_mm512_mask_swizzle_pd(__m512d(in2),
      rmask_10101010, __m512d(in1),_MM_SWIZ_REG_CDAB));
    //cout << "swizzle for pattern 'cdab' with rmask=10101010\n"
    //  << out2 << endl;
  }
  inline void trans2x2_new(double in1_[8], double in2_[8])
  {
    F64vec8& in1 = (F64vec8&)in1_[0];
    F64vec8& in2 = (F64vec8&)in2_[0];
    // compiler should be smart enough to recognize
    // that only one intermediate buffer is needed:
    //trans2x2(in1, in2, in1, in2);
    // but maybe this is not the case, so copy the data
    // that would be used after being rewritten:
    const F64vec8 buff1 = in1;
    trans2x2(in1, in2, buff1, in2);
  }
  inline void trans2x2(double in1_[8], double in2_[8])
  {
    // copy the data that we will first overwrite
    F64vec8& data1 = (F64vec8&)in1_[0];
    F64vec8& data2 = (F64vec8&)in2_[0];
    const F64vec8 buff1 = data1;
    
    // data1: copy odds from evens of data2,
    //   i.e.,
    // for(int i=1; i<8; i+=2) data1[i] = data2[i-1];
    //
    const __mmask16 rmask_01010101 = _mm512_int2mask(0xaa); // mask=10101010
    // replace unmasked data1 with swizzled data2
    // (hopefully compiler will see that this can be done in
    // a single vector instruction)
    data1 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data1),
      rmask_01010101, __m512d(data2),_MM_SWIZ_REG_CDAB));
    //cout << "swizzle for pattern 'cdab' with rmask=01010101\n"
    //  << data1 << endl;
    
    // data2: copy evens from odds of buff1,
    //   i.e.,
    // for(int i=0; i<7; i+=2) data2[i] = buff1[i+1];
    //
    const __mmask16 rmask_10101010 = _mm512_int2mask(0x55); // mask=01010101
    // replace unmasked data2 with swizzled buff1
    data2 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data2),
      rmask_10101010, __m512d(buff1),_MM_SWIZ_REG_CDAB));
    //cout << "swizzle for pattern 'cdab' with rmask=10101010\n"
    //  << data2 << endl;
  }
  
  // transpose 1x2 elements of each 2x4 block in the two rows of a 2x8 matrix
  // used to transpose the 2x2 elements of the 4x4 blocks of the 8x8 matrix
  //
  // This transpose changes data in the form
  //   [a0 a1 a2 a3 a4 a5 a6 a7]
  //   [b0 b1 b2 b3 b4 b5 b6 b7]
  // into the form
  //   [a0 a1 b0 b1 a4 a5 b4 b5]
  //   [a2 a3 b2 b3 a6 a7 b6 b7]
  //
  inline void trans2x4(F64vec8& out1, F64vec8& out2, F64vec8 in1, F64vec8 in2)
  {
    // data1: copy odd 1x2 elements from even 1x2 elements of in2, i.e.,
    //   static const bool odd[] = {0,0,1,1,0,0,1,1};
    //   for(int i=0; i<8; i++) out1[i] = odd[i] ? in2[i-2] : in1[i];
    //
    // replace unmasked in1 with swizzled in2
    const __mmask16 rmask_00110011 = _mm512_int2mask(0xcc); // mask = 11001100
    out1 = F64vec8(_mm512_mask_swizzle_pd(__m512d(in1),
      rmask_00110011, __m512d(in2),_MM_SWIZ_REG_BADC));
  
    // out2: copy even 1x2 elements from odd 1x2 elements of in1, i.e.,
    //   static const bool even[] = {1,1,0,0,1,1,0,0};
    //   for(int i=0; i<8; i++) out2[i] = even[i] ? in1[i+2] : in2[i];
    //
    const __mmask16 rmask_11001100 = _mm512_int2mask(0x33); // mask=00110011
    out2 = F64vec8(_mm512_mask_swizzle_pd(__m512d(in2),
      rmask_11001100, __m512d(in1),_MM_SWIZ_REG_BADC));
  }
  inline void trans2x4(F64vec8& data1, F64vec8& data2)
  {
    // copy the data that we will first overwrite
    const F64vec8 buff1 = data1;
  
    // data1: copy odd 1x2 elements from even 1x2 elements of data2,
    //   i.e.,
    // for(int i=2; i<8; i++) if((i/2)%2) data1[i] = data2[i-2];
    //
    // replace unmasked data1 with swizzled data2
    const __mmask16 rmask_00110011 = _mm512_int2mask(0xcc); // mask = 11001100
    data1 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data1),
      rmask_00110011, __m512d(data2),_MM_SWIZ_REG_BADC));
    //cout << "swizzle for pattern 'badc' with rmask=00110011\n"
    //  << outp1 << endl;
  
    // data2: copy even 1x2 elements from odd 1x2 from odds of buff1,
    //   i.e.,
    // for(int i=0; i<6; i++) if(!((i/2)%2)) data2[i] = buff1[i+2];
    //
    const __mmask16 rmask_11001100 = _mm512_int2mask(0x33); // mask=00110011
    data2 = F64vec8(_mm512_mask_swizzle_pd(__m512d(data2),
      rmask_11001100, __m512d(buff1),_MM_SWIZ_REG_BADC));
    //cout << "swizzle for pattern 'badc' with rmask=11001100\n"
    //  << outp2 << endl;
  }
  inline void trans2x4_new(F64vec8& data1, F64vec8& data2)
  {
    // copy the data that would be used after being rewritten:
    const F64vec8 buff1 = data1;
    trans2x4(data1, data2, buff1, data2);
  }
  inline void trans2x4(double in1[8], double in2[8])
  {
    F64vec8& data1 = (F64vec8&)in1[0];
    F64vec8& data2 = (F64vec8&)in2[0];
    trans2x4(data1,data2);
  }
  
  // transpose 1x4 elements of each 2x4 block in the two rows of a 2x8 matrix
  // used to transpose the 4x4 elements of the 8x8 matrix
  //
  // This transpose changes data in the form
  //   [a0 a1 a2 a3 a4 a5 a6 a7]
  //   [b0 b1 b2 b3 b4 b5 b6 b7]
  // into the form
  //   [a0 a1 a2 a3 b0 b1 b2 b3]
  //   [a4 a5 a6 a7 b4 b5 b6 b7]
  //
  // for swizzle intel supports 256-bit lanes with 64-bit
  // elements (as well as 128-bit lanes with 32-bit elements),
  // but for shuffle intel only supports 128-bit lanes with
  // 32-bit elements, so we must cast to 32-bit data types, do
  // the shuffle, and cast back.
  //
  inline void trans2x8(F64vec8& out1_, F64vec8& out2_, F64vec8 in1_, F64vec8 in2_)
  {
    __m512& out1 = (__m512&)out1_[0];
    __m512& out2 = (__m512&)out2_[0];
    __m512& in1 = (__m512&)in1_[0];
    __m512& in2 = (__m512&)in2_[0];
    // out1: copy high (odd) 1x4 element from low (even) 1x4 element of in2:
    //   static const bool mask[] = {0,0,0,0,1,1,1,1};
    //   for(int i=0; i<8; i++) out1[i] = mask[i] ? in2[i-4] : in1[i];
    //
    // replace unmasked in1 with shuffled in2
    const __mmask16 rmask_00001111 = _mm512_int2mask(0xff00); //mask=11110000
    out1 = _mm512_mask_permute4f128_ps(in1, rmask_00001111, in2,_MM_PERM_BADC);
  
    // out2: copy low (even) 1x4 element from high (odd) 1x4 element of in1,
    //   static const bool mask[] = {1,1,1,1,0,0,0,0};
    //   for(int i=0; i<8; i++) out2[i] = mask[i] ? in1[i+4] : in2[i];
    //
    // replace unmasked in2 with shuffled in1
    const __mmask16 rmask_11110000 = _mm512_int2mask(0x00ff); //mask=00001111
    out2 = _mm512_mask_permute4f128_ps(in2, rmask_11110000, in1,_MM_PERM_BADC);
  }
  inline void trans2x8_new(double in1_[8], double in2_[8])
  {
    F64vec8& in1 = (F64vec8&)in1_[0];
    F64vec8& in2 = (F64vec8&)in2_[0];
    // copy the data that would be used after being rewritten:
    const F64vec8 buff1 = in1;
    trans2x8(in1, in2, buff1, in2);
  }
  inline void trans2x8(double in1[8], double in2[8])
  {
    __m512& data1 = (__m512&)in1[0];
    __m512& data2 = (__m512&)in2[0];
    // copy the data that we will first overwrite
    const __m512 buff1 = data1;
  
    // data1: copy high (odd) 1x4 element from low (even) 1x4 element of data2,
    //   i.e.,
    // for(int i=0; i<4; i++) data1[i+4] = data2[i];
    //
    // replace unmasked data1 with shuffled data2
    const __mmask16 rmask_00001111 = _mm512_int2mask(0xff00); //mask=11110000
    data1 = _mm512_mask_permute4f128_ps(data1, rmask_00001111, data2,_MM_PERM_BADC);
  
    // data2: copy low (even) 1x4 element from high (odd) 1x4 element of buff1,
    //   i.e.,
    // for(int i=0; i<4; i++) data2[i] = buff1[i+4];
    //
    // replace unmasked data2 with shuffled data1
    const __mmask16 rmask_11110000 = _mm512_int2mask(0x00ff); //mask=00001111
    data2 = _mm512_mask_permute4f128_ps(data2, rmask_11110000, buff1,_MM_PERM_BADC);
  }
  
  // transpose with blocked algorithm in
  // 24=8*3 = 8*log_2(8) vector instructions
  // (ignoring 12=4*3 copies made to
  // buffers to reduce number of registers needed)
  //
  inline void transpose_8x8_double_new(double in_[8][8])
  {
    ASSUME_ALIGNED(in_);
    F64vec8* data = (F64vec8*) in_[0];
    // using this auxiliary buffer should decrease the number of
    // instructions needed by 8 but should increase the number of
    // registers used by 7 compared to the implementation in
    // transpose_8x8_double_old
    F64vec8 buff[8];
    // 1. transpose each 2x2 block.
    for(int i=0; i<8; i+=2)
      trans2x2(buff[i], buff[i+1], data[i], data[i+1]);
    // 2. transpose each 4x4 block of 2x2 elements
    trans2x4(buff[0], buff[2]);
    trans2x4(buff[1], buff[3]);
    trans2x4(buff[4], buff[6]);
    trans2x4(buff[5], buff[7]);
    // 3. swap lower left and upper right 4x4 elements
    for(int i=0; i<4; i+=1)
      trans2x8(data[i], data[i+4], buff[i], buff[i+4]);
  }

  inline void transpose_8x8_double(double data[8][8])
  {
    ASSUME_ALIGNED(data);
    // 1. transpose each 2x2 block.
    for(int i=0; i<8; i+=2)
      trans2x2(data[i], data[i+1]);
    // 2. transpose each 4x4 block of 2x2 elements
    trans2x4(data[0], data[2]);
    trans2x4(data[1], data[3]);
    trans2x4(data[4], data[6]);
    trans2x4(data[5], data[7]);
    // 3. swap lower left and upper right 4x4 elements
    for(int i=0; i<4; i+=1)
      trans2x8(data[i], data[i+4]);
  }

  // transpose 8x8 data from in to out, leaving in unmodified
  inline void transpose_8x8_double(F64vec8* out[8], const F64vec8 in[8])
  {
    F64vec8 buff[8];
    // 1. transpose each 2x2 block.
    for(int i=0; i<8; i+=2)
      trans2x2(buff[i], buff[i+1], in[i], in[i+1]);
    // 2. transpose each 4x4 block of 2x2 elements
    trans2x4(buff[0], buff[2]);
    trans2x4(buff[1], buff[3]);
    trans2x4(buff[4], buff[6]);
    trans2x4(buff[5], buff[7]);
    // 3. swap lower left and upper right 4x4 elements
    for(int i=0; i<4; i+=1)
      trans2x8(*(out[i]), *(out[i+4]), buff[i], buff[i+4]);
  }

  // transpose 8x8 data from in to out, leaving in unmodified
  inline void transpose_8x8_double(F64vec8 out[8], F64vec8* in[8])
  {
    F64vec8 buff[8];
    // 1. transpose each 2x2 block.
    for(int i=0; i<8; i+=2)
      trans2x2(buff[i], buff[i+1], *(in[i]), *(in[i+1]));
    // 2. transpose each 4x4 block of 2x2 elements
    trans2x4(buff[0], buff[2]);
    trans2x4(buff[1], buff[3]);
    trans2x4(buff[4], buff[6]);
    trans2x4(buff[5], buff[7]);
    // 3. swap lower left and upper right 4x4 elements
    for(int i=0; i<4; i+=1)
      trans2x8(out[i], out[i+4], buff[i], buff[i+4]);
  }

  /*** end methods for fast transpose ***/

  typedef I32vec16 Ivec;
  typedef F64vec8 Dvec;

//#elif defined(__AVX__)
//  // see http://software.intel.com/en-us/node/462716
//  #include <dvec.h>
//  typedef I32vec8 Ivec;
//  typedef F64vec4 Dvec;
//#elif defined(__SSE2__)
#else
  class Ivec8
  {
    int a[8];
  };
  class Dvec4
  {
    double a[4];
   public:
    Dvec4(double x, double y, double z, double t)
    { a[0]=x; a[1]=y; a[2]=z; a[3]=t; }
    Dvec4(double s)
    { Dvec4(s,s,s,s); }
  };
  inline Dvec4 initialize_vector(double x, double y, double z, double t)
  {
    return Dvec4(t,z,y,x);
  }
  //inline Ivec8 round_down(Dvec4 dvec)
  //{
  //  Ivec8 ret;
  //  for(int i=0;i<4;i++)
  //    ret[i] = floor(dvec[i]);
  //}
  //typedef Ivec8 Ivec;
  typedef Dvec4 Dvec;
#endif

#endif
