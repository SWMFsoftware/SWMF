#ifndef mic_particles_h
#define mic_particles_h
#include "debug.h"
#if defined(__MIC__)
  #include "mic_basics.h"
  #include "arraysfwd.h"

  inline F64vec8 sample_field_mic(
    F64vec8 weights,
    F64vec8 field_components[8])
  {
    F64vec8 fields = F64vec8(0.);
    #pragma unroll
    for(int c=0; c<8; c++)
      fields += F64vec8(weights[c])*field_components[c];
    return fields;
  }

  inline F64vec8 reciprocal_F64vec8(F64vec8 arg)
  {
    // is there a faster way?
    return F64vec8(1.)/arg;
    // convert to single precision for faster reciprocal
    //
    //F32vec16 arg32 = _mm512_cvtpd_pslo(arg); 
    //F32vec16 arg32inv = F32vec16(1.)/arg32;
    //F64vec8 ret64 = _mm512_cvtpslo_pd(arg32inv);
    //return ret64;
  }
  
  // construct weights for two particles based on their canonical positions
  // 
  // each particle has 8 weights, one for each corner of its cell
  //
  inline void construct_weights_for_2pcls(F64vec8 weights[2], F64vec8 X0)
  {
    const F64vec8 X1 = F64vec8(1.) - X0;
    
    // construct weights for particles
    //
    // need versions with broadcasted component.
    //
    // doing via swizzle works only for the x component and
    // requires the same number of operations as the more general
    // technique, so for simplicity I do them all the same rather
    // than handling the x component with the commented-out lines.
    //const F64vec8 Xx0 = X0.aaaa();
    //const F64vec8 Xx1 = X1.aaaa();
    const __mmask8 rmask_01010101 = _mm512_int2mask(0xaa); // mask=10101010
    const __mmask8 rmask_00110011 = _mm512_int2mask(0xcc); // mask=11001100
    const __mmask8 rmask_00001111 = _mm512_int2mask(0xf0); // mask=11110000
    // construct weights for first particle
    //
    // weights for first particle use low vectors of X0 and X1:
    //
    //   weights[0][0] = X0[0]*X0[1]*X0[2]; // weight000
    //   weights[0][1] = X0[0]*X0[1]*X1[2]; // weight001
    //   weights[0][2] = X0[0]*X1[1]*X0[2]; // weight010
    //   weights[0][3] = X0[0]*X1[1]*X1[2]; // weight011
    //   weights[0][4] = X1[0]*X0[1]*X0[2]; // weight100
    //   weights[0][5] = X1[0]*X0[1]*X1[2]; // weight101
    //   weights[0][6] = X1[0]*X1[1]*X0[2]; // weight110
    //   weights[0][7] = X1[0]*X1[1]*X1[2]; // weight111
    //
    {
      //wx = cat_low_halves(Xx0,Xx1);
      const F64vec8 wx = broadcast_mask_blend(X0[0],X1[0],rmask_00001111);
      const F64vec8 wy = broadcast_mask_blend(X0[1],X1[1],rmask_00110011);
      const F64vec8 wz = broadcast_mask_blend(X0[2],X1[2],rmask_01010101);
      weights[0] = wx*wy*wz;
    }
    // construct weights for second particle
    //
    // weights for second particle use high vectors of X0 and X1:
    //
    //   weights[1][0] = X0[4]*X0[5]*X0[6]; // weight000
    //   weights[1][1] = X0[4]*X0[5]*X1[6]; // weight001
    //   weights[1][2] = X0[4]*X1[5]*X0[6]; // weight010
    //   weights[1][3] = X0[4]*X1[5]*X1[6]; // weight011
    //   weights[1][4] = X1[4]*X0[5]*X0[6]; // weight100
    //   weights[1][5] = X1[4]*X0[5]*X1[6]; // weight101
    //   weights[1][6] = X1[4]*X1[5]*X0[6]; // weight110
    //   weights[1][7] = X1[4]*X1[5]*X1[6]; // weight111
    //
    {
      //wx = cat_hgh_halves(Xx0,Xx1);
      const F64vec8 wx = broadcast_mask_blend(X0[4],X1[4],rmask_00001111);
      const F64vec8 wy = broadcast_mask_blend(X0[5],X1[5],rmask_00110011);
      const F64vec8 wz = broadcast_mask_blend(X0[6],X1[6],rmask_01010101);
      weights[1] = wx*wy*wz;
    }
  }
  
  inline F64vec8 compute_uvg_for_2pcls(F64vec8 uorig, F64vec8 B, F64vec8 E, F64vec8 qdto2mc)
  {
    /* F64vec8::aaaa() etc is (defectively) not declared const */
    /*const*/ F64vec8 Om = qdto2mc*B;
    const F64vec8 Omx = Om.aaaa();
    const F64vec8 Omy = Om.bbbb();
    const F64vec8 Omz = Om.cccc();
    // if pushing more than 2 particles at a time, could make
    // use of reduction intrinsics (_mm512_mask_reduce_add_pd)
    // and could calculate the reciprocal more efficiently.
    const F64vec8 omsq_p1 = F64vec8(1.) + Omx*Omx + Omy*Omy + Omz*Omz;
    // is there a faster way to implement the reciprocal?
    const F64vec8 denom = reciprocal_F64vec8(omsq_p1);
    // solve the position equation
    /*const*/ F64vec8 ut = uorig + qdto2mc*E;
    const F64vec8 udotOm =
        ut.aaaa()*Omx
      + ut.bbbb()*Omy
      + ut.cccc()*Omz;
  
    // solve the velocity equation 
    //
    return (ut + cross_product(ut,Om)+udotOm*Om)*denom;
  }

  inline void get_field_components_for_cell(
    F64vec8 field_components0[8],
    F64vec8 field_components1[8],
    const_arr4_double fieldForPcls,
    I32vec16 cx_)
  {
    // interface to the right of cell
    const I32vec16 ix_ = cx_ + I32vec16(1);
  
    int* cx = (int*)&cx_;
    int* ix = (int*)&ix_;
  
    // compute field_components array for particles
    //
    // compute the 1-dimensional subscripts for
    // the 8 corners of the mesh cell
    //
    {
      // broadcast coordinates
      const I32vec16 iX0 = I32vec16(ix[0]);
      const I32vec16 iY0 = I32vec16(ix[1]);
      const I32vec16 iZ0 = I32vec16(ix[2]);
      const I32vec16 cX0 = I32vec16(cx[0]);
      const I32vec16 cY0 = I32vec16(cx[1]);
      const I32vec16 cZ0 = I32vec16(cx[2]);
      const I32vec16 iX1 = I32vec16(ix[0+4]);
      const I32vec16 iY1 = I32vec16(ix[1+4]);
      const I32vec16 iZ1 = I32vec16(ix[2+4]);
      const I32vec16 cX1 = I32vec16(cx[0+4]);
      const I32vec16 cY1 = I32vec16(cx[1+4]);
      const I32vec16 cZ1 = I32vec16(cx[2+4]);
      // for each particle we want:
      //
      //   field_components[0] = fieldForPcls[iX][iY][iZ]; // field000
      //   field_components[1] = fieldForPcls[iX][iY][cZ]; // field001
      //   field_components[2] = fieldForPcls[iX][cY][iZ]; // field010
      //   field_components[3] = fieldForPcls[iX][cY][cZ]; // field011
      //   field_components[4] = fieldForPcls[cX][iY][iZ]; // field100
      //   field_components[5] = fieldForPcls[cX][iY][cZ]; // field101
      //   field_components[6] = fieldForPcls[cX][cY][iZ]; // field110
      //   field_components[7] = fieldForPcls[cX][cY][cZ]; // field111
      //
      // broadcast coordinates of first and second particle to low and high end of I32vec16
      const __mmask16 rmask_00ff = _mm512_int2mask(0xff00); // mask=1111111100000000
      const I32vec16 iX = _mm512_mask_blend_epi32(rmask_00ff,iX0,iX1);
      const I32vec16 iY = _mm512_mask_blend_epi32(rmask_00ff,iY0,iY1);
      const I32vec16 iZ = _mm512_mask_blend_epi32(rmask_00ff,iZ0,iZ1);
      const I32vec16 cX = _mm512_mask_blend_epi32(rmask_00ff,cX0,cX1);
      const I32vec16 cY = _mm512_mask_blend_epi32(rmask_00ff,cY0,cY1);
      const I32vec16 cZ = _mm512_mask_blend_epi32(rmask_00ff,cZ0,cZ1);
      // now use masks to form all possible combinations of coordinates
      const __mmask16 rmask_00001111x2 = _mm512_int2mask(0xf0f0); // mask=11110000x2
      const __mmask16 rmask_00110011x2 = _mm512_int2mask(0xcccc); // mask=11001100x2
      const __mmask16 rmask_01010101x2 = _mm512_int2mask(0xaaaa); // mask=10101010x2
      // x, y, and z subscripts of the 8 corners for both particles
      const I32vec16 sX = _mm512_mask_blend_epi32(rmask_00001111x2,iX,cX);
      const I32vec16 sY = _mm512_mask_blend_epi32(rmask_00110011x2,iY,cY);
      const I32vec16 sZ = _mm512_mask_blend_epi32(rmask_01010101x2,iZ,cZ);
      // compute the starting 1-dimensional index for the 8 corners
      const int nxn = fieldForPcls.dim1();
      const int nyn = fieldForPcls.dim2();
      const int nzn = fieldForPcls.dim3();
      const int nnn = fieldForPcls.dim4();
      const I32vec16 nXn(nxn);
      const I32vec16 nYn(nyn);
      const I32vec16 nZn(nzn);
      // compiler could see whether this would fail, so should cost nothing
      assert(nnn==8); // needed for field_components to be aligned
      const I32vec16 nNn(nnn);
      // unfortunately I32vec16 does not support this syntax,
      //const MyI32vec16 subs = ((sX*nYn+sY)*nZn+sZ)*nNn;
      const I32vec16 s1 = fused_multiply_add(sX,nYn,sY);
      const I32vec16 s2 = fused_multiply_add(s1,nZn,sZ);
      const I32vec16 subs = multiply(s2,nNn);
      //printexpr(subs);
      int const*const subs_arr = (int*)&subs;
      // access underlying 1D array of fieldForPcls
      const double* fieldForPcls1d = fieldForPcls.get_arr();
      #pragma unroll
      for(int i=0; i<8; i++)
      {
        field_components0[i] = *(F64vec8*)&fieldForPcls1d[subs_arr[i]];
      }
      #pragma unroll
      for(int i=0; i<8; i++)
      {
        field_components1[i] = *(F64vec8*)&fieldForPcls1d[subs_arr[i+8]];
      }
    }
  }

#endif
#endif
