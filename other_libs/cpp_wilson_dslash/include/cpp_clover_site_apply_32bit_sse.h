#ifndef CPP_CLOVER_SITE_APPLY_32_BIT_SSE_H
#define CPP_CLOVER_SITE_APPLY_32_BIT_SSE_H

#include <cpp_clover_types.h>
#include <cpp_dslash_types.h>

#include <xmmintrin.h>


namespace CPlusPlusClover { 

  // Use 32bit types
  using namespace Clover32BitTypes;
using namespace Dslash32BitTypes;

  namespace CPlusPlusClover32Bit { 


    inline 
    void cloverSiteApply(FourSpinor result, const CloverTerm clov, const FourSpinor src) 
    {

      // Unrolled version 3. 
      // Took unrolled version 2 and wrote out in real() and imag() 
      // parts. Rearranged so that all the reals follow each other
      //  in the output so that we can write linearly
      __m128 r;
      __m128 s;
      typedef union { 
	float floats[4];
	__m128 vec;
      } V4F;

      __m128 c;
      
      V4F mask1111 __attribute__((unused)) = { 1.0, 1.0, 1.0, 1.0 };
      V4F mask1001 __attribute__((unused)) = { 1.0, -1.0,-1.0, 1.0 };
      V4F mask1010 __attribute__((unused)) = { 1.0, -1.0,1.0,-1.0 };
      V4F mask0101 __attribute__((unused)) = { -1.0, 1.0,-1.0, 1.0 };
      __m128 tmp, tmp1;

      // result[0][0][0]  = clov[0].diag[0] * src[0][0][0];
      // result[0][0][1]  = clov[0].diag[0] * src[0][0][1];
      // result[0][1][0]  = clov[0].diag[1] * src[0][1][0];
      // result[0][1][1]  = clov[0].diag[1] * src[0][1][1];
      s = _mm_load_ps((float *)&(src[0][0][0]));

      //  c 
      //   diag[0] diag[1] junk junk
      c = _mm_loadl_pi(tmp, (__m64 *)&(clov[0].diag[0]));
      //  tmp 
      //   junk junk diag[0] diag[1]
      tmp = _mm_movelh_ps(tmp,c);      

      // Want c = diag[0] diag[0] diag[1] diag[1] 
      c   = _mm_shuffle_ps(c, tmp, 0xf0);
      r = _mm_mul_ps(c,s);

      //result[0][0][0] += clov[0].off_diag[0][0] * src[0][1][0];
      //result[0][0][1] += clov[0].off_diag[0][0] * src[0][1][1];
      //result[0][1][0] += clov[0].off_diag[0][0] * src[0][0][0];
      //result[0][1][1] += clov[0].off_diag[0][0] * src[0][0][1];
      s = _mm_loadl_pi(s, (__m64 *)&(src[0][1][0]));
      s = _mm_loadh_pi(s, (__m64 *)&(src[0][0][0]));
      c = _mm_load_ps1((float *)&(clov[0].off_diag[0][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[0][0][0] += clov[0].off_diag[0][1] * src[0][1][1];
      // result[0][0][1] -= clov[0].off_diag[0][1] * src[0][1][0];
      // result[0][1][0] -= clov[0].off_diag[0][1] * src[0][0][1];
      // result[0][1][1] += clov[0].off_diag[0][1] * src[0][0][0];   
      c = _mm_load_ps1((float *)&(clov[0].off_diag[0][1]));
      s = _mm_shuffle_ps(s,s,0xb1); 
      c = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1001.vec, c);
      r = _mm_add_ps(r, tmp);

      // result[0][0][0] += clov[0].off_diag[1][0] * src[0][2][0];
      // result[0][0][1] -= clov[0].off_diag[1][1] * src[0][2][0];
      // result[0][1][0] += clov[0].off_diag[2][0] * src[0][2][0];
      // result[0][1][1] -= clov[0].off_diag[2][1] * src[0][2][0];
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[1][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[2][0]));
      s = _mm_load_ps1((float *)&(src[0][2][0]));
      tmp1 = _mm_mul_ps(c, s);
      tmp = _mm_mul_ps(mask1010.vec,tmp1);
      r = _mm_add_ps(r,tmp);

      // result[0][0][0] += clov[0].off_diag[1][1] * src[0][2][1];
      // result[0][0][1] += clov[0].off_diag[1][0] * src[0][2][1];
      // result[0][1][0] += clov[0].off_diag[2][1] * src[0][2][1];
      // result[0][1][1] += clov[0].off_diag[2][0] * src[0][2][1];
      s = _mm_load_ps1((float *)&(src[0][2][1]));
      c = _mm_shuffle_ps(c,c,0xb1); 
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[0][0][0] += clov[0].off_diag[3][0] * src[1][0][0];
      // result[0][0][1] -= clov[0].off_diag[3][1] * src[1][0][0];
      // result[0][1][0] += clov[0].off_diag[4][0] * src[1][0][0];
      // result[0][1][1] -= clov[0].off_diag[4][1] * src[1][0][0];
      s = _mm_load_ps1((float *)&(src[1][0][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[3]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[4]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec, tmp1);
      r = _mm_add_ps(r,tmp);

      // result[0][0][0] += clov[0].off_diag[3][1] * src[1][0][1];
      // result[0][0][1] += clov[0].off_diag[3][0] * src[1][0][1];
      // result[0][1][0] += clov[0].off_diag[4][1] * src[1][0][1];
      // result[0][1][1] += clov[0].off_diag[4][0] * src[1][0][1];
      s = _mm_load_ps1((float *)&(src[1][0][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[0][0][0] += clov[0].off_diag[6][0] * src[1][1][0];
      // result[0][0][1] -= clov[0].off_diag[6][1] * src[1][1][0];
      // result[0][1][0] += clov[0].off_diag[7][0] * src[1][1][0];
      // result[0][1][1] -= clov[0].off_diag[7][1] * src[1][1][0];
      s = _mm_load_ps1((float *)&(src[1][1][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[6]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[7]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec, tmp1);
      r = _mm_add_ps(r,tmp);

      // result[0][0][0] += clov[0].off_diag[6][1] * src[1][1][1];
      // result[0][0][1] += clov[0].off_diag[6][0] * src[1][1][1];
      // result[0][1][0] += clov[0].off_diag[7][1] * src[1][1][1];
      // result[0][1][1] += clov[0].off_diag[7][0] * src[1][1][1];
      s = _mm_load_ps1((float *)&(src[1][1][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);


      // result[0][0][0] += clov[0].off_diag[10][0] * src[1][2][0];
      // result[0][0][1] -= clov[0].off_diag[10][1] * src[1][2][0];
      // result[0][1][0] += clov[0].off_diag[11][0] * src[1][2][0];
      // result[0][1][1] -= clov[0].off_diag[11][1] * src[1][2][0];

      s = _mm_load_ps1((float *)&(src[1][2][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[10]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[11]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec, tmp1);
      r = _mm_add_ps(r,tmp);


      // result[0][0][0] += clov[0].off_diag[10][1] * src[1][2][1];
      // result[0][0][1] += clov[0].off_diag[10][0] * src[1][2][1];
      // result[0][1][0] += clov[0].off_diag[11][1] * src[1][2][1];
      // result[0][1][1] += clov[0].off_diag[11][0] * src[1][2][1];
      
      s = _mm_load_ps1((float *)&(src[1][2][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);
      _mm_store_ps((float *)&(result[0][0][0]), r);

      // result[0][2][0] = clov[0].diag[2] * src[0][2][0];
      // result[0][2][1] = clov[0].diag[2] * src[0][2][1];
      // result[1][0][0] = clov[0].diag[3] * src[1][0][0];
      // result[1][0][1] = clov[0].diag[3] * src[1][0][1];
      s = _mm_load_ps((float *)&(src[0][2][0]));
      
      // c:
      //   diag[2] diag[3] junk junk
      c = _mm_loadl_pi(tmp, (__m64 *)&(clov[0].diag[2]));
      //  tmp 
      //   junk junk diag[2] diag[3]
      tmp = _mm_movelh_ps(tmp,c);      

      // Want c = diag[2] diag[2] diag[3] diag[3] 
      c   = _mm_shuffle_ps(c, tmp, 0xf0);
      r = _mm_mul_ps(c,s);


      //result[0][2][0] += clov[0].off_diag[ 1][0] * src[0][0][0];
      //result[0][2][1] += clov[0].off_diag[ 1][1] * src[0][0][0];
      //result[1][0][0] += clov[0].off_diag[ 3][0] * src[0][0][0];
      //result[1][0][1] += clov[0].off_diag[ 3][1] * src[0][0][0];
      s = _mm_load_ps1((float *)&(src[0][0][0]));
      c = _mm_loadl_pi(c, (__m64 *)&(clov[0].off_diag[1][0]));
      c = _mm_loadh_pi(c, (__m64 *)&(clov[0].off_diag[3][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);	

      // result[0][2][0] -= clov[0].off_diag[ 1][1] * src[0][0][1];
      // result[0][2][1] += clov[0].off_diag[ 1][0] * src[0][0][1];
      // result[1][0][0] -= clov[0].off_diag[ 3][1] * src[0][0][1];
      // result[1][0][1] += clov[0].off_diag[ 3][0] * src[0][0][1];
      s = _mm_load_ps1((float *)&(src[0][0][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);

      // result[0][2][0] += clov[0].off_diag[ 2][0] * src[0][1][0];
      // result[0][2][1] += clov[0].off_diag[ 2][1] * src[0][1][0];
      // result[1][0][0] += clov[0].off_diag[ 4][0] * src[0][1][0];
      // result[1][0][1] += clov[0].off_diag[ 4][1] * src[0][1][0];
      s = _mm_load_ps1((float *)&(src[0][1][0]));
      c = _mm_loadl_pi(c, (__m64 *)&(clov[0].off_diag[2][0]));
      c = _mm_loadh_pi(c, (__m64 *)&(clov[0].off_diag[4][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);	

      //result[0][2][0] -= clov[0].off_diag[ 2][1] * src[0][1][1];
      //result[0][2][1] += clov[0].off_diag[ 2][0] * src[0][1][1];
      //result[1][0][0] -= clov[0].off_diag[ 4][1] * src[0][1][1];
      //result[1][0][1] += clov[0].off_diag[ 4][0] * src[0][1][1];
      s = _mm_load_ps1((float *)&(src[0][1][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);


      // result[0][2][0] += clov[0].off_diag[5][0] * src[1][0][0];
      // result[0][2][1] += clov[0].off_diag[5][0] * src[1][0][1];
      // result[1][0][0] += clov[0].off_diag[5][0] * src[0][2][0];
      // result[1][0][1] += clov[0].off_diag[5][0] * src[0][2][1];
      c = _mm_load_ps1((float *)&(clov[0].off_diag[5][0]));
      s = _mm_loadl_pi(s,(__m64 *)&(src[1][0][0])); 
      s = _mm_loadh_pi(s,(__m64 *)&(src[0][2][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[0][2][0] += clov[0].off_diag[5][1] * src[1][0][1];
      // result[0][2][1] -= clov[0].off_diag[5][1] * src[1][0][0];
      // result[1][0][0] -= clov[0].off_diag[5][1] * src[0][2][1];
      // result[1][0][1] += clov[0].off_diag[5][1] * src[0][2][0];
      c = _mm_load_ps1((float *)&(clov[0].off_diag[5][1]));
      s = _mm_shuffle_ps(s,s,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1001.vec, tmp1);
      r = _mm_add_ps(r,tmp);

      // result[0][2][0] += clov[0].off_diag[8][0]  * src[1][1][0];
      // result[0][2][1] -= clov[0].off_diag[8][1]  * src[1][1][0];
      // result[1][0][0] += clov[0].off_diag[ 9][0] * src[1][1][0];
      // result[1][0][1] -= clov[0].off_diag[ 9][1] * src[1][1][0];

      s = _mm_load_ps1((float *)&(src[1][1][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[8][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[9][0]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec,tmp1);
      r = _mm_add_ps(r,tmp);

      // result[0][2][0] += clov[0].off_diag[8][1]  * src[1][1][1];
      // result[0][2][1] += clov[0].off_diag[8][0]  * src[1][1][1];
      // result[1][0][0] += clov[0].off_diag[ 9][1] * src[1][1][1];
      // result[1][0][1] += clov[0].off_diag[ 9][0] * src[1][1][1];
      s = _mm_load_ps1((float *)&(src[1][1][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r, tmp);


      // result[0][2][0] += clov[0].off_diag[12][0] * src[1][2][0];
      // result[0][2][1] -= clov[0].off_diag[12][1] * src[1][2][0];
      // result[1][0][0] += clov[0].off_diag[13][0] * src[1][2][0];
      // result[1][0][1] -= clov[0].off_diag[13][1] * src[1][2][0];
      s = _mm_load_ps1((float *)&(src[1][2][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[12][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[13][0]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec,tmp1);
      r = _mm_add_ps(r,tmp);

      // result[0][2][0] += clov[0].off_diag[12][1] * src[1][2][1];
      // result[0][2][1] += clov[0].off_diag[12][0] * src[1][2][1];
      // result[1][0][0] += clov[0].off_diag[13][1] * src[1][2][1];
      // result[1][0][1] += clov[0].off_diag[13][0] * src[1][2][1];

      s = _mm_load_ps1((float *)&(src[1][2][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r, tmp);

      _mm_store_ps((float *)&(result[0][2][0]),r);	       

      //      result[1][1][0]  = clov[0].diag[ 4] * src[1][1][0];
      //      result[1][1][1]  = clov[0].diag[ 4] * src[1][1][1];
      //      result[1][2][0]  = clov[0].diag[ 5] * src[1][2][0];
      //      result[1][2][1]  = clov[0].diag[ 5] * src[1][2][1];
      s = _mm_load_ps((float *)&(src[1][1][0]));

      // c:
      //   diag[4] diag[5] junk junk
      c = _mm_loadl_pi(tmp, (__m64 *)&(clov[0].diag[4]));
      //  tmp 
      //   junk junk diag[4] diag[5]
      tmp = _mm_movelh_ps(tmp,c);      

      // Want c = diag[4] diag[4] diag[5] diag[5] 
      c   = _mm_shuffle_ps(c, tmp, 0xf0);
      r = _mm_mul_ps(c,s);
      


      // result[1][1][0] += clov[0].off_diag[ 6][0] * src[0][0][0];
      // result[1][1][1] += clov[0].off_diag[ 6][1] * src[0][0][0];
      // result[1][2][0] += clov[0].off_diag[10][0] * src[0][0][0];
      // result[1][2][1] += clov[0].off_diag[10][1] * src[0][0][0];
      s = _mm_load_ps1((float *)&(src[0][0][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[6][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[10][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[1][1][0] -= clov[0].off_diag[ 6][1] * src[0][0][1];
      // result[1][1][1] += clov[0].off_diag[ 6][0] * src[0][0][1];
      // result[1][2][0] -= clov[0].off_diag[10][1] * src[0][0][1];
      // result[1][2][1] += clov[0].off_diag[10][0] * src[0][0][1];
      s = _mm_load_ps1((float *)&(src[0][0][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);

      // result[1][1][0] += clov[0].off_diag[ 7][0] * src[0][1][0];
      // result[1][1][1] += clov[0].off_diag[ 7][1] * src[0][1][0];
      // result[1][2][0] += clov[0].off_diag[11][0] * src[0][1][0];
      // result[1][2][1] += clov[0].off_diag[11][1] * src[0][1][0];
      s = _mm_load_ps1((float *)&(src[0][1][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[7][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[11][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[1][1][0] -= clov[0].off_diag[ 7][1] * src[0][1][1];
      // result[1][1][1] += clov[0].off_diag[ 7][0] * src[0][1][1];
      // result[1][2][0] -= clov[0].off_diag[11][1] * src[0][1][1];
      // result[1][2][1] += clov[0].off_diag[11][0] * src[0][1][1];
      s = _mm_load_ps1((float *)&(src[0][1][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);



      //result[1][1][0] += clov[0].off_diag[ 8][0] * src[0][2][0];
      //result[1][1][1] += clov[0].off_diag[ 8][1] * src[0][2][0];
      //result[1][2][0] += clov[0].off_diag[12][0] * src[0][2][0];
      //result[1][2][1] += clov[0].off_diag[12][1] * src[0][2][0];
      s = _mm_load_ps1((float *)&(src[0][2][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[8][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[12][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);


      // result[1][1][0] -= clov[0].off_diag[ 8][1] * src[0][2][1];
      // result[1][1][1] += clov[0].off_diag[ 8][0] * src[0][2][1];
      // result[1][2][0] -= clov[0].off_diag[12][1] * src[0][2][1];
      // result[1][2][1] += clov[0].off_diag[12][0] * src[0][2][1];
      s = _mm_load_ps1((float *)&(src[0][2][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);

      // result[1][1][0] += clov[0].off_diag[ 9][0] * src[1][0][0];
      // result[1][1][1] += clov[0].off_diag[ 9][1] * src[1][0][0];
      // result[1][2][0] += clov[0].off_diag[13][0] * src[1][0][0];
      // result[1][2][1] += clov[0].off_diag[13][1] * src[1][0][0];
      s = _mm_load_ps1((float *)&(src[1][0][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[0].off_diag[9][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[0].off_diag[13][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);


      // result[1][1][0] -= clov[0].off_diag[ 9][1] * src[1][0][1];
      // result[1][1][1] += clov[0].off_diag[ 9][0] * src[1][0][1];
      // result[1][2][0] -= clov[0].off_diag[13][1] * src[1][0][1];
      // result[1][2][1] += clov[0].off_diag[13][0] * src[1][0][1];
      s = _mm_load_ps1((float *)&(src[1][0][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);

      // result[1][1][0] += clov[0].off_diag[14][0] * src[1][2][0];
      // result[1][1][1] += clov[0].off_diag[14][0] * src[1][2][1];
      // result[1][2][0] += clov[0].off_diag[14][0] * src[1][1][0];
      // result[1][2][1] += clov[0].off_diag[14][0] * src[1][1][1];
      c = _mm_load_ps1((float *)&(clov[0].off_diag[14][0]));
      s = _mm_loadl_pi(s,(__m64 *)&(src[1][2][0]));
      s = _mm_loadh_pi(s,(__m64 *)&(src[1][1][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);
      
      // result[1][1][0] += clov[0].off_diag[14][1] * src[1][2][1];
      // result[1][1][1] -= clov[0].off_diag[14][1] * src[1][2][0];
      // result[1][2][0] -= clov[0].off_diag[14][1] * src[1][1][1];
      // result[1][2][1] += clov[0].off_diag[14][1] * src[1][1][0];
      c = _mm_load_ps1((float *)&(clov[0].off_diag[14][1]));
      s = _mm_shuffle_ps(s,s, 0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1001.vec, tmp1);
      r = _mm_add_ps(r,tmp);
     
      _mm_store_ps((float *)&(result[1][1][0]),r);	             


      // result[2][0][0]  = clov[1].diag[0]  * src[2][0][0];
      // result[2][0][1]  = clov[1].diag[0]  * src[2][0][1];
      // result[2][1][0]  = clov[1].diag[1]  * src[2][1][0];
      // result[2][1][1]  = clov[1].diag[1]  * src[2][1][1];
      s = _mm_load_ps((float *)&(src[2][0][0]));
           //  c 
      //   diag[0] diag[1] junk junk
      c = _mm_loadl_pi(tmp, (__m64 *)&(clov[1].diag[0]));
      //  tmp 
      //   junk junk diag[0] diag[1]
      tmp = _mm_movelh_ps(tmp,c);      

      // Want c = diag[0] diag[0] diag[1] diag[1] 
      c   = _mm_shuffle_ps(c, tmp, 0xf0);
      r = _mm_mul_ps(c,s);



      // result[2][0][0] += clov[1].off_diag[0][0]  * src[2][1][0];
      // result[2][0][1] += clov[1].off_diag[0][0]  * src[2][1][1];
      // result[2][1][0] += clov[1].off_diag[0][0]  * src[2][0][0];
      // result[2][1][1] += clov[1].off_diag[0][0]  * src[2][0][1];

      s = _mm_loadl_pi(s, (__m64 *)&(src[2][1][0]));
      s = _mm_loadh_pi(s, (__m64 *)&(src[2][0][0]));
      c = _mm_load_ps1((float *)&(clov[1].off_diag[0][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[2][0][0] += clov[1].off_diag[0][1]  * src[2][1][1];
      // result[2][0][1] -= clov[1].off_diag[0][1]  * src[2][1][0];
      // result[2][1][0] -= clov[1].off_diag[0][1]  * src[2][0][1];
      // result[2][1][1] += clov[1].off_diag[0][1]  * src[2][0][0];
      c = _mm_load_ps1((float *)&(clov[1].off_diag[0][1]));
      s = _mm_shuffle_ps(s,s,0xb1); 
      c = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1001.vec, c);
      r = _mm_add_ps(r, tmp);

      // result[2][0][0] += clov[1].off_diag[1][0]  * src[2][2][0];
      // result[2][0][1] -= clov[1].off_diag[1][1]  * src[2][2][0];
      // result[2][1][0] += clov[1].off_diag[2][0]  * src[2][2][0];
      // result[2][1][1] -= clov[1].off_diag[2][1]  * src[2][2][0];
      c = _mm_load_ps((float *)&(clov[1].off_diag[1][0]));
      s = _mm_load_ps1((float *)&(src[2][2][0]));
      tmp1 = _mm_mul_ps(c, s);
      tmp = _mm_mul_ps(mask1010.vec,tmp1);
      r = _mm_add_ps(r,tmp);

      // result[2][0][0] += clov[1].off_diag[1][1]  * src[2][2][1];
      // result[2][0][1] += clov[1].off_diag[1][0]  * src[2][2][1];
      // result[2][1][0] += clov[1].off_diag[2][1]  * src[2][2][1];
      // result[2][1][1] += clov[1].off_diag[2][0]  * src[2][2][1];
      s = _mm_load_ps1((float *)&(src[2][2][1]));
      c = _mm_shuffle_ps(c,c,0xb1); 
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);


      // result[2][0][0] += clov[1].off_diag[3][0]  * src[3][0][0];
      // result[2][0][1] -= clov[1].off_diag[3][1]  * src[3][0][0];
      // result[2][1][0] += clov[1].off_diag[4][0]  * src[3][0][0];
      // result[2][1][1] -= clov[1].off_diag[4][1]  * src[3][0][0];
      s = _mm_load_ps1((float *)&(src[3][0][0]));
      c = _mm_load_ps((float *)&(clov[1].off_diag[3][0]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec, tmp1);
      r = _mm_add_ps(r,tmp);

      // result[2][0][0] += clov[1].off_diag[3][1]  * src[3][0][1];
      // result[2][0][1] += clov[1].off_diag[3][0]  * src[3][0][1];
      // result[2][1][0] += clov[1].off_diag[4][1]  * src[3][0][1];
      // result[2][1][1] += clov[1].off_diag[4][0]  * src[3][0][1];
      s = _mm_load_ps1((float *)&(src[3][0][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[2][0][0] += clov[1].off_diag[6][0]  * src[3][1][0];
      // result[2][0][1] -= clov[1].off_diag[6][1]  * src[3][1][0];
      // result[2][1][0] += clov[1].off_diag[7][0]  * src[3][1][0];
      // result[2][1][1] -= clov[1].off_diag[7][1]  * src[3][1][0];
      s = _mm_load_ps1((float *)&(src[3][1][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[1].off_diag[6]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[1].off_diag[7]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec, tmp1);
      r = _mm_add_ps(r,tmp);

      // result[2][0][0] += clov[1].off_diag[6][1]  * src[3][1][1];
      // result[2][0][1] += clov[1].off_diag[6][0]  * src[3][1][1];
      // result[2][1][0] += clov[1].off_diag[7][1]  * src[3][1][1];
      // result[2][1][1] += clov[1].off_diag[7][0]  * src[3][1][1];
      s = _mm_load_ps1((float *)&(src[3][1][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[2][0][0] += clov[1].off_diag[10][0] * src[3][2][0];
      // result[2][0][1] -= clov[1].off_diag[10][1] * src[3][2][0];
      // result[2][1][0] += clov[1].off_diag[11][0] * src[3][2][0];
      // result[2][1][1] -= clov[1].off_diag[11][1] * src[3][2][0];
      s = _mm_load_ps1((float *)&(src[3][2][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[1].off_diag[10][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[1].off_diag[11][0]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec, tmp1);
      r = _mm_add_ps(r,tmp);

      // result[2][0][0] += clov[1].off_diag[10][1] * src[3][2][1];
      // result[2][0][1] += clov[1].off_diag[10][0] * src[3][2][1];
      // result[2][1][0] += clov[1].off_diag[11][1] * src[3][2][1];
      // result[2][1][1] += clov[1].off_diag[11][0] * src[3][2][1];
       s = _mm_load_ps1((float *)&(src[3][2][1]));
       c = _mm_shuffle_ps(c,c,0xb1);
       tmp = _mm_mul_ps(c,s);
       r = _mm_add_ps(r,tmp);
      _mm_store_ps((float *)&(result[2][0][0]), r);

      // result[2][2][0]  = clov[1].diag[ 2]  * src[2][2][0];
      // result[2][2][1]  = clov[1].diag[ 2]  * src[2][2][1];
      // result[3][0][0]  = clov[1].diag[ 3]  * src[3][0][0];
      // result[3][0][1]  = clov[1].diag[ 3]  * src[3][0][1];
      s = _mm_load_ps((float *)&(src[2][2][0]));
      
      // c:
      //   diag[2] diag[3] junk junk
      c = _mm_loadl_pi(tmp, (__m64 *)&(clov[1].diag[2]));
      //  tmp 
      //   junk junk diag[2] diag[3]
      tmp = _mm_movelh_ps(tmp,c);      

      // Want c = diag[2] diag[2] diag[3] diag[3] 
      c   = _mm_shuffle_ps(c, tmp, 0xf0);
      r = _mm_mul_ps(c,s);

      // result[2][2][0] += clov[1].off_diag[ 1][0]  * src[2][0][0];
      // result[2][2][1] += clov[1].off_diag[ 1][1]  * src[2][0][0];
      // result[3][0][0] += clov[1].off_diag[ 3][0]  * src[2][0][0];
      // result[3][0][1] += clov[1].off_diag[ 3][1]  * src[2][0][0];
      s = _mm_load_ps1((float *)&(src[2][0][0]));
      c = _mm_loadl_pi(c, (__m64 *)&(clov[1].off_diag[1][0]));
      c = _mm_loadh_pi(c, (__m64 *)&(clov[1].off_diag[3][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);	

      // result[2][2][0] -= clov[1].off_diag[ 1][1]  * src[2][0][1];
      // result[2][2][1] += clov[1].off_diag[ 1][0]  * src[2][0][1];
      // result[3][0][0] -= clov[1].off_diag[ 3][1]  * src[2][0][1];
      // result[3][0][1] += clov[1].off_diag[ 3][0]  * src[2][0][1];
      s = _mm_load_ps1((float *)&(src[2][0][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);

      // result[2][2][0] += clov[1].off_diag[ 2][0]  * src[2][1][0];
      // result[2][2][1] += clov[1].off_diag[ 2][1]  * src[2][1][0];
      // result[3][0][0] += clov[1].off_diag[ 4][0]  * src[2][1][0];
      // result[3][0][1] += clov[1].off_diag[ 4][1]  * src[2][1][0];
      s = _mm_load_ps1((float *)&(src[2][1][0]));
      c = _mm_loadl_pi(c, (__m64 *)&(clov[1].off_diag[2][0]));
      c = _mm_loadh_pi(c, (__m64 *)&(clov[1].off_diag[4][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);	

      // result[2][2][0] -= clov[1].off_diag[ 2][1]  * src[2][1][1];
      // result[2][2][1] += clov[1].off_diag[ 2][0]  * src[2][1][1];
      // result[3][0][0] -= clov[1].off_diag[ 4][1]  * src[2][1][1];
      // result[3][0][1] += clov[1].off_diag[ 4][0]  * src[2][1][1];
      s = _mm_load_ps1((float *)&(src[2][1][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);


      // result[2][2][0] += clov[1].off_diag[5][0]   * src[3][0][0];
      // result[2][2][1] += clov[1].off_diag[5][0]   * src[3][0][1];
      // result[3][0][0] += clov[1].off_diag[ 5][0]  * src[2][2][0];
      // result[3][0][1] += clov[1].off_diag[ 5][0]  * src[2][2][1];
      c = _mm_load_ps1((float *)&(clov[1].off_diag[5][0]));
      s = _mm_loadl_pi(s,(__m64 *)&(src[3][0][0])); 
      s = _mm_loadh_pi(s,(__m64 *)&(src[2][2][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[2][2][0] += clov[1].off_diag[5][1]   * src[3][0][1];
      // result[2][2][1] -= clov[1].off_diag[5][1]   * src[3][0][0];
      // result[3][0][0] -= clov[1].off_diag[ 5][1]  * src[2][2][1];
      // result[3][0][1] += clov[1].off_diag[ 5][1]  * src[2][2][0];
      c = _mm_load_ps1((float *)&(clov[1].off_diag[5][1]));
      s = _mm_shuffle_ps(s,s,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1001.vec, tmp1);
      r = _mm_add_ps(r,tmp);


      // result[2][2][0] += clov[1].off_diag[8][0]   * src[3][1][0];
      // result[2][2][1] -= clov[1].off_diag[8][1]   * src[3][1][0];
      // result[3][0][0] += clov[1].off_diag[ 9][0]  * src[3][1][0];
      // result[3][0][1] -= clov[1].off_diag[ 9][1]  * src[3][1][0];
      s = _mm_load_ps1((float *)&(src[3][1][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[1].off_diag[8][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[1].off_diag[9][0]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec,tmp1);
      r = _mm_add_ps(r,tmp);

      // result[2][2][0] += clov[1].off_diag[8][1]   * src[3][1][1];
      // result[2][2][1] += clov[1].off_diag[8][0]   * src[3][1][1];
      // result[3][0][0] += clov[1].off_diag[ 9][1]  * src[3][1][1];
      // result[3][0][1] += clov[1].off_diag[ 9][0]  * src[3][1][1];
      s = _mm_load_ps1((float *)&(src[3][1][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r, tmp);

      // result[2][2][0] += clov[1].off_diag[12][0]  * src[3][2][0];
      // result[2][2][1] -= clov[1].off_diag[12][1]  * src[3][2][0];
      // result[3][0][0] += clov[1].off_diag[13][0]  * src[3][2][0];
      // result[3][0][1] -= clov[1].off_diag[13][1]  * src[3][2][0];
      s = _mm_load_ps1((float *)&(src[3][2][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[1].off_diag[12][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[1].off_diag[13][0]));
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1010.vec,tmp1);
      r = _mm_add_ps(r,tmp);

      // result[2][2][0] += clov[1].off_diag[12][1]  * src[3][2][1];
      // result[2][2][1] += clov[1].off_diag[12][0]  * src[3][2][1];
      // result[3][0][0] += clov[1].off_diag[13][1]  * src[3][2][1];
      // result[3][0][1] += clov[1].off_diag[13][0]  * src[3][2][1];
      s = _mm_load_ps1((float *)&(src[3][2][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r, tmp);
      _mm_store_ps((float *)&(result[2][2][0]),r);

      // result[3][1][0]  = clov[1].diag[ 4]  * src[3][1][0];
      // result[3][1][1]  = clov[1].diag[ 4]  * src[3][1][1];
      // result[3][2][0]  = clov[1].diag[ 5]  * src[3][2][0];
      // result[3][2][1]  = clov[1].diag[ 5]  * src[3][2][1];
      s = _mm_load_ps((float *)&(src[3][1][0]));

      // c:
      //   diag[4] diag[5] junk junk
      c = _mm_loadl_pi(tmp, (__m64 *)&(clov[1].diag[4]));
      //  tmp 
      //   junk junk diag[4] diag[5]
      tmp = _mm_movelh_ps(tmp,c);      

      // Want c = diag[4] diag[4] diag[5] diag[5] 
      c   = _mm_shuffle_ps(c, tmp, 0xf0);
      r = _mm_mul_ps(c,s);

      // result[3][1][0] += clov[1].off_diag[ 6][0]  * src[2][0][0];
      // result[3][1][1] += clov[1].off_diag[ 6][1]  * src[2][0][0];
      // result[3][2][0] += clov[1].off_diag[10][0]  * src[2][0][0];
      // result[3][2][1] += clov[1].off_diag[10][1]  * src[2][0][0];
      s = _mm_load_ps1((float *)&(src[2][0][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[1].off_diag[6][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[1].off_diag[10][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[3][1][0] -= clov[1].off_diag[ 6][1]  * src[2][0][1];
      // result[3][1][1] += clov[1].off_diag[ 6][0]  * src[2][0][1];
      // result[3][2][0] -= clov[1].off_diag[10][1]  * src[2][0][1];
      // result[3][2][1] += clov[1].off_diag[10][0]  * src[2][0][1];
      s = _mm_load_ps1((float *)&(src[2][0][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);

      // result[3][1][0] += clov[1].off_diag[ 7][0]  * src[2][1][0];
      // result[3][1][1] += clov[1].off_diag[ 7][1]  * src[2][1][0];
      // result[3][2][0] += clov[1].off_diag[11][0]  * src[2][1][0];
      // result[3][2][1] += clov[1].off_diag[11][1]  * src[2][1][0];
      s = _mm_load_ps1((float *)&(src[2][1][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[1].off_diag[7][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[1].off_diag[11][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[3][1][0] -= clov[1].off_diag[ 7][1]  * src[2][1][1];
      // result[3][1][1] += clov[1].off_diag[ 7][0]  * src[2][1][1];
      // result[3][2][0] -= clov[1].off_diag[11][1]  * src[2][1][1];
      // result[3][2][1] += clov[1].off_diag[11][0]  * src[2][1][1];
      s = _mm_load_ps1((float *)&(src[2][1][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);

      // result[3][1][0] += clov[1].off_diag[ 8][0]  * src[2][2][0];
      // result[3][1][1] += clov[1].off_diag[ 8][1]  * src[2][2][0];
      // result[3][2][0] += clov[1].off_diag[12][0]  * src[2][2][0];
      // result[3][2][1] += clov[1].off_diag[12][1]  * src[2][2][0];
      s = _mm_load_ps1((float *)&(src[2][2][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[1].off_diag[8][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[1].off_diag[12][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[3][1][0] -= clov[1].off_diag[ 8][1]  * src[2][2][1];
      // result[3][1][1] += clov[1].off_diag[ 8][0]  * src[2][2][1];
      // result[3][2][0] -= clov[1].off_diag[12][1]  * src[2][2][1];
      // result[3][2][1] += clov[1].off_diag[12][0]  * src[2][2][1];
      s = _mm_load_ps1((float *)&(src[2][2][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);

      // result[3][1][0] += clov[1].off_diag[ 9][0]  * src[3][0][0];
      // result[3][1][1] += clov[1].off_diag[ 9][1]  * src[3][0][0];
      // result[3][2][0] += clov[1].off_diag[13][0]  * src[3][0][0];
      // result[3][2][1] += clov[1].off_diag[13][1]  * src[3][0][0];
      s = _mm_load_ps1((float *)&(src[3][0][0]));
      c = _mm_loadl_pi(c,(__m64 *)&(clov[1].off_diag[9][0]));
      c = _mm_loadh_pi(c,(__m64 *)&(clov[1].off_diag[13][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[3][1][0] -= clov[1].off_diag[ 9][1]  * src[3][0][1];
      // result[3][1][1] += clov[1].off_diag[ 9][0]  * src[3][0][1];
      // result[3][2][0] -= clov[1].off_diag[13][1]  * src[3][0][1];
      // result[3][2][1] += clov[1].off_diag[13][0]  * src[3][0][1];
      s = _mm_load_ps1((float *)&(src[3][0][1]));
      c = _mm_shuffle_ps(c,c,0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask0101.vec, tmp1);
      r = _mm_add_ps(r, tmp);

      // result[3][1][0] += clov[1].off_diag[14][0]  * src[3][2][0];
      // result[3][1][1] += clov[1].off_diag[14][0]  * src[3][2][1];
      // result[3][2][0] += clov[1].off_diag[14][0]  * src[3][1][0];
      // result[3][2][1] += clov[1].off_diag[14][0]  * src[3][1][1];
     c = _mm_load_ps1((float *)&(clov[1].off_diag[14][0]));
      s = _mm_loadl_pi(s,(__m64 *)&(src[3][2][0]));
      s = _mm_loadh_pi(s,(__m64 *)&(src[3][1][0]));
      tmp = _mm_mul_ps(c,s);
      r = _mm_add_ps(r,tmp);

      // result[3][1][0] += clov[1].off_diag[14][1]  * src[3][2][1];
      // result[3][1][1] -= clov[1].off_diag[14][1]  * src[3][2][0];
      // result[3][2][0] -= clov[1].off_diag[14][1]  * src[3][1][1];
      // result[3][2][1] += clov[1].off_diag[14][1]  * src[3][1][0];
      c = _mm_load_ps1((float *)&(clov[1].off_diag[14][1]));
      s = _mm_shuffle_ps(s,s, 0xb1);
      tmp1 = _mm_mul_ps(c,s);
      tmp = _mm_mul_ps(mask1001.vec, tmp1);
      r = _mm_add_ps(r,tmp);
     
      _mm_store_ps((float *)&(result[3][1][0]),r);	             

    }
  }
}

#endif
