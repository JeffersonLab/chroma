#ifndef CPP_DSLASH_PARSCALAR_DECOMP_HVV_64BIT_SSE2_H
#define CPP_DSLASH_PARSCALAR_DECOMP_HVV_64BIT_SSE2_H

#include <xmmintrin.h>
#include <cpp_dslash_types.h>
#include <sse_sign_64bit.h>

using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;

namespace CPlusPlusWilsonDslash { 
  namespace  DslashParscalar64Bit { 

    			
inline
void decomp_hvv_gamma0_plus( FourSpinor src, 
			     GaugeMatrix u,
			     HalfSpinor dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;

  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};

  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0  */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);


  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 1 */
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);

}

inline
void decomp_hvv_gamma1_plus( FourSpinor src, 
			     GaugeMatrix u,
			     HalfSpinor dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};

  
  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);


  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);



}

inline
void decomp_hvv_gamma2_plus( FourSpinor src, 
			     GaugeMatrix u,
			     HalfSpinor dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0 ,0x0 }};

  /* Projection: Munge components 0 & 2 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);

  /* Components 1 & 3 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);

}


inline
void decomp_hvv_gamma3_plus( FourSpinor src, 
			     GaugeMatrix u,
			     HalfSpinor dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};

   /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);


}


inline
void decomp_hvv_gamma0_minus( FourSpinor src, 
			      GaugeMatrix u,
			      HalfSpinor dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};

  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0 */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 1 */
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);
 
}

inline
void decomp_hvv_gamma1_minus( FourSpinor src, 
			      GaugeMatrix u,
			      HalfSpinor dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};
  
  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0 */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 1 */
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);

}

inline
void decomp_hvv_gamma2_minus( FourSpinor src, 
			      GaugeMatrix u,
			      HalfSpinor dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0, 0x0 }};

  /* Projection: Munge components 0 & 2 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0 */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);

  /* Components 1 & 3 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);
 
}


inline
void decomp_hvv_gamma3_minus( FourSpinor src, 
			      GaugeMatrix u,
			      HalfSpinor dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  __m128d xmm6 ALIGN;
  __m128d xmm7 ALIGN;
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0,0x0 }};
  
   /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  /* Store component 0 */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);


  /* Multiply */
  /* Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0]);  /* _c11_re */
  xmm6 = _mm_load_sd(&u[0][1][0]);  /* _c21_re */
  xmm4 = _mm_load_sd(&u[1][0][0]);  /* _c12_re */
  xmm7 = _mm_load_sd(&u[1][2][0]);  /* _c32_re */
  xmm5 = _mm_load_sd(&u[2][0][0]);  /* _c13_re */

  xmm3 = _mm_unpacklo_pd(xmm3, xmm3); 
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm4 = _mm_unpacklo_pd(xmm4, xmm4);
  
  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd(xmm5, xmm5);
  
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);
  
  xmm6 = _mm_load_sd(&u[2][1][0]); /* _c23_re */
  xmm7 = _mm_load_sd(&u[0][2][0]); /* _c31_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);
  
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0]); /* _c22_re */
  xmm7 = _mm_load_sd(&u[2][2][0]); /* _c33_re */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1]); /* _c11_im */
  xmm7 = _mm_load_sd(&u[1][1][1]); /* _c22_im */

  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[2][2][1] ); /* c33im */
  xmm7 = _mm_load_sd( &u[1][0][1] ); /* c12im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm2, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd( &u[0][1][1] ); /* c21im */
  xmm7 = _mm_load_sd( &u[2][0][1] ); /* c13im */
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd( &u[0][2][1] ); /* c31im */
  xmm6 = _mm_load_sd( &u[2][1][1] ); /* c23im */
  xmm7 = _mm_load_sd( &u[1][2][1] ); /* c32im */

  xmm0 = _mm_unpacklo_pd(xmm0, xmm0);  
  xmm6 = _mm_unpacklo_pd(xmm6, xmm6);
  xmm7 = _mm_unpacklo_pd(xmm7, xmm7);

  xmm0 = _mm_mul_pd(xmm2, xmm0);
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);

  xmm3 = _mm_add_pd(xmm0, xmm3);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);


  /* Store component 1 */
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);

}

  }
}

#endif
