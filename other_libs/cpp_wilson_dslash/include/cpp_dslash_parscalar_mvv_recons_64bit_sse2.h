#ifndef CPP_DSLASH_PARSCALAR_MVV_RECONS_64BIT_SSE2_H
#define CPP_DSLASH_PARSCALAR_MVV_RECONS_64BIT_SSE2_H

#include <xmmintrin.h>
#include <cpp_dslash_types.h>
#include <sse_sign_64bit.h>

using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;

namespace CPlusPlusWilsonDslash { 
  namespace  DslashParscalar64Bit { 
  
inline void mvv_recons_gamma0_plus( HalfSpinor src, 
			     GaugeMatrix u,
			    FourSpinor dst)
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


  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */

  /* _sse_store_up 0 */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  /* sse_store_up 1 */
  _mm_store_pd(&dst[3][0][0], xmm3);
  _mm_store_pd(&dst[3][1][0], xmm4);
  _mm_store_pd(&dst[3][2][0], xmm5);


  /* Now deal with components 2-4 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  
  
  /* End multiply */
  
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);

  /* recons */
  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  xmm3 = _mm_xor_pd( sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd( sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd( sse_sgn.vector, xmm5);

  /* Store lower component */
  _mm_store_pd(&dst[2][0][0], xmm3);
  _mm_store_pd(&dst[2][1][0], xmm4);
  _mm_store_pd(&dst[2][2][0], xmm5);


}

inline void mvv_recons_gamma1_plus_add( HalfSpinor src, 
				 GaugeMatrix u,
				FourSpinor dst)
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


  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */



  /* Load up partial sum (in dest) */
  xmm0 = _mm_load_pd(&dst[0][0][0]);
  xmm1 = _mm_load_pd(&dst[0][1][0]);
  xmm2 = _mm_load_pd(&dst[0][2][0]);

  /* Accumulate */
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* _sse_store 0 */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Load partial sum component 3 */
  xmm0 = _mm_load_pd(&dst[3][0][0]);
  xmm1 = _mm_load_pd(&dst[3][1][0]);
  xmm2 = _mm_load_pd(&dst[3][2][0]);

  /* Accumulate */
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* _sse_store 3 */
  _mm_store_pd(&dst[3][0][0], xmm0);
  _mm_store_pd(&dst[3][1][0], xmm1);
  _mm_store_pd(&dst[3][2][0], xmm2);


  /* Now deal with components 2-4 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  
  
  /* End multiply */
  xmm0 = _mm_load_pd(&dst[1][0][0]);
  xmm1 = _mm_load_pd(&dst[1][1][0]);
  xmm2 = _mm_load_pd(&dst[1][2][0]);
  
  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

  xmm0 = _mm_load_pd(&dst[2][0][0]);
  xmm1 = _mm_load_pd(&dst[2][1][0]);
  xmm2 = _mm_load_pd(&dst[2][2][0]);
  
  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  _mm_store_pd(&dst[2][0][0], xmm0);
  _mm_store_pd(&dst[2][1][0], xmm1);
  _mm_store_pd(&dst[2][2][0], xmm2);

}

inline void mvv_recons_gamma2_plus_add( HalfSpinor src, 
				 GaugeMatrix u,
				FourSpinor dst)
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


  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */



  /* Load up partial sum (in dest) */
  xmm0 = _mm_load_pd(&dst[0][0][0]);
  xmm1 = _mm_load_pd(&dst[0][1][0]);
  xmm2 = _mm_load_pd(&dst[0][2][0]);

  /* Accumulate */
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* _sse_store 0 */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Load partial sum component 3 */
  xmm0 = _mm_load_pd(&dst[2][0][0]);
  xmm1 = _mm_load_pd(&dst[2][1][0]);
  xmm2 = _mm_load_pd(&dst[2][2][0]);

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* _sse_store 3 */
  _mm_store_pd(&dst[2][0][0], xmm0);
  _mm_store_pd(&dst[2][1][0], xmm1);
  _mm_store_pd(&dst[2][2][0], xmm2);


  /* Now deal with components 2-4 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  
  
  /* End multiply */
  xmm0 = _mm_load_pd(&dst[1][0][0]);
  xmm1 = _mm_load_pd(&dst[1][1][0]);
  xmm2 = _mm_load_pd(&dst[1][2][0]);
  
  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

  xmm0 = _mm_load_pd(&dst[3][0][0]);
  xmm1 = _mm_load_pd(&dst[3][1][0]);
  xmm2 = _mm_load_pd(&dst[3][2][0]);
  

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  _mm_store_pd(&dst[3][0][0], xmm0);
  _mm_store_pd(&dst[3][1][0], xmm1);
  _mm_store_pd(&dst[3][2][0], xmm2);

}

inline void mvv_recons_gamma2_plus_add_store( HalfSpinor src, 
				       GaugeMatrix u,
				       FourSpinor sum,
				       FourSpinor dst)
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


  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */



  /* Load up partial sum (in dest) */
  xmm0 = _mm_load_pd(&sum[0][0][0]);
  xmm1 = _mm_load_pd(&sum[0][1][0]);
  xmm2 = _mm_load_pd(&sum[0][2][0]);

  /* Accumulate */
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* _sse_store 0 */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Load partial sum component 3 */
  xmm0 = _mm_load_pd(&sum[2][0][0]);
  xmm1 = _mm_load_pd(&sum[2][1][0]);
  xmm2 = _mm_load_pd(&sum[2][2][0]);

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* _sse_store 3 */
  _mm_store_pd(&dst[2][0][0], xmm0);
  _mm_store_pd(&dst[2][1][0], xmm1);
  _mm_store_pd(&dst[2][2][0], xmm2);


  /* Now deal with components 2-4 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  
  
  /* End multiply */
  xmm0 = _mm_load_pd(&sum[1][0][0]);
  xmm1 = _mm_load_pd(&sum[1][1][0]);
  xmm2 = _mm_load_pd(&sum[1][2][0]);
  
  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

  xmm0 = _mm_load_pd(&sum[3][0][0]);
  xmm1 = _mm_load_pd(&sum[3][1][0]);
  xmm2 = _mm_load_pd(&sum[3][2][0]);
  

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  _mm_store_pd(&dst[3][0][0], xmm0);
  _mm_store_pd(&dst[3][1][0], xmm1);
  _mm_store_pd(&dst[3][2][0], xmm2);

}




inline void mvv_recons_gamma3_plus_add_store( HalfSpinor src, 
			     GaugeMatrix u,
			     FourSpinor sum,
			    FourSpinor dst)
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

  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);


  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  /* End multiply */

  xmm0 = _mm_load_pd(&sum[0][0][0]);
  xmm1 = _mm_load_pd(&sum[0][1][0]);
  xmm2 = _mm_load_pd(&sum[0][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  xmm0 = _mm_load_pd(&sum[2][0][0]);
  xmm1 = _mm_load_pd(&sum[2][1][0]);
  xmm2 = _mm_load_pd(&sum[2][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  _mm_store_pd(&dst[2][0][0], xmm0);
  _mm_store_pd(&dst[2][1][0], xmm1);
  _mm_store_pd(&dst[2][2][0], xmm2);


  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */
  xmm0 = _mm_load_pd(&sum[1][0][0]);
  xmm1 = _mm_load_pd(&sum[1][1][0]);
  xmm2 = _mm_load_pd(&sum[1][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

  xmm0 = _mm_load_pd(&sum[3][0][0]);
  xmm1 = _mm_load_pd(&sum[3][1][0]);
  xmm2 = _mm_load_pd(&sum[3][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  _mm_store_pd(&dst[3][0][0], xmm0);
  _mm_store_pd(&dst[3][1][0], xmm1);
  _mm_store_pd(&dst[3][2][0], xmm2);


}



inline void mvv_recons_gamma0_minus( HalfSpinor src, 
			     GaugeMatrix u,
			    FourSpinor dst)
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
  SSEMask2 conj = {{1,-1 }};


  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */

  /* _sse_store_up 0 */
  _mm_store_pd(&dst[0][0][0], xmm3);
  _mm_store_pd(&dst[0][1][0], xmm4);
  _mm_store_pd(&dst[0][2][0], xmm5);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);

  xmm3 = _mm_mul_pd(conj.vector, xmm3);
  xmm4 = _mm_mul_pd(conj.vector, xmm4);
  xmm5 = _mm_mul_pd(conj.vector, xmm5);

  /* sse_store_up 1 */
  _mm_store_pd(&dst[3][0][0], xmm3);
  _mm_store_pd(&dst[3][1][0], xmm4);
  _mm_store_pd(&dst[3][2][0], xmm5);


  /* Now deal with components 2-4 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  
  
  /* End multiply */
  
  _mm_store_pd(&dst[1][0][0], xmm3);
  _mm_store_pd(&dst[1][1][0], xmm4);
  _mm_store_pd(&dst[1][2][0], xmm5);

  /* recons */
  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  xmm3 = _mm_mul_pd( conj.vector, xmm3);
  xmm4 = _mm_mul_pd( conj.vector, xmm4);
  xmm5 = _mm_mul_pd( conj.vector, xmm5);

  /* Store lower component */
  _mm_store_pd(&dst[2][0][0], xmm3);
  _mm_store_pd(&dst[2][1][0], xmm4);
  _mm_store_pd(&dst[2][2][0], xmm5);


}

inline void mvv_recons_gamma1_minus_add( HalfSpinor src, 
				 GaugeMatrix u,
				FourSpinor dst)
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


  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */



  /* Load up partial sum (in dest) */
  xmm0 = _mm_load_pd(&dst[0][0][0]);
  xmm1 = _mm_load_pd(&dst[0][1][0]);
  xmm2 = _mm_load_pd(&dst[0][2][0]);

  /* Accumulate */
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* _sse_store 0 */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Load partial sum component 3 */
  xmm0 = _mm_load_pd(&dst[3][0][0]);
  xmm1 = _mm_load_pd(&dst[3][1][0]);
  xmm2 = _mm_load_pd(&dst[3][2][0]);

  /* Accumulate */
  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* _sse_store 3 */
  _mm_store_pd(&dst[3][0][0], xmm0);
  _mm_store_pd(&dst[3][1][0], xmm1);
  _mm_store_pd(&dst[3][2][0], xmm2);


  /* Now deal with components 2-4 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  
  
  /* End multiply */
  xmm0 = _mm_load_pd(&dst[1][0][0]);
  xmm1 = _mm_load_pd(&dst[1][1][0]);
  xmm2 = _mm_load_pd(&dst[1][2][0]);
  
  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

  xmm0 = _mm_load_pd(&dst[2][0][0]);
  xmm1 = _mm_load_pd(&dst[2][1][0]);
  xmm2 = _mm_load_pd(&dst[2][2][0]);
  
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&dst[2][0][0], xmm0);
  _mm_store_pd(&dst[2][1][0], xmm1);
  _mm_store_pd(&dst[2][2][0], xmm2);

}


inline void mvv_recons_gamma2_minus_add( HalfSpinor src, 
				 GaugeMatrix u,
				FourSpinor dst)
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


  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */



  /* Load up partial sum (in dest) */
  xmm0 = _mm_load_pd(&dst[0][0][0]);
  xmm1 = _mm_load_pd(&dst[0][1][0]);
  xmm2 = _mm_load_pd(&dst[0][2][0]);

  /* Accumulate */
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* _sse_store 0 */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Load partial sum component 3 */
  xmm0 = _mm_load_pd(&dst[2][0][0]);
  xmm1 = _mm_load_pd(&dst[2][1][0]);
  xmm2 = _mm_load_pd(&dst[2][2][0]);

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* _sse_store 3 */
  _mm_store_pd(&dst[2][0][0], xmm0);
  _mm_store_pd(&dst[2][1][0], xmm1);
  _mm_store_pd(&dst[2][2][0], xmm2);


  /* Now deal with components 2-4 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  
  
  /* End multiply */
  xmm0 = _mm_load_pd(&dst[1][0][0]);
  xmm1 = _mm_load_pd(&dst[1][1][0]);
  xmm2 = _mm_load_pd(&dst[1][2][0]);
  
  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

  xmm0 = _mm_load_pd(&dst[3][0][0]);
  xmm1 = _mm_load_pd(&dst[3][1][0]);
  xmm2 = _mm_load_pd(&dst[3][2][0]);
  

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&dst[3][0][0], xmm0);
  _mm_store_pd(&dst[3][1][0], xmm1);
  _mm_store_pd(&dst[3][2][0], xmm2);

}


inline void mvv_recons_gamma2_minus_add_store( HalfSpinor src, 
					GaugeMatrix u,
					FourSpinor sum,
					FourSpinor dst)
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


  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */



  /* Load up partial sum (in dest) */
  xmm0 = _mm_load_pd(&sum[0][0][0]);
  xmm1 = _mm_load_pd(&sum[0][1][0]);
  xmm2 = _mm_load_pd(&sum[0][2][0]);

  /* Accumulate */
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* _sse_store 0 */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Load partial sum component 3 */
  xmm0 = _mm_load_pd(&sum[2][0][0]);
  xmm1 = _mm_load_pd(&sum[2][1][0]);
  xmm2 = _mm_load_pd(&sum[2][2][0]);

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* _sse_store 3 */
  _mm_store_pd(&dst[2][0][0], xmm0);
  _mm_store_pd(&dst[2][1][0], xmm1);
  _mm_store_pd(&dst[2][2][0], xmm2);


  /* Now deal with components 2-4 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  
  
  /* End multiply */
  xmm0 = _mm_load_pd(&sum[1][0][0]);
  xmm1 = _mm_load_pd(&sum[1][1][0]);
  xmm2 = _mm_load_pd(&sum[1][2][0]);
  
  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

  xmm0 = _mm_load_pd(&sum[3][0][0]);
  xmm1 = _mm_load_pd(&sum[3][1][0]);
  xmm2 = _mm_load_pd(&sum[3][2][0]);
  

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&dst[3][0][0], xmm0);
  _mm_store_pd(&dst[3][1][0], xmm1);
  _mm_store_pd(&dst[3][2][0], xmm2);

}




inline void mvv_recons_gamma3_minus_add_store( HalfSpinor src, 
			     GaugeMatrix u,
			     FourSpinor sum,
			    FourSpinor dst)
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

  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);


  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
  /* End multiply */

  xmm0 = _mm_load_pd(&sum[0][0][0]);
  xmm1 = _mm_load_pd(&sum[0][1][0]);
  xmm2 = _mm_load_pd(&sum[0][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  xmm0 = _mm_load_pd(&sum[2][0][0]);
  xmm1 = _mm_load_pd(&sum[2][1][0]);
  xmm2 = _mm_load_pd(&sum[2][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&dst[2][0][0], xmm0);
  _mm_store_pd(&dst[2][1][0], xmm1);
  _mm_store_pd(&dst[2][2][0], xmm2);


  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  /* SU3 Multiply */
  xmm3 = _mm_load_sd(&u[0][0][0] );
  xmm6 = _mm_load_sd(&u[1][0][0] );
  xmm4 = _mm_load_sd(&u[0][1][0] );
  xmm7 = _mm_load_sd(&u[2][1][0] );
  xmm5 = _mm_load_sd(&u[0][2][0] );

  xmm3 = _mm_unpacklo_pd( xmm3, xmm3 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm4 = _mm_unpacklo_pd( xmm4, xmm4 );

  xmm3 = _mm_mul_pd(xmm0, xmm3);
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm5 = _mm_unpacklo_pd( xmm5, xmm5 );
  xmm4 = _mm_mul_pd(xmm0, xmm4);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_mul_pd(xmm0, xmm5);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][2][0] );
  xmm7 = _mm_load_sd(&u[2][0][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm5 = _mm_add_pd(xmm6, xmm5);
  xmm3 = _mm_add_pd(xmm7, xmm3);

  xmm6 = _mm_load_sd(&u[1][1][0] );
  xmm7 = _mm_load_sd(&u[2][2][0] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm2, xmm7);
  xmm4 = _mm_add_pd(xmm6, xmm4);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm6 = _mm_load_sd(&u[0][0][1] );
  xmm7 = _mm_load_sd(&u[1][1][1] );

  xmm0 = _mm_shuffle_pd(xmm0, xmm0, 0x1);
  xmm1 = _mm_shuffle_pd(xmm1, xmm1, 0x1);
  xmm2 = _mm_shuffle_pd(xmm2, xmm2, 0x1);
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );
  
  xmm0 = _mm_xor_pd(sse_sgn.vector, xmm0);
  xmm1 = _mm_xor_pd(sse_sgn.vector, xmm1);
  xmm2 = _mm_xor_pd(sse_sgn.vector, xmm2);

  xmm6 = _mm_mul_pd(xmm0, xmm6);
  xmm7 = _mm_mul_pd(xmm1, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm4 = _mm_add_pd(xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[2][2][1] );
  xmm7 = _mm_load_sd(&u[0][1][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd( xmm2, xmm6);
  xmm7 = _mm_mul_pd( xmm0, xmm7);
  xmm5 = _mm_add_pd( xmm6, xmm5);
  xmm4 = _mm_add_pd( xmm7, xmm4);

  xmm6 = _mm_load_sd(&u[1][0][1] );
  xmm7 = _mm_load_sd(&u[0][2][1] );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm6 = _mm_mul_pd(xmm1, xmm6);
  xmm7 = _mm_mul_pd(xmm0, xmm7);
  xmm3 = _mm_add_pd(xmm6, xmm3);
  xmm5 = _mm_add_pd(xmm7, xmm5);

  xmm0 = _mm_load_sd(&u[2][0][1] );
  xmm6 = _mm_load_sd(&u[1][2][1] );
  xmm7 = _mm_load_sd(&u[2][1][1] );

  xmm0 = _mm_unpacklo_pd( xmm0, xmm0 );
  xmm6 = _mm_unpacklo_pd( xmm6, xmm6 );
  xmm7 = _mm_unpacklo_pd( xmm7, xmm7 );

  xmm0 = _mm_mul_pd( xmm2, xmm0 );
  xmm6 = _mm_mul_pd( xmm1, xmm6 );
  xmm7 = _mm_mul_pd( xmm2, xmm7 );
  xmm3 = _mm_add_pd( xmm0, xmm3 );
  xmm5 = _mm_add_pd( xmm6, xmm5 );
  xmm4 = _mm_add_pd( xmm7, xmm4 );
    
  /* End multiply */
  xmm0 = _mm_load_pd(&sum[1][0][0]);
  xmm1 = _mm_load_pd(&sum[1][1][0]);
  xmm2 = _mm_load_pd(&sum[1][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

  xmm0 = _mm_load_pd(&sum[3][0][0]);
  xmm1 = _mm_load_pd(&sum[3][1][0]);
  xmm2 = _mm_load_pd(&sum[3][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&dst[3][0][0], xmm0);
  _mm_store_pd(&dst[3][1][0], xmm1);
  _mm_store_pd(&dst[3][2][0], xmm2);


}


  }
}

#endif
