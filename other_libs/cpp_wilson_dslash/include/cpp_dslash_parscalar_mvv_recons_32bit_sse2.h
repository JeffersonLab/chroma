#ifndef CPP_DSLASH_PARSCALAR_MVV_RECONS_32BIT_SSE2_H
#define CPP_DSLASH_PARSCALAR_MVV_RECONS_32BIT_SSE2_H

#include <xmmintrin.h>
#include <cpp_dslash_types.h>

using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;

#include <sse_sign_32bit.h>
namespace CPlusPlusWilsonDslash { 

  namespace  DslashParscalar32Bit { 

 
  
inline
void mvv_recons_gamma0_plus(HalfSpinor src, 
			    GaugeMatrix u,
			    HalfSpinor upper_sum, HalfSpinor lower_sum)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */
  _mm_store_ps(&upper_sum[0][0][0],xmm3);
  _mm_store_ps(&upper_sum[1][0][0],xmm4);
  _mm_store_ps(&upper_sum[2][0][0],xmm5);

  /* Recons */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);
  
  xmm3 = _mm_xor_ps(signs13.vector, xmm3);
  xmm4 = _mm_xor_ps(signs13.vector, xmm4);
  xmm5 = _mm_xor_ps(signs13.vector, xmm5);
  
  /* Store up */
  _mm_store_ps(&lower_sum[0][0][0],xmm3);
  _mm_store_ps(&lower_sum[1][0][0],xmm4);
  _mm_store_ps(&lower_sum[2][0][0],xmm5);
  
}

inline
void mvv_recons_gamma1_plus_add(HalfSpinor src, 
				GaugeMatrix u,
				HalfSpinor upper_sum, 
				HalfSpinor lower_sum)
{

  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;


  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  _mm_store_ps( &upper_sum[0][0][0],xmm0 );
  _mm_store_ps( &upper_sum[1][0][0],xmm1 );
  _mm_store_ps( &upper_sum[2][0][0],xmm2 );

  /* Load lower sum project and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);
  
  xmm3 = _mm_xor_ps(signs12.vector, xmm3);
  xmm4 = _mm_xor_ps(signs12.vector, xmm4);
  xmm5 = _mm_xor_ps(signs12.vector, xmm5);

  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  _mm_store_ps( &lower_sum[0][0][0],xmm0 );
  _mm_store_ps( &lower_sum[1][0][0],xmm1 );
  _mm_store_ps( &lower_sum[2][0][0],xmm2 );
  

}

inline
void mvv_recons_gamma2_plus_add(  HalfSpinor src, 
				  GaugeMatrix u,
				HalfSpinor upper_sum, 
				HalfSpinor lower_sum)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */


  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  _mm_store_ps( &upper_sum[0][0][0],xmm0 );
  _mm_store_ps( &upper_sum[1][0][0],xmm1 );
  _mm_store_ps( &upper_sum[2][0][0],xmm2 );

  /* Load lower sum project and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);
  
  xmm3 = _mm_xor_ps(signs14.vector, xmm3);
  xmm4 = _mm_xor_ps(signs14.vector, xmm4);
  xmm5 = _mm_xor_ps(signs14.vector, xmm5);

  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  _mm_store_ps( &lower_sum[0][0][0],xmm0 );
  _mm_store_ps( &lower_sum[1][0][0],xmm1 );
  _mm_store_ps( &lower_sum[2][0][0],xmm2 );
  


}

inline
void mvv_recons_gamma2_plus_add_store(  HalfSpinor src, 
				        GaugeMatrix u,
				        HalfSpinor upper_sum, 
				        HalfSpinor lower_sum,
				FourSpinor dst)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;


  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */


  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  /* Load lower sum project and accumulate */
  xmm6 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm7 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm8 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);
  
  xmm3 = _mm_xor_ps(signs14.vector, xmm3);
  xmm4 = _mm_xor_ps(signs14.vector, xmm4);
  xmm5 = _mm_xor_ps(signs14.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);


  /* Try not deswizzling */
  _mm_store_ps(&dst[0][0][0], xmm0);
  _mm_store_ps(&dst[0][2][0], xmm1);
  _mm_store_ps(&dst[1][1][0], xmm2);
  _mm_store_ps(&dst[2][0][0], xmm6);
  _mm_store_ps(&dst[2][2][0], xmm7);
  _mm_store_ps(&dst[3][1][0], xmm8);


}

inline
void mvv_recons_gamma3_plus_add_store(  HalfSpinor src, 
			      GaugeMatrix u,
			      HalfSpinor upper_sum, 
			      HalfSpinor lower_sum,
			    FourSpinor dst)
{

  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */


  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);



  /* Load lower sum and accumulate */
  xmm6 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm7 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm8 = _mm_load_ps( &lower_sum[2][0][0] );

  /* Recons -- sse_vector sub */
  xmm6 = _mm_sub_ps( xmm6, xmm3 );
  xmm7 = _mm_sub_ps( xmm7, xmm4 );
  xmm8 = _mm_sub_ps( xmm8, xmm5 );



  /* Try not deswizzling */
  _mm_store_ps(&dst[0][0][0], xmm0);
  _mm_store_ps(&dst[0][2][0], xmm1);
  _mm_store_ps(&dst[1][1][0], xmm2);
  _mm_store_ps(&dst[2][0][0], xmm6);
  _mm_store_ps(&dst[2][2][0], xmm7);
  _mm_store_ps(&dst[3][1][0], xmm8);




}



inline
void mvv_recons_gamma0_minus(  HalfSpinor src, 
			      GaugeMatrix u,
			    HalfSpinor upper_sum, HalfSpinor lower_sum)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;


  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }}
; SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */


  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */
  _mm_store_ps(&upper_sum[0][0][0],xmm3);
  _mm_store_ps(&upper_sum[1][0][0],xmm4);
  _mm_store_ps(&upper_sum[2][0][0],xmm5);

  /* Recons */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);
  
  xmm3 = _mm_xor_ps(signs24.vector, xmm3);
  xmm4 = _mm_xor_ps(signs24.vector, xmm4);
  xmm5 = _mm_xor_ps(signs24.vector, xmm5);
  
  /* Store up */
  _mm_store_ps(&lower_sum[0][0][0],xmm3);
  _mm_store_ps(&lower_sum[1][0][0],xmm4);
  _mm_store_ps(&lower_sum[2][0][0],xmm5);
  
}

inline
void mvv_recons_gamma1_minus_add(  HalfSpinor src, 
				  GaugeMatrix u,
				HalfSpinor upper_sum, 
				HalfSpinor lower_sum)
{

  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */


  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  _mm_store_ps( &upper_sum[0][0][0],xmm0 );
  _mm_store_ps( &upper_sum[1][0][0],xmm1 );
  _mm_store_ps( &upper_sum[2][0][0],xmm2 );


  /* Load lower sum project and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);
  
  xmm3 = _mm_xor_ps(signs34.vector, xmm3);
  xmm4 = _mm_xor_ps(signs34.vector, xmm4);
  xmm5 = _mm_xor_ps(signs34.vector, xmm5);

  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  _mm_store_ps( &lower_sum[0][0][0],xmm0 );
  _mm_store_ps( &lower_sum[1][0][0],xmm1 );
  _mm_store_ps( &lower_sum[2][0][0],xmm2 );
  

}

inline
void mvv_recons_gamma2_minus_add(  HalfSpinor src, 
				  GaugeMatrix u,
				HalfSpinor upper_sum, 
				HalfSpinor lower_sum)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

   /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  _mm_store_ps( &upper_sum[0][0][0],xmm0 );
  _mm_store_ps( &upper_sum[1][0][0],xmm1 );
  _mm_store_ps( &upper_sum[2][0][0],xmm2 );

  /* Load lower sum project and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);
  
  xmm3 = _mm_xor_ps(signs23.vector, xmm3);
  xmm4 = _mm_xor_ps(signs23.vector, xmm4);
  xmm5 = _mm_xor_ps(signs23.vector, xmm5);

  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  _mm_store_ps( &lower_sum[0][0][0],xmm0 );
  _mm_store_ps( &lower_sum[1][0][0],xmm1 );
  _mm_store_ps( &lower_sum[2][0][0],xmm2 );
  


}

inline
void mvv_recons_gamma2_minus_add_store(  HalfSpinor src, 
				         GaugeMatrix u,
				         HalfSpinor upper_sum, 
				         HalfSpinor lower_sum,
				       FourSpinor dst)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

   /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );


  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);


  /* Load lower sum project and accumulate */
  xmm6 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm7 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm8 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);
  
  xmm3 = _mm_xor_ps(signs23.vector, xmm3);
  xmm4 = _mm_xor_ps(signs23.vector, xmm4);
  xmm5 = _mm_xor_ps(signs23.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);


  /* Try not deswizzling */
  _mm_store_ps(&dst[0][0][0], xmm0);
  _mm_store_ps(&dst[0][2][0], xmm1);
  _mm_store_ps(&dst[1][1][0], xmm2);
  _mm_store_ps(&dst[2][0][0], xmm6);
  _mm_store_ps(&dst[2][2][0], xmm7);
  _mm_store_ps(&dst[3][1][0], xmm8);


}

inline
void mvv_recons_gamma3_minus_add_store(  HalfSpinor src, 
			      GaugeMatrix u,
			      HalfSpinor upper_sum, 
			      HalfSpinor lower_sum,
			    FourSpinor dst)
{

  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  __m128 xmm8;
  __m128 xmm9;
  __m128 xmm10;
  __m128 xmm11;
  __m128 xmm12;

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */


  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm10 = _mm_load_ss(&u[1][1][0]);
  xmm8 = _mm_load_ss(&u[1][2][0]);
  xmm9 = _mm_load_ss(&u[2][0][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm11 = _mm_load_ss(&u[2][2][0]);

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  
  xmm8 = _mm_shuffle_ps(xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps(xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm1, xmm8);
  xmm9 = _mm_mul_ps(xmm2, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm3 = _mm_add_ps(xmm9, xmm3);

  
  xmm10 = _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm2, xmm11);
  xmm4 = _mm_add_ps(xmm10, xmm4);
  xmm5 = _mm_add_ps(xmm11, xmm5);


  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm9 = _mm_load_ss( &u[0][1][1] );
  xmm11 = _mm_load_ss(&u[0][2][1] );
  xmm10 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );

  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs13.vector, xmm0);
  xmm1 = _mm_xor_ps(signs13.vector, xmm1);
  xmm2 = _mm_xor_ps(signs13.vector, xmm2);

  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);

  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm12 = _mm_load_ss( &u[2][0][1] );

  xmm4 = _mm_add_ps(xmm7,xmm4);

  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm8 = _mm_load_ss( &u[2][2][1] );

  xmm8 = _mm_shuffle_ps( xmm8, xmm8, 0x0);
  xmm9 = _mm_shuffle_ps( xmm9, xmm9, 0x0);
  xmm8 = _mm_mul_ps(xmm2, xmm8);
  xmm9 = _mm_mul_ps(xmm0, xmm9);
  xmm5 = _mm_add_ps(xmm8, xmm5);
  xmm4 = _mm_add_ps(xmm9, xmm4);


  xmm10 =  _mm_shuffle_ps( xmm10, xmm10, 0x0);
  xmm11 = _mm_shuffle_ps( xmm11, xmm11, 0x0);
  xmm10 = _mm_mul_ps(xmm1, xmm10);
  xmm11 = _mm_mul_ps(xmm0, xmm11);
  xmm3 = _mm_add_ps(xmm10, xmm3);
  xmm5 = _mm_add_ps(xmm11, xmm5);




  xmm0 = _mm_shuffle_ps(xmm12, xmm12, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);


  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);



  /* Load lower sum and accumulate */
  xmm6 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm7 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm8 = _mm_load_ps( &lower_sum[2][0][0] );

  /* Recons -- sse_vector sub */
  xmm6 = _mm_add_ps( xmm6, xmm3 );
  xmm7 = _mm_add_ps( xmm7, xmm4 );
  xmm8 = _mm_add_ps( xmm8, xmm5 );

  /* Try not deswizzling */
  _mm_store_ps(&dst[0][0][0], xmm0);
  _mm_store_ps(&dst[0][2][0], xmm1);
  _mm_store_ps(&dst[1][1][0], xmm2);
  _mm_store_ps(&dst[2][0][0], xmm6);
  _mm_store_ps(&dst[2][2][0], xmm7);
  _mm_store_ps(&dst[3][1][0], xmm8);



}




  }
}

#endif
