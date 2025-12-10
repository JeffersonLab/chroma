#ifndef CPP_DSLASH_PARSCALAR_DECOMP_HVV_32BIT_SSE2_H
#define CPP_DSLASH_PARSCALAR_DECOMP_HVV_32BIT_SSE2_H

#include <xmmintrin.h>
#include <cpp_dslash_types.h>
using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;

#include <sse_sign_32bit.h>
namespace CPlusPlusWilsonDslash { 

  namespace  DslashParscalar32Bit { 



inline
void decomp_hvv_gamma0_plus( FourSpinor src, 
			     GaugeMatrix u,
			    HalfSpinor dst)
{
  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
  SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

  /* Now the decomposition: gamma0_plus */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(xmm3, signs13.vector);
  xmm4 = _mm_xor_ps(xmm4, signs13.vector);
  xmm5 = _mm_xor_ps(xmm5, signs13.vector);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);

  /* MAT HVV BEGIN */
  /* HALF VECTOR in xmm0,1,2 on entry */
  /* Result in      xmm3,4,5 on exit */


  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
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
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
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


  /* Store up */
  _mm_store_ps(&dst[0][0][0],xmm3);
  _mm_store_ps(&dst[1][0][0],xmm4);
  _mm_store_ps(&dst[2][0][0],xmm5);
  
}

inline
void decomp_hvv_gamma1_plus( FourSpinor src, 
			     GaugeMatrix u,
			    HalfSpinor dst)
{
  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
  SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);


  /* Now the decomposition: gamma0_plus */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(xmm3, signs12.vector);
  xmm4 = _mm_xor_ps(xmm4, signs12.vector);
  xmm5 = _mm_xor_ps(xmm5, signs12.vector);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);

  /* MAT HVV BEGIN */
  /* HALF VECTOR in xmm0,1,2 on entry */
  /* Result in      xmm3,4,5 on exit */


  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
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
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
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


  /* Store up */
  _mm_store_ps(&dst[0][0][0],xmm3);
  _mm_store_ps(&dst[1][0][0],xmm4);
  _mm_store_ps(&dst[2][0][0],xmm5);
  
}

inline
void decomp_hvv_gamma2_plus( FourSpinor src, 
			     GaugeMatrix u,
			    HalfSpinor dst)
{
  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
  SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};


  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);


  /* Now the decomposition: gamma0_plus */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(xmm3, signs14.vector);
  xmm4 = _mm_xor_ps(xmm4, signs14.vector);
  xmm5 = _mm_xor_ps(xmm5, signs14.vector);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);

  /* MAT HVV BEGIN */
  /* HALF VECTOR in xmm0,1,2 on entry */
  /* Result in      xmm3,4,5 on exit */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
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
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
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


  /* Store up */
  _mm_store_ps(&dst[0][0][0],xmm3);
  _mm_store_ps(&dst[1][0][0],xmm4);
  _mm_store_ps(&dst[2][0][0],xmm5);
  
}


inline
void decomp_hvv_gamma3_plus( FourSpinor src, 
			     GaugeMatrix u,
			    HalfSpinor dst)
{
  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
  SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);

  /* MAT HVV BEGIN */
  /* HALF VECTOR in xmm0,1,2 on entry */
  /* Result in      xmm3,4,5 on exit */
  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
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
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
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


  /* Store up */
  _mm_store_ps(&dst[0][0][0],xmm3);
  _mm_store_ps(&dst[1][0][0],xmm4);
  _mm_store_ps(&dst[2][0][0],xmm5);
  
}


inline
void decomp_hvv_gamma0_minus( FourSpinor src, 
			     GaugeMatrix u,
			    HalfSpinor dst)
{

  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
  SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;

  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

  /* Now the decomposition: gamma0_minus */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(xmm3, signs24.vector);
  xmm4 = _mm_xor_ps(xmm4, signs24.vector);
  xmm5 = _mm_xor_ps(xmm5, signs24.vector);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);

  /* MAT HVV BEGIN */
  /* HALF VECTOR in xmm0,1,2 on entry */
  /* Result in      xmm3,4,5 on exit */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
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
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
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


  /* Store up */
  _mm_store_ps(&dst[0][0][0],xmm3);
  _mm_store_ps(&dst[1][0][0],xmm4);
  _mm_store_ps(&dst[2][0][0],xmm5);
  
}

inline
void decomp_hvv_gamma1_minus( FourSpinor src, 
			     GaugeMatrix u,
			    HalfSpinor dst)
{
  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
  SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};


  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

  /* Now the decomposition: gamma0_minus */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(xmm3, signs34.vector);
  xmm4 = _mm_xor_ps(xmm4, signs34.vector);
  xmm5 = _mm_xor_ps(xmm5, signs34.vector);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);

  /* MAT HVV BEGIN */
  /* HALF VECTOR in xmm0,1,2 on entry */
  /* Result in      xmm3,4,5 on exit */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
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
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
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


  /* Store up */
  _mm_store_ps(&dst[0][0][0],xmm3);
  _mm_store_ps(&dst[1][0][0],xmm4);
  _mm_store_ps(&dst[2][0][0],xmm5);
  
}

inline
void decomp_hvv_gamma2_minus( FourSpinor src, 
			     GaugeMatrix u,
			    HalfSpinor dst)
{
  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
  SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);

  /* Now the decomposition: gamma0_minus */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(xmm3, signs23.vector);
  xmm4 = _mm_xor_ps(xmm4, signs23.vector);
  xmm5 = _mm_xor_ps(xmm5, signs23.vector);

  xmm0 = _mm_add_ps(xmm0, xmm3);
  xmm1 = _mm_add_ps(xmm1, xmm4);
  xmm2 = _mm_add_ps(xmm2, xmm5);

  /* MAT HVV BEGIN */
  /* HALF VECTOR in xmm0,1,2 on entry */
  /* Result in      xmm3,4,5 on exit */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
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
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
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


  /* Store up */
  _mm_store_ps(&dst[0][0][0],xmm3);
  _mm_store_ps(&dst[1][0][0],xmm4);
  _mm_store_ps(&dst[2][0][0],xmm5);
  
}


inline
void decomp_hvv_gamma3_minus( FourSpinor src, 
			     GaugeMatrix u,
			    HalfSpinor dst)
{
  SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};
  SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};


  /* Space for upper components */
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;

  /* Space for lower components */
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 xmm6;
  __m128 xmm7;


  __m128 t1; 
  __m128 t2; 


  xmm0 = _mm_load_ps(&src[0][0][0]);
  xmm2 = _mm_load_ps(&src[0][2][0]);
  xmm6 = _mm_load_ps(&src[1][1][0]);
  
  xmm3 = _mm_load_ps(&src[2][0][0]);
  xmm5 = _mm_load_ps(&src[2][2][0]);
  xmm7 = _mm_load_ps(&src[3][1][0]);

  xmm1 = _mm_xor_ps(xmm1,xmm1); // This should zero 
  xmm4 = _mm_xor_ps(xmm4,xmm4);

  xmm1 = _mm_movelh_ps(xmm1,xmm6);
  xmm4 = _mm_movelh_ps(xmm4,xmm7);

  xmm1 = _mm_movehl_ps(xmm1, xmm0);
  xmm4 = _mm_movehl_ps(xmm4, xmm3);


  xmm0 = _mm_shuffle_ps(xmm0, xmm2, 0xe4);
  xmm3 = _mm_shuffle_ps(xmm3, xmm5, 0xe4);

  xmm2 = _mm_shuffle_ps(xmm2, xmm6, 0xe4);
  xmm5 = _mm_shuffle_ps(xmm5, xmm7, 0xe4);


  xmm0 = _mm_sub_ps(xmm0, xmm3);
  xmm1 = _mm_sub_ps(xmm1, xmm4);
  xmm2 = _mm_sub_ps(xmm2, xmm5);

  /* MAT HVV BEGIN */
  /* HALF VECTOR in xmm0,1,2 on entry */
  /* Result in      xmm3,4,5 on exit */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[0][1][0]);
  xmm4 = _mm_load_ss(&u[1][0][0]);
  xmm7 = _mm_load_ss(&u[1][2][0]);
  xmm5 = _mm_load_ss(&u[2][0][0]);
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
  xmm6 = _mm_load_ss(&u[2][1][0]);
  xmm7 = _mm_load_ss(&u[0][2][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_xor_ps(signs24.vector, xmm0);
  xmm1 = _mm_xor_ps(signs24.vector, xmm1);
  xmm2 = _mm_xor_ps(signs24.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[1][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[0][1][1] );
  xmm7 = _mm_load_ss(&u[2][0][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[0][2][1] );
  xmm6 = _mm_load_ss( &u[2][1][1] );
  xmm7 = _mm_load_ss( &u[1][2][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
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


  /* Store up */
  _mm_store_ps(&dst[0][0][0],xmm3);
  _mm_store_ps(&dst[1][0][0],xmm4);
  _mm_store_ps(&dst[2][0][0],xmm5);
  
}


  }
}

#endif
