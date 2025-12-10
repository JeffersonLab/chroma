#ifndef CPP_DSLASH_PARSCALAR_RECONS_32BIT_SSE2_H
#define CPP_DSLASH_PARSCALAR_RECONS_32BIT_SSE2_H

#include <xmmintrin.h>
#include <cpp_dslash_types.h>
using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;

#include <sse_sign_32bit.h>

namespace CPlusPlusWilsonDslash { 

  namespace  DslashParscalar32Bit { 


inline
void recons_4dir_plus( HalfSpinor hs0,
		       HalfSpinor hs1,
		       HalfSpinor hs2,
		       HalfSpinor hs3,
		      FourSpinor spinor)
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

  SSESign signs24 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};


  /* The partial sum is stored in right order */
  xmm0 = _mm_load_ps(&spinor[0][0][0]);
  xmm1 = _mm_load_ps(&spinor[0][2][0]);
  xmm2 = _mm_load_ps(&spinor[1][1][0]);

  xmm6 = _mm_load_ps(&spinor[2][0][0]);
  xmm7 = _mm_load_ps(&spinor[2][2][0]);
  xmm8 = _mm_load_ps(&spinor[3][1][0]);

  /* Top components. Don't need to reruct just
     accumulate */

  /* half spinor from dir 0 */
  /* index order is [color][spin][reim] */
  xmm3 = _mm_load_ps(&hs0[0][0][0]);
  xmm4 = _mm_load_ps(&hs0[1][0][0]); 
  xmm5 = _mm_load_ps(&hs0[2][0][0]);
 
  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs24.vector, xmm3);
  xmm4 = _mm_xor_ps(signs24.vector, xmm4);
  xmm5 = _mm_xor_ps(signs24.vector, xmm5);
  
  xmm6 = _mm_add_ps( xmm3, xmm6 );
  xmm7 = _mm_add_ps( xmm4, xmm7 );
  xmm8 = _mm_add_ps( xmm5, xmm8 );

  /* Half  spinor from dir 1 */
  xmm3 = _mm_load_ps(&hs1[0][0][0]);
  xmm4 = _mm_load_ps(&hs1[1][0][0]);
  xmm5 = _mm_load_ps(&hs1[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs34.vector, xmm3);
  xmm4 = _mm_xor_ps(signs34.vector, xmm4);
  xmm5 = _mm_xor_ps(signs34.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Half spinor from dir 2 */
  xmm3 = _mm_load_ps(&hs2[0][0][0]);
  xmm4 = _mm_load_ps(&hs2[1][0][0]);
  xmm5 = _mm_load_ps(&hs2[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs23.vector, xmm3);
  xmm4 = _mm_xor_ps(signs23.vector, xmm4);
  xmm5 = _mm_xor_ps(signs23.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);


  /* Half Spinor from dir 3 */
  xmm3 = _mm_load_ps(&hs3[0][0][0]);
  xmm4 = _mm_load_ps(&hs3[1][0][0]);
  xmm5 = _mm_load_ps(&hs3[2][0][0]);

  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Accumulate */
  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Deswizzle and store */
  /* Store top half of result spinor */
  /* Currently we have:                */
  /* xmm0  = [ 101 | 100 | 001 | 000 ] */
  /* xmm1  = [ 111 | 110 | 011 | 010 ] */
  /* xmm2  = [ 121 | 120 | 021 | 020 ] */

  /* Want to std::map to                    */
  /*  R0   = [ 011 | 010 | 001 | 000 ] */
  /*  R1   = [ 101 | 100 | 021 | 020 ] */
  /*  R2   = [ 121 | 120 | 111 | 110 ] */

  /* Then store as r0, r1, r2 */
  /* xmm0 has the correct low two numbers to be R0 */
  /* xmm2 has the correct low two numbers to be R1 */

  /* Let us get all of xmm2 into xmm3 -> leaves xmm3 to be R1 */
  xmm9 = _mm_xor_ps(xmm9, xmm9);
  xmm10 = _mm_xor_ps(xmm10, xmm10);
  xmm11 = _mm_xor_ps(xmm11, xmm11);

  xmm9 = _mm_add_ps(xmm9, xmm0);
  xmm9 = _mm_movelh_ps(xmm9, xmm1);

  xmm10 = _mm_add_ps(xmm10, xmm2);
  xmm10 = _mm_shuffle_ps(xmm10, xmm0, 0xe4);

  xmm11 = _mm_add_ps(xmm11, xmm2);
  xmm11 = _mm_movehl_ps(xmm11, xmm1);

  xmm0 = _mm_xor_ps(xmm0,xmm0);
  xmm1 = _mm_xor_ps(xmm1,xmm1);
  xmm2 = _mm_xor_ps(xmm2,xmm2);

  xmm0 = _mm_add_ps(xmm0,xmm6);
  xmm0 = _mm_movelh_ps(xmm0,xmm7);

  xmm1 = _mm_add_ps(xmm1, xmm8);
  xmm1 = _mm_shuffle_ps(xmm1, xmm6, 0xe4);
  
  xmm2 = _mm_add_ps(xmm2, xmm8);
  xmm2 = _mm_movehl_ps(xmm2, xmm7);

  _mm_store_ps((float *)&spinor[0][0][0], xmm9);
  _mm_store_ps((float *)&spinor[0][2][0], xmm10);
  _mm_store_ps((float *)&spinor[1][1][0], xmm11);

  _mm_store_ps((float *)&spinor[2][0][0], xmm0);
  _mm_store_ps((float *)&spinor[2][2][0], xmm1);
  _mm_store_ps((float *)&spinor[3][1][0], xmm2);

}

inline
void recons_3dir_plus( HalfSpinor hs0,
		       HalfSpinor hs1,
		       HalfSpinor hs2,
		       FourSpinor spinor)
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

  SSESign signs24 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

  xmm0 = _mm_load_ps(&spinor[0][0][0]);
  xmm1 = _mm_load_ps(&spinor[0][2][0]);
  xmm2 = _mm_load_ps(&spinor[1][1][0]);
  xmm6 = _mm_load_ps(&spinor[2][0][0]);
  xmm7 = _mm_load_ps(&spinor[2][2][0]);
  xmm8 = _mm_load_ps(&spinor[3][1][0]);

  /* Top components. Don't need to reruct just
     accumulate */

  /* half spinor from dir 0 */
  /* index order is [color][spin][reim] */
  xmm3 = _mm_load_ps(&hs0[0][0][0]);
  xmm4 = _mm_load_ps(&hs0[1][0][0]); 
  xmm5 = _mm_load_ps(&hs0[2][0][0]);
 
  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs24.vector, xmm3);
  xmm4 = _mm_xor_ps(signs24.vector, xmm4);
  xmm5 = _mm_xor_ps(signs24.vector, xmm5);
  
  xmm6 = _mm_add_ps( xmm3, xmm6 );
  xmm7 = _mm_add_ps( xmm4, xmm7 );
  xmm8 = _mm_add_ps( xmm5, xmm8 );


  /* Half  spinor from dir 1 */
  xmm3 = _mm_load_ps(&hs1[0][0][0]);
  xmm4 = _mm_load_ps(&hs1[1][0][0]);
  xmm5 = _mm_load_ps(&hs1[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs34.vector, xmm3);
  xmm4 = _mm_xor_ps(signs34.vector, xmm4);
  xmm5 = _mm_xor_ps(signs34.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Half spinor from dir 2 */
  xmm3 = _mm_load_ps(&hs2[0][0][0]);
  xmm4 = _mm_load_ps(&hs2[1][0][0]);
  xmm5 = _mm_load_ps(&hs2[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs23.vector, xmm3);
  xmm4 = _mm_xor_ps(signs23.vector, xmm4);
  xmm5 = _mm_xor_ps(signs23.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);


  /* Deswizzle and store */
  /* Store top half of result spinor */
  /* Currently we have:                */
  /* xmm0  = [ 101 | 100 | 001 | 000 ] */
  /* xmm1  = [ 111 | 110 | 011 | 010 ] */
  /* xmm2  = [ 121 | 120 | 021 | 020 ] */

  /* Want to std::map to                    */
  /*  R0   = [ 011 | 010 | 001 | 000 ] */
  /*  R1   = [ 101 | 100 | 021 | 020 ] */
  /*  R2   = [ 121 | 120 | 111 | 110 ] */

  /* Then store as r0, r1, r2 */
  /* xmm0 has the correct low two numbers to be R0 */
  /* xmm2 has the correct low two numbers to be R1 */

  /* Let us get all of xmm2 into xmm3 -> leaves xmm3 to be R1 */
  xmm9 = _mm_xor_ps(xmm9, xmm9);
  xmm10 = _mm_xor_ps(xmm10, xmm10);
  xmm11 = _mm_xor_ps(xmm11, xmm11);

  xmm9 = _mm_add_ps(xmm9, xmm0);
  xmm9 = _mm_movelh_ps(xmm9, xmm1);

  xmm10 = _mm_add_ps(xmm10, xmm2);
  xmm10 = _mm_shuffle_ps(xmm10, xmm0, 0xe4);

  xmm11 = _mm_add_ps(xmm11, xmm2);
  xmm11 = _mm_movehl_ps(xmm11, xmm1);

  xmm0 = _mm_xor_ps(xmm0,xmm0);
  xmm1 = _mm_xor_ps(xmm1,xmm1);
  xmm2 = _mm_xor_ps(xmm2,xmm2);

  xmm0 = _mm_add_ps(xmm0,xmm6);
  xmm0 = _mm_movelh_ps(xmm0,xmm7);

  xmm1 = _mm_add_ps(xmm1, xmm8);
  xmm1 = _mm_shuffle_ps(xmm1, xmm6, 0xe4);
  
  xmm2 = _mm_add_ps(xmm2, xmm8);
  xmm2 = _mm_movehl_ps(xmm2, xmm7);

  _mm_store_ps((float *)&spinor[0][0][0], xmm9);
  _mm_store_ps((float *)&spinor[0][2][0], xmm10);
  _mm_store_ps((float *)&spinor[1][1][0], xmm11);

  _mm_store_ps((float *)&spinor[2][0][0], xmm0);
  _mm_store_ps((float *)&spinor[2][2][0], xmm1);
  _mm_store_ps((float *)&spinor[3][1][0], xmm2);

  
}

inline
void recons_4dir_minus( HalfSpinor hs0,
		        HalfSpinor hs1,
		        HalfSpinor hs2,
		        HalfSpinor hs3,
		       FourSpinor spinor)
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

  SSESign signs13 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};  


  xmm0 = _mm_load_ps(&spinor[0][0][0]);
  xmm1 = _mm_load_ps(&spinor[0][2][0]);
  xmm2 = _mm_load_ps(&spinor[1][1][0]);
  xmm6 = _mm_load_ps(&spinor[2][0][0]);
  xmm7 = _mm_load_ps(&spinor[2][2][0]);
  xmm8 = _mm_load_ps(&spinor[3][1][0]);


  /* Top components. Don't need to reconstruct just
     accumulate */

  /* half spinor from dir 0 */
  /* index order is [color][spin][reim] */
  xmm3 = _mm_load_ps(&hs0[0][0][0]);
  xmm4 = _mm_load_ps(&hs0[1][0][0]); 
  xmm5 = _mm_load_ps(&hs0[2][0][0]);
 
  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs13.vector, xmm3);
  xmm4 = _mm_xor_ps(signs13.vector, xmm4);
  xmm5 = _mm_xor_ps(signs13.vector, xmm5);
  
  xmm6 = _mm_add_ps( xmm3, xmm6 );
  xmm7 = _mm_add_ps( xmm4, xmm7 );
  xmm8 = _mm_add_ps( xmm5, xmm8 );


  /* Half  spinor from dir 1 */
  xmm3 = _mm_load_ps(&hs1[0][0][0]);
  xmm4 = _mm_load_ps(&hs1[1][0][0]);
  xmm5 = _mm_load_ps(&hs1[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs12.vector, xmm3);
  xmm4 = _mm_xor_ps(signs12.vector, xmm4);
  xmm5 = _mm_xor_ps(signs12.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);


  /* Half spinor from dir 2 */
  xmm3 = _mm_load_ps(&hs2[0][0][0]);
  xmm4 = _mm_load_ps(&hs2[1][0][0]);
  xmm5 = _mm_load_ps(&hs2[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs14.vector, xmm3);
  xmm4 = _mm_xor_ps(signs14.vector, xmm4);
  xmm5 = _mm_xor_ps(signs14.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Half Spinor from dir 3 */
  xmm3 = _mm_load_ps(&hs3[0][0][0]);
  xmm4 = _mm_load_ps(&hs3[1][0][0]);
  xmm5 = _mm_load_ps(&hs3[2][0][0]);

  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Accumulate */
  xmm6 = _mm_sub_ps(xmm6, xmm3);
  xmm7 = _mm_sub_ps(xmm7, xmm4);
  xmm8 = _mm_sub_ps(xmm8, xmm5);

  /* Deswizzle and store */
  /* Store top half of result spinor */
  /* Currently we have:                */
  /* xmm0  = [ 101 | 100 | 001 | 000 ] */
  /* xmm1  = [ 111 | 110 | 011 | 010 ] */
  /* xmm2  = [ 121 | 120 | 021 | 020 ] */

  /* Want to std::map to                    */
  /*  R0   = [ 011 | 010 | 001 | 000 ] */
  /*  R1   = [ 101 | 100 | 021 | 020 ] */
  /*  R2   = [ 121 | 120 | 111 | 110 ] */

  /* Then store as r0, r1, r2 */
  /* xmm0 has the correct low two numbers to be R0 */
  /* xmm2 has the correct low two numbers to be R1 */

  /* Let us get all of xmm2 into xmm3 -> leaves xmm3 to be R1 */
  xmm9 = _mm_xor_ps(xmm9, xmm9);
  xmm10 = _mm_xor_ps(xmm10, xmm10);
  xmm11 = _mm_xor_ps(xmm11, xmm11);

  xmm9 = _mm_add_ps(xmm9, xmm0);
  xmm9 = _mm_movelh_ps(xmm9, xmm1);

  xmm10 = _mm_add_ps(xmm10, xmm2);
  xmm10 = _mm_shuffle_ps(xmm10, xmm0, 0xe4);

  xmm11 = _mm_add_ps(xmm11, xmm2);
  xmm11 = _mm_movehl_ps(xmm11, xmm1);

  xmm0 = _mm_xor_ps(xmm0,xmm0);
  xmm1 = _mm_xor_ps(xmm1,xmm1);
  xmm2 = _mm_xor_ps(xmm2,xmm2);

  xmm0 = _mm_add_ps(xmm0,xmm6);
  xmm0 = _mm_movelh_ps(xmm0,xmm7);

  xmm1 = _mm_add_ps(xmm1, xmm8);
  xmm1 = _mm_shuffle_ps(xmm1, xmm6, 0xe4);
  
  xmm2 = _mm_add_ps(xmm2, xmm8);
  xmm2 = _mm_movehl_ps(xmm2, xmm7);

  _mm_store_ps((float *)&spinor[0][0][0], xmm9);
  _mm_store_ps((float *)&spinor[0][2][0], xmm10);
  _mm_store_ps((float *)&spinor[1][1][0], xmm11);

  _mm_store_ps((float *)&spinor[2][0][0], xmm0);
  _mm_store_ps((float *)&spinor[2][2][0], xmm1);
  _mm_store_ps((float *)&spinor[3][1][0], xmm2);



}

inline
void recons_3dir_minus( HalfSpinor hs0,
		        HalfSpinor hs1,
		        HalfSpinor hs2,
			FourSpinor spinor)
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

  SSESign signs13 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};  


  xmm0 = _mm_load_ps(&spinor[0][0][0]);
  xmm1 = _mm_load_ps(&spinor[0][2][0]);
  xmm2 = _mm_load_ps(&spinor[1][1][0]);
  xmm6 = _mm_load_ps(&spinor[2][0][0]);
  xmm7 = _mm_load_ps(&spinor[2][2][0]);
  xmm8 = _mm_load_ps(&spinor[3][1][0]);

  /* Top components. Don't need to reconstruct just
     accumulate */

  /* half spinor from dir 0 */
  /* index order is [color][spin][reim] */
  xmm3 = _mm_load_ps(&hs0[0][0][0]);
  xmm4 = _mm_load_ps(&hs0[1][0][0]); 
  xmm5 = _mm_load_ps(&hs0[2][0][0]);
 
  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);


  /* recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs13.vector, xmm3);
  xmm4 = _mm_xor_ps(signs13.vector, xmm4);
  xmm5 = _mm_xor_ps(signs13.vector, xmm5);
  
  xmm6 = _mm_add_ps( xmm3, xmm6 );
  xmm7 = _mm_add_ps( xmm4, xmm7 );
  xmm8 = _mm_add_ps( xmm5, xmm8 );

  /* Half  spinor from dir 1 */
  xmm3 = _mm_load_ps(&hs1[0][0][0]);
  xmm4 = _mm_load_ps(&hs1[1][0][0]);
  xmm5 = _mm_load_ps(&hs1[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs12.vector, xmm3);
  xmm4 = _mm_xor_ps(signs12.vector, xmm4);
  xmm5 = _mm_xor_ps(signs12.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Half spinor from dir 2 */
  xmm3 = _mm_load_ps(&hs2[0][0][0]);
  xmm4 = _mm_load_ps(&hs2[1][0][0]);
  xmm5 = _mm_load_ps(&hs2[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs14.vector, xmm3);
  xmm4 = _mm_xor_ps(signs14.vector, xmm4);
  xmm5 = _mm_xor_ps(signs14.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Let us get all of xmm2 into xmm3 -> leaves xmm3 to be R1 */
  xmm9 = _mm_xor_ps(xmm9, xmm9);
  xmm10 = _mm_xor_ps(xmm10, xmm10);
  xmm11 = _mm_xor_ps(xmm11, xmm11);

  xmm9 = _mm_add_ps(xmm9, xmm0);
  xmm9 = _mm_movelh_ps(xmm9, xmm1);

  xmm10 = _mm_add_ps(xmm10, xmm2);
  xmm10 = _mm_shuffle_ps(xmm10, xmm0, 0xe4);

  xmm11 = _mm_add_ps(xmm11, xmm2);
  xmm11 = _mm_movehl_ps(xmm11, xmm1);

  xmm0 = _mm_xor_ps(xmm0,xmm0);
  xmm1 = _mm_xor_ps(xmm1,xmm1);
  xmm2 = _mm_xor_ps(xmm2,xmm2);

  xmm0 = _mm_add_ps(xmm0,xmm6);
  xmm0 = _mm_movelh_ps(xmm0,xmm7);

  xmm1 = _mm_add_ps(xmm1, xmm8);
  xmm1 = _mm_shuffle_ps(xmm1, xmm6, 0xe4);
  
  xmm2 = _mm_add_ps(xmm2, xmm8);
  xmm2 = _mm_movehl_ps(xmm2, xmm7);

  _mm_store_ps((float *)&spinor[0][0][0], xmm9);
  _mm_store_ps((float *)&spinor[0][2][0], xmm10);
  _mm_store_ps((float *)&spinor[1][1][0], xmm11);

  _mm_store_ps((float *)&spinor[2][0][0], xmm0);
  _mm_store_ps((float *)&spinor[2][2][0], xmm1);
  _mm_store_ps((float *)&spinor[3][1][0], xmm2);


}


  }
}

#endif
