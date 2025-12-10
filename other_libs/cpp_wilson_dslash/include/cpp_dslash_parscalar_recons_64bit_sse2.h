#ifndef CPP_DSLASH_PARSCALAR_RECONS_64BIT_SSE2_H
#define CPP_DSLASH_PARSCALAR_RECONS_64BIT_SSE2_H

#include <xmmintrin.h>
#include <cpp_dslash_types.h>
#include <sse_sign_64bit.h>

using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;

namespace CPlusPlusWilsonDslash { 
  namespace  DslashParscalar64Bit { 
  

inline void recons_4dir_plus( HalfSpinor hs0,
		       HalfSpinor hs1,
		       HalfSpinor hs2,
		       HalfSpinor hs3,
		      FourSpinor spinor)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;  
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;

  SSEMask sse_sgn ALIGN = {{0x0, 0x80000000, 0x0,0x0 }};

  /* Component 0 */
  xmm3 = _mm_load_pd( &hs0[0][0][0] );
  xmm4 = _mm_load_pd( &hs0[0][1][0] );
  xmm5 = _mm_load_pd( &hs0[0][2][0] );

  xmm0 = _mm_load_pd( &spinor[0][0][0] );
  xmm1 = _mm_load_pd( &spinor[0][1][0] );
  xmm2 = _mm_load_pd( &spinor[0][2][0] );

  
  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[0][0][0] );
  xmm4 = _mm_load_pd( &hs1[0][1][0] );
  xmm5 = _mm_load_pd( &hs1[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[0][0][0] );
  xmm4 = _mm_load_pd( &hs2[0][1][0] );
  xmm5 = _mm_load_pd( &hs2[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs3[0][0][0] );
  xmm4 = _mm_load_pd( &hs3[0][1][0] );
  xmm5 = _mm_load_pd( &hs3[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[0][0][0], xmm0);
  _mm_store_pd(&spinor[0][1][0], xmm1);
  _mm_store_pd(&spinor[0][2][0], xmm2);

  /* Component 1 */

  xmm3 = _mm_load_pd( &hs0[1][0][0] );
  xmm4 = _mm_load_pd( &hs0[1][1][0] );
  xmm5 = _mm_load_pd( &hs0[1][2][0] );

  xmm0 = _mm_load_pd( &spinor[1][0][0] );
  xmm1 = _mm_load_pd( &spinor[1][1][0] );
  xmm2 = _mm_load_pd( &spinor[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[1][0][0] );
  xmm4 = _mm_load_pd( &hs1[1][1][0] );
  xmm5 = _mm_load_pd( &hs1[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[1][0][0] );
  xmm4 = _mm_load_pd( &hs2[1][1][0] );
  xmm5 = _mm_load_pd( &hs2[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs3[1][0][0] );
  xmm4 = _mm_load_pd( &hs3[1][1][0] );
  xmm5 = _mm_load_pd( &hs3[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[1][0][0], xmm0);
  _mm_store_pd(&spinor[1][1][0], xmm1);
  _mm_store_pd(&spinor[1][2][0], xmm2);

  /* Component 2 */

  xmm3 = _mm_load_pd( &hs0[1][0][0] );
  xmm4 = _mm_load_pd( &hs0[1][1][0] );
  xmm5 = _mm_load_pd( &hs0[1][2][0] );

  xmm0 = _mm_load_pd( &spinor[2][0][0] );
  xmm1 = _mm_load_pd( &spinor[2][1][0] );
  xmm2 = _mm_load_pd( &spinor[2][2][0] );


  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1); 

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[1][0][0] );
  xmm4 = _mm_load_pd( &hs1[1][1][0] );
  xmm5 = _mm_load_pd( &hs1[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[0][0][0] );
  xmm4 = _mm_load_pd( &hs2[0][1][0] );
  xmm5 = _mm_load_pd( &hs2[0][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);


  xmm3 = _mm_load_pd( &hs3[0][0][0] );
  xmm4 = _mm_load_pd( &hs3[0][1][0] );
  xmm5 = _mm_load_pd( &hs3[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[2][0][0], xmm0);
  _mm_store_pd(&spinor[2][1][0], xmm1);
  _mm_store_pd(&spinor[2][2][0], xmm2);


  xmm3 = _mm_load_pd( &hs0[0][0][0] );
  xmm4 = _mm_load_pd( &hs0[0][1][0] );
  xmm5 = _mm_load_pd( &hs0[0][2][0] );

  xmm0 = _mm_load_pd( &spinor[3][0][0] );
  xmm1 = _mm_load_pd( &spinor[3][1][0] );
  xmm2 = _mm_load_pd( &spinor[3][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);


  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[0][0][0] );
  xmm4 = _mm_load_pd( &hs1[0][1][0] );
  xmm5 = _mm_load_pd( &hs1[0][2][0] );

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[1][0][0] );
  xmm4 = _mm_load_pd( &hs2[1][1][0] );
  xmm5 = _mm_load_pd( &hs2[1][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  xmm3 = _mm_load_pd( &hs3[1][0][0] );
  xmm4 = _mm_load_pd( &hs3[1][1][0] );
  xmm5 = _mm_load_pd( &hs3[1][2][0] );

  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);

  _mm_store_pd(&spinor[3][0][0], xmm0);
  _mm_store_pd(&spinor[3][1][0], xmm1);
  _mm_store_pd(&spinor[3][2][0], xmm2);
 


}

inline void recons_3dir_plus( HalfSpinor hs0,
		       HalfSpinor hs1,
		       HalfSpinor hs2,
		       FourSpinor spinor)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;  
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;

  SSEMask sse_sgn ALIGN = {{0x0, 0x80000000, 0x0,0x0 }};

  /* Component 0 */
  xmm3 = _mm_load_pd( &hs0[0][0][0] );
  xmm4 = _mm_load_pd( &hs0[0][1][0] );
  xmm5 = _mm_load_pd( &hs0[0][2][0] );

  xmm0 = _mm_load_pd( &spinor[0][0][0] );
  xmm1 = _mm_load_pd( &spinor[0][1][0] );
  xmm2 = _mm_load_pd( &spinor[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[0][0][0] );
  xmm4 = _mm_load_pd( &hs1[0][1][0] );
  xmm5 = _mm_load_pd( &hs1[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[0][0][0] );
  xmm4 = _mm_load_pd( &hs2[0][1][0] );
  xmm5 = _mm_load_pd( &hs2[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[0][0][0], xmm0);
  _mm_store_pd(&spinor[0][1][0], xmm1);
  _mm_store_pd(&spinor[0][2][0], xmm2);

  /* Component 1 */

  xmm3 = _mm_load_pd( &hs0[1][0][0] );
  xmm4 = _mm_load_pd( &hs0[1][1][0] );
  xmm5 = _mm_load_pd( &hs0[1][2][0] );

  xmm0 = _mm_load_pd( &spinor[1][0][0] );
  xmm1 = _mm_load_pd( &spinor[1][1][0] );
  xmm2 = _mm_load_pd( &spinor[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[1][0][0] );
  xmm4 = _mm_load_pd( &hs1[1][1][0] );
  xmm5 = _mm_load_pd( &hs1[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[1][0][0] );
  xmm4 = _mm_load_pd( &hs2[1][1][0] );
  xmm5 = _mm_load_pd( &hs2[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[1][0][0], xmm0);
  _mm_store_pd(&spinor[1][1][0], xmm1);
  _mm_store_pd(&spinor[1][2][0], xmm2);

  /* Component 2 */

  xmm3 = _mm_load_pd( &hs0[1][0][0] );
  xmm4 = _mm_load_pd( &hs0[1][1][0] );
  xmm5 = _mm_load_pd( &hs0[1][2][0] );

  xmm0 = _mm_load_pd( &spinor[2][0][0] );
  xmm1 = _mm_load_pd( &spinor[2][1][0] );
  xmm2 = _mm_load_pd( &spinor[2][2][0] );


  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1); 

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[1][0][0] );
  xmm4 = _mm_load_pd( &hs1[1][1][0] );
  xmm5 = _mm_load_pd( &hs1[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[0][0][0] );
  xmm4 = _mm_load_pd( &hs2[0][1][0] );
  xmm5 = _mm_load_pd( &hs2[0][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);



  _mm_store_pd(&spinor[2][0][0], xmm0);
  _mm_store_pd(&spinor[2][1][0], xmm1);
  _mm_store_pd(&spinor[2][2][0], xmm2);


  xmm3 = _mm_load_pd( &hs0[0][0][0] );
  xmm4 = _mm_load_pd( &hs0[0][1][0] );
  xmm5 = _mm_load_pd( &hs0[0][2][0] );

  xmm0 = _mm_load_pd( &spinor[3][0][0] );
  xmm1 = _mm_load_pd( &spinor[3][1][0] );
  xmm2 = _mm_load_pd( &spinor[3][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);


  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[0][0][0] );
  xmm4 = _mm_load_pd( &hs1[0][1][0] );
  xmm5 = _mm_load_pd( &hs1[0][2][0] );

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[1][0][0] );
  xmm4 = _mm_load_pd( &hs2[1][1][0] );
  xmm5 = _mm_load_pd( &hs2[1][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm3, xmm0);
  xmm1 = _mm_add_pd(xmm4, xmm1);
  xmm2 = _mm_add_pd(xmm5, xmm2);


  _mm_store_pd(&spinor[3][0][0], xmm0);
  _mm_store_pd(&spinor[3][1][0], xmm1);
  _mm_store_pd(&spinor[3][2][0], xmm2);
 


}

inline void recons_4dir_minus( HalfSpinor hs0,
		        HalfSpinor hs1,
		        HalfSpinor hs2,
		        HalfSpinor hs3,
		       FourSpinor spinor)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;  
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;


  SSEMask sse_sgn ALIGN = {{0x0, 0x80000000, 0x0,0x0 }};

  /* Component 0 */
  xmm3 = _mm_load_pd( &hs0[0][0][0] );
  xmm4 = _mm_load_pd( &hs0[0][1][0] );
  xmm5 = _mm_load_pd( &hs0[0][2][0] );

  xmm0 = _mm_load_pd( &spinor[0][0][0] );
  xmm1 = _mm_load_pd( &spinor[0][1][0] );
  xmm2 = _mm_load_pd( &spinor[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[0][0][0] );
  xmm4 = _mm_load_pd( &hs1[0][1][0] );
  xmm5 = _mm_load_pd( &hs1[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[0][0][0] );
  xmm4 = _mm_load_pd( &hs2[0][1][0] );
  xmm5 = _mm_load_pd( &hs2[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs3[0][0][0] );
  xmm4 = _mm_load_pd( &hs3[0][1][0] );
  xmm5 = _mm_load_pd( &hs3[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[0][0][0], xmm0);
  _mm_store_pd(&spinor[0][1][0], xmm1);
  _mm_store_pd(&spinor[0][2][0], xmm2);

  /* Component 1 */

  xmm3 = _mm_load_pd( &hs0[1][0][0] );
  xmm4 = _mm_load_pd( &hs0[1][1][0] );
  xmm5 = _mm_load_pd( &hs0[1][2][0] );

  xmm0 = _mm_load_pd( &spinor[1][0][0] );
  xmm1 = _mm_load_pd( &spinor[1][1][0] );
  xmm2 = _mm_load_pd( &spinor[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[1][0][0] );
  xmm4 = _mm_load_pd( &hs1[1][1][0] );
  xmm5 = _mm_load_pd( &hs1[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[1][0][0] );
  xmm4 = _mm_load_pd( &hs2[1][1][0] );
  xmm5 = _mm_load_pd( &hs2[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs3[1][0][0] );
  xmm4 = _mm_load_pd( &hs3[1][1][0] );
  xmm5 = _mm_load_pd( &hs3[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[1][0][0], xmm0);
  _mm_store_pd(&spinor[1][1][0], xmm1);
  _mm_store_pd(&spinor[1][2][0], xmm2);

  /* Component 2 */

  xmm3 = _mm_load_pd( &hs0[1][0][0] );
  xmm4 = _mm_load_pd( &hs0[1][1][0] );
  xmm5 = _mm_load_pd( &hs0[1][2][0] );

  xmm0 = _mm_load_pd( &spinor[2][0][0] );
  xmm1 = _mm_load_pd( &spinor[2][1][0] );
  xmm2 = _mm_load_pd( &spinor[2][2][0] );


  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1); 

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[1][0][0] );
  xmm4 = _mm_load_pd( &hs1[1][1][0] );
  xmm5 = _mm_load_pd( &hs1[1][2][0] );

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[0][0][0] );
  xmm4 = _mm_load_pd( &hs2[0][1][0] );
  xmm5 = _mm_load_pd( &hs2[0][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  xmm3 = _mm_load_pd( &hs3[0][0][0] );
  xmm4 = _mm_load_pd( &hs3[0][1][0] );
  xmm5 = _mm_load_pd( &hs3[0][2][0] );

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[2][0][0], xmm0);
  _mm_store_pd(&spinor[2][1][0], xmm1);
  _mm_store_pd(&spinor[2][2][0], xmm2);

  /* Component 3 */
  xmm3 = _mm_load_pd( &hs0[0][0][0] );
  xmm4 = _mm_load_pd( &hs0[0][1][0] );
  xmm5 = _mm_load_pd( &hs0[0][2][0] );

  xmm0 = _mm_load_pd( &spinor[3][0][0] );
  xmm1 = _mm_load_pd( &spinor[3][1][0] );
  xmm2 = _mm_load_pd( &spinor[3][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);


  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[0][0][0] );
  xmm4 = _mm_load_pd( &hs1[0][1][0] );
  xmm5 = _mm_load_pd( &hs1[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[1][0][0] );
  xmm4 = _mm_load_pd( &hs2[1][1][0] );
  xmm5 = _mm_load_pd( &hs2[1][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs3[1][0][0] );
  xmm4 = _mm_load_pd( &hs3[1][1][0] );
  xmm5 = _mm_load_pd( &hs3[1][2][0] );

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  _mm_store_pd(&spinor[3][0][0], xmm0);
  _mm_store_pd(&spinor[3][1][0], xmm1);
  _mm_store_pd(&spinor[3][2][0], xmm2);
 


}

inline void recons_3dir_minus( HalfSpinor hs0,
		        HalfSpinor hs1,
		        HalfSpinor hs2,
			FourSpinor spinor)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;  
  __m128d xmm2 ALIGN;
  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;


  SSEMask sse_sgn ALIGN = {{0x0, 0x80000000, 0x0,0x0 }};

  /* Component 0 */
  xmm3 = _mm_load_pd( &hs0[0][0][0] );
  xmm4 = _mm_load_pd( &hs0[0][1][0] );
  xmm5 = _mm_load_pd( &hs0[0][2][0] );

  xmm0 = _mm_load_pd( &spinor[0][0][0] );
  xmm1 = _mm_load_pd( &spinor[0][1][0] );
  xmm2 = _mm_load_pd( &spinor[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[0][0][0] );
  xmm4 = _mm_load_pd( &hs1[0][1][0] );
  xmm5 = _mm_load_pd( &hs1[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[0][0][0] );
  xmm4 = _mm_load_pd( &hs2[0][1][0] );
  xmm5 = _mm_load_pd( &hs2[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[0][0][0], xmm0);
  _mm_store_pd(&spinor[0][1][0], xmm1);
  _mm_store_pd(&spinor[0][2][0], xmm2);

  /* Component 1 */

  xmm3 = _mm_load_pd( &hs0[1][0][0] );
  xmm4 = _mm_load_pd( &hs0[1][1][0] );
  xmm5 = _mm_load_pd( &hs0[1][2][0] );

  xmm0 = _mm_load_pd( &spinor[1][0][0] );
  xmm1 = _mm_load_pd( &spinor[1][1][0] );
  xmm2 = _mm_load_pd( &spinor[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[1][0][0] );
  xmm4 = _mm_load_pd( &hs1[1][1][0] );
  xmm5 = _mm_load_pd( &hs1[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[1][0][0] );
  xmm4 = _mm_load_pd( &hs2[1][1][0] );
  xmm5 = _mm_load_pd( &hs2[1][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[1][0][0], xmm0);
  _mm_store_pd(&spinor[1][1][0], xmm1);
  _mm_store_pd(&spinor[1][2][0], xmm2);

  /* Component 2 */

  xmm3 = _mm_load_pd( &hs0[1][0][0] );
  xmm4 = _mm_load_pd( &hs0[1][1][0] );
  xmm5 = _mm_load_pd( &hs0[1][2][0] );

  xmm0 = _mm_load_pd( &spinor[2][0][0] );
  xmm1 = _mm_load_pd( &spinor[2][1][0] );
  xmm2 = _mm_load_pd( &spinor[2][2][0] );


  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1); 

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[1][0][0] );
  xmm4 = _mm_load_pd( &hs1[1][1][0] );
  xmm5 = _mm_load_pd( &hs1[1][2][0] );

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[0][0][0] );
  xmm4 = _mm_load_pd( &hs2[0][1][0] );
  xmm5 = _mm_load_pd( &hs2[0][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[2][0][0], xmm0);
  _mm_store_pd(&spinor[2][1][0], xmm1);
  _mm_store_pd(&spinor[2][2][0], xmm2);

  /* Component 3 */
  xmm3 = _mm_load_pd( &hs0[0][0][0] );
  xmm4 = _mm_load_pd( &hs0[0][1][0] );
  xmm5 = _mm_load_pd( &hs0[0][2][0] );

  xmm0 = _mm_load_pd( &spinor[3][0][0] );
  xmm1 = _mm_load_pd( &spinor[3][1][0] );
  xmm2 = _mm_load_pd( &spinor[3][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);


  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);


  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs1[0][0][0] );
  xmm4 = _mm_load_pd( &hs1[0][1][0] );
  xmm5 = _mm_load_pd( &hs1[0][2][0] );

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  xmm3 = _mm_load_pd( &hs2[1][0][0] );
  xmm4 = _mm_load_pd( &hs2[1][1][0] );
  xmm5 = _mm_load_pd( &hs2[1][2][0] );

  xmm3 = _mm_shuffle_pd(xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd(xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd(xmm5, xmm5, 0x1);

  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);


  _mm_store_pd(&spinor[3][0][0], xmm0);
  _mm_store_pd(&spinor[3][1][0], xmm1);
  _mm_store_pd(&spinor[3][2][0], xmm2);
 


}



  }
}

#endif
