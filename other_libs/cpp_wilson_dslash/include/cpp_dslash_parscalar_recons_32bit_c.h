#ifndef CPP_DSLASH_PARSCALAR_RECONS_32BIT_C_H
#define CPP_DSLASH_PARSCALAR_RECONS_32BIT_C_H

#include <cpp_dslash_types.h>
using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;


namespace CPlusPlusWilsonDslash { 

  namespace  DslashParscalar32Bit { 


inline
void recons_4dir_plus( HalfSpinor hs0,
		       HalfSpinor hs1,
		       HalfSpinor hs2,
		       HalfSpinor hs3,
		      FourSpinor spinor)
{

  HalfSpinor upper_sum;
  upper_sum[0][0][0] = spinor[0][0][0];
  upper_sum[0][0][1] = spinor[0][0][1];
  upper_sum[0][1][0] = spinor[0][1][0];
  upper_sum[0][1][1] = spinor[0][1][1];
  upper_sum[1][0][0] = spinor[0][2][0];
  upper_sum[1][0][1] = spinor[0][2][1];
  upper_sum[1][1][0] = spinor[1][0][0];
  upper_sum[1][1][1] = spinor[1][0][1];
  upper_sum[2][0][0] = spinor[1][1][0];
  upper_sum[2][0][1] = spinor[1][1][1];
  upper_sum[2][1][0] = spinor[1][2][0];
  upper_sum[2][1][1] = spinor[1][2][1];

  HalfSpinor lower_sum;
  lower_sum[0][0][0] = spinor[2][0][0];
  lower_sum[0][0][1] = spinor[2][0][1];
  lower_sum[0][1][0] = spinor[2][1][0];
  lower_sum[0][1][1] = spinor[2][1][1];
  lower_sum[1][0][0] = spinor[2][2][0];
  lower_sum[1][0][1] = spinor[2][2][1];
  lower_sum[1][1][0] = spinor[3][0][0];
  lower_sum[1][1][1] = spinor[3][0][1];
  lower_sum[2][0][0] = spinor[3][1][0];
  lower_sum[2][0][1] = spinor[3][1][1];
  lower_sum[2][1][0] = spinor[3][2][0];
  lower_sum[2][1][1] = spinor[3][2][1];


  /*                              ( 1  0  0 +i)  ( a0 )    ( a0 + i a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1 +i  0)  ( a1 )  = ( a1 + i a2 )
   *                    0         ( 0 -i  1  0)  ( a2 )    ( a2 - i a1 )
   *                              (-i  0  0  1)  ( a3 )    ( a3 - i a0 )
   */

  // Top components are just an add
  upper_sum[0][0][0] += hs0[0][0][0];
  upper_sum[0][0][1] += hs0[0][0][1];
  upper_sum[0][1][0] += hs0[0][1][0];
  upper_sum[0][1][1] += hs0[0][1][1];
  upper_sum[1][0][0] += hs0[1][0][0];
  upper_sum[1][0][1] += hs0[1][0][1];
  upper_sum[1][1][0] += hs0[1][1][0];
  upper_sum[1][1][1] += hs0[1][1][1];
  upper_sum[2][0][0] += hs0[2][0][0];
  upper_sum[2][0][1] += hs0[2][0][1];
  upper_sum[2][1][0] += hs0[2][1][0];
  upper_sum[2][1][1] += hs0[2][1][1];
  
  /*
   *      ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
   *      ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
   */
  lower_sum[0][0][0] += hs0[0][1][1];
  lower_sum[0][0][1] -= hs0[0][1][0];
  lower_sum[0][1][0] += hs0[0][0][1];
  lower_sum[0][1][1] -= hs0[0][0][0];
  
  
  lower_sum[1][0][0] += hs0[1][1][1];
  lower_sum[1][0][1] -= hs0[1][1][0];
  lower_sum[1][1][0] += hs0[1][0][1];
  lower_sum[1][1][1] -= hs0[1][0][0];
  
  
  lower_sum[2][0][0] += hs0[2][1][1];
  lower_sum[2][0][1] -= hs0[2][1][0];
  lower_sum[2][1][0] += hs0[2][0][1];
  lower_sum[2][1][1] -= hs0[2][0][0];

  /*                              ( 1  0  0 -1)  ( a0 )    ( a0 - a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  1  0)  ( a1 )  = ( a1 + a2 )
   *                    1         ( 0  1  1  0)  ( a2 )    ( a2 + a1 )
   *                              (-1  0  0  1)  ( a3 )    ( a3 - a0 )
   */
  // Top components are just an add
  upper_sum[0][0][0] += hs1[0][0][0];
  upper_sum[0][0][1] += hs1[0][0][1];
  upper_sum[0][1][0] += hs1[0][1][0];
  upper_sum[0][1][1] += hs1[0][1][1];
  upper_sum[1][0][0] += hs1[1][0][0];
  upper_sum[1][0][1] += hs1[1][0][1];
  upper_sum[1][1][0] += hs1[1][1][0];
  upper_sum[1][1][1] += hs1[1][1][1];
  upper_sum[2][0][0] += hs1[2][0][0];
  upper_sum[2][0][1] += hs1[2][0][1];
  upper_sum[2][1][0] += hs1[2][1][0];
  upper_sum[2][1][1] += hs1[2][1][1];



  /*
   *      ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
   *      ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
   */
  lower_sum[0][0][0] += hs1[0][1][0];
  lower_sum[0][0][1] += hs1[0][1][1];
  lower_sum[0][1][0] -= hs1[0][0][0];
  lower_sum[0][1][1] -= hs1[0][0][1];
  
  lower_sum[1][0][0] += hs1[1][1][0];
  lower_sum[1][0][1] += hs1[1][1][1];
  lower_sum[1][1][0] -= hs1[1][0][0];
  lower_sum[1][1][1] -= hs1[1][0][1];
  
  
  lower_sum[2][0][0] += hs1[2][1][0];
  lower_sum[2][0][1] += hs1[2][1][1];
  lower_sum[2][1][0] -= hs1[2][0][0];
  lower_sum[2][1][1] -= hs1[2][0][1];
  
  /*                              ( 1  0  i  0)  ( a0 )    ( a0 + i a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0 -i)  ( a1 )  = ( a1 - i a3 )
   *                    2         (-i  0  1  0)  ( a2 )    ( a2 - i a0 )
   *                              ( 0  i  0  1)  ( a3 )    ( a3 + i a1 )
   */
  // Top components are just an add
  upper_sum[0][0][0] += hs2[0][0][0];
  upper_sum[0][0][1] += hs2[0][0][1];
  upper_sum[0][1][0] += hs2[0][1][0];
  upper_sum[0][1][1] += hs2[0][1][1];
  upper_sum[1][0][0] += hs2[1][0][0];
  upper_sum[1][0][1] += hs2[1][0][1];
  upper_sum[1][1][0] += hs2[1][1][0];
  upper_sum[1][1][1] += hs2[1][1][1];
  upper_sum[2][0][0] += hs2[2][0][0];
  upper_sum[2][0][1] += hs2[2][0][1];
  upper_sum[2][1][0] += hs2[2][1][0];
  upper_sum[2][1][1] += hs2[2][1][1];

  /*
   *  ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *  ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
   */
  
  lower_sum[0][0][0] += hs2[0][0][1];
  lower_sum[0][0][1] -= hs2[0][0][0];
  lower_sum[0][1][0] -= hs2[0][1][1];
  lower_sum[0][1][1] += hs2[0][1][0];
  
  lower_sum[1][0][0] += hs2[1][0][1];
  lower_sum[1][0][1] -= hs2[1][0][0];
  lower_sum[1][1][0] -= hs2[1][1][1];
  lower_sum[1][1][1] += hs2[1][1][0];
  
  lower_sum[2][0][0] += hs2[2][0][1];
  lower_sum[2][0][1] -= hs2[2][0][0];
  lower_sum[2][1][0] -= hs2[2][1][1];
  lower_sum[2][1][1] += hs2[2][1][0];

  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
   *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
   *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
   */
  
  // The top components are just an add, but we can store at the 
  // same time. NB we do need to appropriately deswizzle here.
  // Top components are just an add - rearranged to store in linear order

  spinor[0][0][0] = upper_sum[0][0][0] + hs3[0][0][0];
  spinor[0][0][1] = upper_sum[0][0][1] + hs3[0][0][1];
  spinor[0][1][0] = upper_sum[1][0][0] + hs3[1][0][0];
  spinor[0][1][1] = upper_sum[1][0][1] + hs3[1][0][1];
  spinor[0][2][0] = upper_sum[2][0][0] + hs3[2][0][0];
  spinor[0][2][1] = upper_sum[2][0][1] + hs3[2][0][1];

  spinor[1][0][0] = upper_sum[0][1][0] + hs3[0][1][0];
  spinor[1][0][1] = upper_sum[0][1][1] + hs3[0][1][1];
  spinor[1][1][0] = upper_sum[1][1][0] + hs3[1][1][0];
  spinor[1][1][1] = upper_sum[1][1][1] + hs3[1][1][1];
  spinor[1][2][0] = upper_sum[2][1][0] + hs3[2][1][0];
  spinor[1][2][1] = upper_sum[2][1][1] + hs3[2][1][1];

  /*
   *      ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
   *      ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
   */
  // rearranged to store in linear order
  spinor[2][0][0] = lower_sum[0][0][0]+ hs3[0][0][0];
  spinor[2][0][1] = lower_sum[0][0][1]+ hs3[0][0][1];
  spinor[2][1][0] = lower_sum[1][0][0]+ hs3[1][0][0];
  spinor[2][1][1] = lower_sum[1][0][1]+ hs3[1][0][1];
  spinor[2][2][0] = lower_sum[2][0][0]+ hs3[2][0][0];
  spinor[2][2][1] = lower_sum[2][0][1]+ hs3[2][0][1]; 

  spinor[3][0][0] = lower_sum[0][1][0]+ hs3[0][1][0];
  spinor[3][0][1] = lower_sum[0][1][1]+ hs3[0][1][1];
  spinor[3][1][0] = lower_sum[1][1][0]+ hs3[1][1][0];
  spinor[3][1][1] = lower_sum[1][1][1]+ hs3[1][1][1];
  spinor[3][2][0] = lower_sum[2][1][0]+ hs3[2][1][0];
  spinor[3][2][1] = lower_sum[2][1][1]+ hs3[2][1][1];
}


inline
void recons_3dir_plus( HalfSpinor hs0,
		       HalfSpinor hs1,
		       HalfSpinor hs2,
		       FourSpinor spinor)
{
 HalfSpinor upper_sum;
  upper_sum[0][0][0] = spinor[0][0][0];
  upper_sum[0][0][1] = spinor[0][0][1];
  upper_sum[0][1][0] = spinor[0][1][0];
  upper_sum[0][1][1] = spinor[0][1][1];
  upper_sum[1][0][0] = spinor[0][2][0];
  upper_sum[1][0][1] = spinor[0][2][1];
  upper_sum[1][1][0] = spinor[1][0][0];
  upper_sum[1][1][1] = spinor[1][0][1];
  upper_sum[2][0][0] = spinor[1][1][0];
  upper_sum[2][0][1] = spinor[1][1][1];
  upper_sum[2][1][0] = spinor[1][2][0];
  upper_sum[2][1][1] = spinor[1][2][1];

  HalfSpinor lower_sum;
  lower_sum[0][0][0] = spinor[2][0][0];
  lower_sum[0][0][1] = spinor[2][0][1];
  lower_sum[0][1][0] = spinor[2][1][0];
  lower_sum[0][1][1] = spinor[2][1][1];
  lower_sum[1][0][0] = spinor[2][2][0];
  lower_sum[1][0][1] = spinor[2][2][1];
  lower_sum[1][1][0] = spinor[3][0][0];
  lower_sum[1][1][1] = spinor[3][0][1];
  lower_sum[2][0][0] = spinor[3][1][0];
  lower_sum[2][0][1] = spinor[3][1][1];
  lower_sum[2][1][0] = spinor[3][2][0];
  lower_sum[2][1][1] = spinor[3][2][1];


  /*                              ( 1  0  0 +i)  ( a0 )    ( a0 + i a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1 +i  0)  ( a1 )  = ( a1 + i a2 )
   *                    0         ( 0 -i  1  0)  ( a2 )    ( a2 - i a1 )
   *                              (-i  0  0  1)  ( a3 )    ( a3 - i a0 )
   */

  // Top components are just an add
  upper_sum[0][0][0] += hs0[0][0][0];
  upper_sum[0][0][1] += hs0[0][0][1];
  upper_sum[0][1][0] += hs0[0][1][0];
  upper_sum[0][1][1] += hs0[0][1][1];
  upper_sum[1][0][0] += hs0[1][0][0];
  upper_sum[1][0][1] += hs0[1][0][1];
  upper_sum[1][1][0] += hs0[1][1][0];
  upper_sum[1][1][1] += hs0[1][1][1];
  upper_sum[2][0][0] += hs0[2][0][0];
  upper_sum[2][0][1] += hs0[2][0][1];
  upper_sum[2][1][0] += hs0[2][1][0];
  upper_sum[2][1][1] += hs0[2][1][1];
  
  /*
   *      ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
   *      ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
   */
  lower_sum[0][0][0] += hs0[0][1][1];
  lower_sum[0][0][1] -= hs0[0][1][0];
  lower_sum[0][1][0] += hs0[0][0][1];
  lower_sum[0][1][1] -= hs0[0][0][0];
  
  
  lower_sum[1][0][0] += hs0[1][1][1];
  lower_sum[1][0][1] -= hs0[1][1][0];
  lower_sum[1][1][0] += hs0[1][0][1];
  lower_sum[1][1][1] -= hs0[1][0][0];
  
  
  lower_sum[2][0][0] += hs0[2][1][1];
  lower_sum[2][0][1] -= hs0[2][1][0];
  lower_sum[2][1][0] += hs0[2][0][1];
  lower_sum[2][1][1] -= hs0[2][0][0];

  /*                              ( 1  0  0 -1)  ( a0 )    ( a0 - a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  1  0)  ( a1 )  = ( a1 + a2 )
   *                    1         ( 0  1  1  0)  ( a2 )    ( a2 + a1 )
   *                              (-1  0  0  1)  ( a3 )    ( a3 - a0 )
   */
  // Top components are just an add
  upper_sum[0][0][0] += hs1[0][0][0];
  upper_sum[0][0][1] += hs1[0][0][1];
  upper_sum[0][1][0] += hs1[0][1][0];
  upper_sum[0][1][1] += hs1[0][1][1];
  upper_sum[1][0][0] += hs1[1][0][0];
  upper_sum[1][0][1] += hs1[1][0][1];
  upper_sum[1][1][0] += hs1[1][1][0];
  upper_sum[1][1][1] += hs1[1][1][1];
  upper_sum[2][0][0] += hs1[2][0][0];
  upper_sum[2][0][1] += hs1[2][0][1];
  upper_sum[2][1][0] += hs1[2][1][0];
  upper_sum[2][1][1] += hs1[2][1][1];



  /*
   *      ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
   *      ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
   */
  lower_sum[0][0][0] += hs1[0][1][0];
  lower_sum[0][0][1] += hs1[0][1][1];
  lower_sum[0][1][0] -= hs1[0][0][0];
  lower_sum[0][1][1] -= hs1[0][0][1];
  
  lower_sum[1][0][0] += hs1[1][1][0];
  lower_sum[1][0][1] += hs1[1][1][1];
  lower_sum[1][1][0] -= hs1[1][0][0];
  lower_sum[1][1][1] -= hs1[1][0][1];
  
  
  lower_sum[2][0][0] += hs1[2][1][0];
  lower_sum[2][0][1] += hs1[2][1][1];
  lower_sum[2][1][0] -= hs1[2][0][0];
  lower_sum[2][1][1] -= hs1[2][0][1];
  
  /*                              ( 1  0  i  0)  ( a0 )    ( a0 + i a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0 -i)  ( a1 )  = ( a1 - i a3 )
   *                    2         (-i  0  1  0)  ( a2 )    ( a2 - i a0 )
   *                              ( 0  i  0  1)  ( a3 )    ( a3 + i a1 )
   */
  // Top components are just an add
  // rearramnged to store in linear order
  spinor[0][0][0] = upper_sum[0][0][0] + hs2[0][0][0];
  spinor[0][0][1] = upper_sum[0][0][1] + hs2[0][0][1];
  spinor[0][1][0] = upper_sum[1][0][0] + hs2[1][0][0];
  spinor[0][1][1] = upper_sum[1][0][1] + hs2[1][0][1];
  spinor[0][2][0] = upper_sum[2][0][0] + hs2[2][0][0];
  spinor[0][2][1] = upper_sum[2][0][1] + hs2[2][0][1];

  spinor[1][0][0] = upper_sum[0][1][0] + hs2[0][1][0];
  spinor[1][0][1] = upper_sum[0][1][1] + hs2[0][1][1];
  spinor[1][1][0] = upper_sum[1][1][0] + hs2[1][1][0];
  spinor[1][1][1] = upper_sum[1][1][1] + hs2[1][1][1];
  spinor[1][2][0] = upper_sum[2][1][0] + hs2[2][1][0];
  spinor[1][2][1] = upper_sum[2][1][1] + hs2[2][1][1];


  /*
   *  ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *  ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
   */

  // Rearranged to store in linear order
  spinor[2][0][0] = lower_sum[0][0][0] + hs2[0][0][1];
  spinor[2][0][1] = lower_sum[0][0][1] - hs2[0][0][0];
  spinor[2][1][0] = lower_sum[1][0][0] + hs2[1][0][1];
  spinor[2][1][1] = lower_sum[1][0][1] - hs2[1][0][0];
  spinor[2][2][0] = lower_sum[2][0][0] + hs2[2][0][1];
  spinor[2][2][1] = lower_sum[2][0][1] - hs2[2][0][0];

  spinor[3][0][0] = lower_sum[0][1][0] - hs2[0][1][1];
  spinor[3][0][1] = lower_sum[0][1][1] + hs2[0][1][0];
  spinor[3][1][0] = lower_sum[1][1][0] - hs2[1][1][1];
  spinor[3][1][1] = lower_sum[1][1][1] + hs2[1][1][0];
  spinor[3][2][0] = lower_sum[2][1][0] - hs2[2][1][1];
  spinor[3][2][1] = lower_sum[2][1][1] + hs2[2][1][0];
  
}

inline
void recons_4dir_minus( HalfSpinor hs0,
		        HalfSpinor hs1,
		        HalfSpinor hs2,
		        HalfSpinor hs3,
		       FourSpinor spinor)
{


  HalfSpinor upper_sum;
  upper_sum[0][0][0] = spinor[0][0][0];
  upper_sum[0][0][1] = spinor[0][0][1];
  upper_sum[0][1][0] = spinor[0][1][0];
  upper_sum[0][1][1] = spinor[0][1][1];
  upper_sum[1][0][0] = spinor[0][2][0];
  upper_sum[1][0][1] = spinor[0][2][1];
  upper_sum[1][1][0] = spinor[1][0][0];
  upper_sum[1][1][1] = spinor[1][0][1];
  upper_sum[2][0][0] = spinor[1][1][0];
  upper_sum[2][0][1] = spinor[1][1][1];
  upper_sum[2][1][0] = spinor[1][2][0];
  upper_sum[2][1][1] = spinor[1][2][1];

  HalfSpinor lower_sum;
  lower_sum[0][0][0] = spinor[2][0][0];
  lower_sum[0][0][1] = spinor[2][0][1];
  lower_sum[0][1][0] = spinor[2][1][0];
  lower_sum[0][1][1] = spinor[2][1][1];
  lower_sum[1][0][0] = spinor[2][2][0];
  lower_sum[1][0][1] = spinor[2][2][1];
  lower_sum[1][1][0] = spinor[3][0][0];
  lower_sum[1][1][1] = spinor[3][0][1];
  lower_sum[2][0][0] = spinor[3][1][0];
  lower_sum[2][0][1] = spinor[3][1][1];
  lower_sum[2][1][0] = spinor[3][2][0];
  lower_sum[2][1][1] = spinor[3][2][1];



  /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
   *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
   *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )
   */
  
  // Top components are just an add
  upper_sum[0][0][0] += hs0[0][0][0];
  upper_sum[0][0][1] += hs0[0][0][1];

  upper_sum[0][1][0] += hs0[0][1][0];
  upper_sum[0][1][1] += hs0[0][1][1];

  upper_sum[1][0][0] += hs0[1][0][0];
  upper_sum[1][0][1] += hs0[1][0][1];

  upper_sum[1][1][0] += hs0[1][1][0];
  upper_sum[1][1][1] += hs0[1][1][1];

  upper_sum[2][0][0] += hs0[2][0][0];
  upper_sum[2][0][1] += hs0[2][0][1];

  upper_sum[2][1][0] += hs0[2][1][0];
  upper_sum[2][1][1] += hs0[2][1][1];
 
  /*
   * ( b2r + i b2i )  =  ( {a2r - a1i} + i{a2i + a1r} )  =  ( - b1i + i b1r )
   * ( b3r + i b3i )     ( {a3r - a0i} + i{a3i + a0r} )     ( - b0i + i b0r ) 
   */
     
  // col=0
  lower_sum[0][0][0] -= hs0[0][1][1];
  lower_sum[0][0][1] += hs0[0][1][0];
  lower_sum[0][1][0] -= hs0[0][0][1];
  lower_sum[0][1][1] += hs0[0][0][0];
  
  // col=1
  lower_sum[1][0][0] -= hs0[1][1][1];
  lower_sum[1][0][1] += hs0[1][1][0];
  lower_sum[1][1][0] -= hs0[1][0][1];
  lower_sum[1][1][1] += hs0[1][0][0];
  
  // col=2
  lower_sum[2][0][0] -= hs0[2][1][1];
  lower_sum[2][0][1] += hs0[2][1][0];
  lower_sum[2][1][0] -= hs0[2][0][1];
  lower_sum[2][1][1] += hs0[2][0][0];

  /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
   *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
   *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
   */

  // Top components are just an add
  upper_sum[0][0][0] += hs1[0][0][0];
  upper_sum[0][0][1] += hs1[0][0][1];
  upper_sum[0][1][0] += hs1[0][1][0];
  upper_sum[0][1][1] += hs1[0][1][1];
  upper_sum[1][0][0] += hs1[1][0][0];
  upper_sum[1][0][1] += hs1[1][0][1];
  upper_sum[1][1][0] += hs1[1][1][0];
  upper_sum[1][1][1] += hs1[1][1][1];
  upper_sum[2][0][0] += hs1[2][0][0];
  upper_sum[2][0][1] += hs1[2][0][1];
  upper_sum[2][1][0] += hs1[2][1][0];
  upper_sum[2][1][1] += hs1[2][1][1];

  /*   
   * The bottom components of be may be reconstructed using the formula
   *  ( b2r + i b2i )  =  ( {a2r - a1r} + i{a2i - a1i} )  =  ( - b1r - i b1i )
   *  ( b3r + i b3i )     ( {a3r + a0r} + i{a3i + a0i} )     (   b0r + i b0i ) 
   */
  lower_sum[0][0][0] -= hs1[0][1][0];
  lower_sum[0][0][1] -= hs1[0][1][1];
  lower_sum[0][1][0] += hs1[0][0][0];
  lower_sum[0][1][1] += hs1[0][0][1];
  
  lower_sum[1][0][0] -= hs1[1][1][0];
  lower_sum[1][0][1] -= hs1[1][1][1];
  lower_sum[1][1][0] += hs1[1][0][0];
  lower_sum[1][1][1] += hs1[1][0][1];
  
  lower_sum[2][0][0] -= hs1[2][1][0];
  lower_sum[2][0][1] -= hs1[2][1][1];
  lower_sum[2][1][0] += hs1[2][0][0];
  lower_sum[2][1][1] += hs1[2][0][1];

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )
   */

  // Top components are just an add
  upper_sum[0][0][0] += hs2[0][0][0];
  upper_sum[0][0][1] += hs2[0][0][1];
  upper_sum[0][1][0] += hs2[0][1][0];
  upper_sum[0][1][1] += hs2[0][1][1];
  upper_sum[1][0][0] += hs2[1][0][0];
  upper_sum[1][0][1] += hs2[1][0][1];
  upper_sum[1][1][0] += hs2[1][1][0];
  upper_sum[1][1][1] += hs2[1][1][1];
  upper_sum[2][0][0] += hs2[2][0][0];
  upper_sum[2][0][1] += hs2[2][0][1];
  upper_sum[2][1][0] += hs2[2][1][0];
  upper_sum[2][1][1] += hs2[2][1][1];

  /*
   * The bottom components of be may be reconstructed using the formula
   *      ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
   *      ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
   */
  

  lower_sum[0][0][0] -= hs2[0][0][1];
  lower_sum[0][0][1] += hs2[0][0][0];
  lower_sum[0][1][0] += hs2[0][1][1];
  lower_sum[0][1][1] -= hs2[0][1][0];
  
  lower_sum[1][0][0] -= hs2[1][0][1];
  lower_sum[1][0][1] += hs2[1][0][0];
  lower_sum[1][1][0] += hs2[1][1][1];
  lower_sum[1][1][1] -= hs2[1][1][0];
  
  lower_sum[2][0][0] -= hs2[2][0][1];
  lower_sum[2][0][1] += hs2[2][0][0];
  lower_sum[2][1][0] += hs2[2][1][1];
  lower_sum[2][1][1] -= hs2[2][1][0];

  /*    
   *                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
   *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
   *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
   */

  /* Upper part is just a sum, but must deswizzle into 
     normal spinor format. I have rearranged this so that
     the spinor is stored sequentially */
  spinor[0][0][0] = upper_sum[0][0][0] + hs3[0][0][0];
  spinor[0][0][1] = upper_sum[0][0][1] + hs3[0][0][1];
  spinor[0][1][0] = upper_sum[1][0][0] + hs3[1][0][0];
  spinor[0][1][1] = upper_sum[1][0][1] + hs3[1][0][1];
  spinor[0][2][0] = upper_sum[2][0][0] + hs3[2][0][0];
  spinor[0][2][1] = upper_sum[2][0][1] + hs3[2][0][1];

  spinor[1][0][0] = upper_sum[0][1][0] + hs3[0][1][0];
  spinor[1][0][1] = upper_sum[0][1][1] + hs3[0][1][1];
  spinor[1][1][0] = upper_sum[1][1][0] + hs3[1][1][0];
  spinor[1][1][1] = upper_sum[1][1][1] + hs3[1][1][1];
  spinor[1][2][0] = upper_sum[2][1][0] + hs3[2][1][0];
  spinor[1][2][1] = upper_sum[2][1][1] + hs3[2][1][1];
   
  /*
   * ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
   * ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
   */

  /* We must deswizzle the spinor into storage order. I have rearranged
     the elements so that storage is sequential... */
  spinor[2][0][0] =  lower_sum[0][0][0] - hs3[0][0][0];
  spinor[2][0][1] =  lower_sum[0][0][1] - hs3[0][0][1];
  spinor[3][0][0] =  lower_sum[0][1][0] - hs3[0][1][0];
  spinor[3][0][1] =  lower_sum[0][1][1] - hs3[0][1][1];

  spinor[2][1][0] =  lower_sum[1][0][0] - hs3[1][0][0];
  spinor[2][1][1] =  lower_sum[1][0][1] - hs3[1][0][1];
  spinor[3][1][0] =  lower_sum[1][1][0] - hs3[1][1][0];
  spinor[3][1][1] =  lower_sum[1][1][1] - hs3[1][1][1];

  spinor[2][2][0] =  lower_sum[2][0][0] - hs3[2][0][0];
  spinor[2][2][1] =  lower_sum[2][0][1] - hs3[2][0][1];
  spinor[3][2][0] =  lower_sum[2][1][0] - hs3[2][1][0];
  spinor[3][2][1] =  lower_sum[2][1][1] - hs3[2][1][1];
 
}

inline
void recons_3dir_minus( HalfSpinor hs0,
		        HalfSpinor hs1,
		        HalfSpinor hs2,
			FourSpinor spinor)
{
  HalfSpinor upper_sum;
  upper_sum[0][0][0] = spinor[0][0][0];
  upper_sum[0][0][1] = spinor[0][0][1];
  upper_sum[0][1][0] = spinor[0][1][0];
  upper_sum[0][1][1] = spinor[0][1][1];
  upper_sum[1][0][0] = spinor[0][2][0];
  upper_sum[1][0][1] = spinor[0][2][1];
  upper_sum[1][1][0] = spinor[1][0][0];
  upper_sum[1][1][1] = spinor[1][0][1];
  upper_sum[2][0][0] = spinor[1][1][0];
  upper_sum[2][0][1] = spinor[1][1][1];
  upper_sum[2][1][0] = spinor[1][2][0];
  upper_sum[2][1][1] = spinor[1][2][1];

  HalfSpinor lower_sum;
  lower_sum[0][0][0] = spinor[2][0][0];
  lower_sum[0][0][1] = spinor[2][0][1];
  lower_sum[0][1][0] = spinor[2][1][0];
  lower_sum[0][1][1] = spinor[2][1][1];
  lower_sum[1][0][0] = spinor[2][2][0];
  lower_sum[1][0][1] = spinor[2][2][1];
  lower_sum[1][1][0] = spinor[3][0][0];
  lower_sum[1][1][1] = spinor[3][0][1];
  lower_sum[2][0][0] = spinor[3][1][0];
  lower_sum[2][0][1] = spinor[3][1][1];
  lower_sum[2][1][0] = spinor[3][2][0];
  lower_sum[2][1][1] = spinor[3][2][1];



  /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
   *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
   *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )
   */
  
  // Top components are just an add
  upper_sum[0][0][0] += hs0[0][0][0];
  upper_sum[0][0][1] += hs0[0][0][1];
  upper_sum[0][1][0] += hs0[0][1][0];
  upper_sum[0][1][1] += hs0[0][1][1];
  upper_sum[1][0][0] += hs0[1][0][0];
  upper_sum[1][0][1] += hs0[1][0][1];
  upper_sum[1][1][0] += hs0[1][1][0];
  upper_sum[1][1][1] += hs0[1][1][1];
  upper_sum[2][0][0] += hs0[2][0][0];
  upper_sum[2][0][1] += hs0[2][0][1];
  upper_sum[2][1][0] += hs0[2][1][0];
  upper_sum[2][1][1] += hs0[2][1][1];
 
  /*
   * ( b2r + i b2i )  =  ( {a2r - a1i} + i{a2i + a1r} )  =  ( - b1i + i b1r )
   * ( b3r + i b3i )     ( {a3r - a0i} + i{a3i + a0r} )     ( - b0i + i b0r ) 
   */
     
  // col=0
  lower_sum[0][0][0] -= hs0[0][1][1];
  lower_sum[0][0][1] += hs0[0][1][0];
  lower_sum[0][1][0] -= hs0[0][0][1];
  lower_sum[0][1][1] += hs0[0][0][0];
  
  // col=1
  lower_sum[1][0][0] -= hs0[1][1][1];
  lower_sum[1][0][1] += hs0[1][1][0];
  lower_sum[1][1][0] -= hs0[1][0][1];
  lower_sum[1][1][1] += hs0[1][0][0];
  
  // col=2
  lower_sum[2][0][0] -= hs0[2][1][1];
  lower_sum[2][0][1] += hs0[2][1][0];
  lower_sum[2][1][0] -= hs0[2][0][1];
  lower_sum[2][1][1] += hs0[2][0][0];

  /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
   *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
   *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
   */

  // Top components are just an add
  upper_sum[0][0][0] += hs1[0][0][0];
  upper_sum[0][0][1] += hs1[0][0][1];
  upper_sum[0][1][0] += hs1[0][1][0];
  upper_sum[0][1][1] += hs1[0][1][1];
  upper_sum[1][0][0] += hs1[1][0][0];
  upper_sum[1][0][1] += hs1[1][0][1];
  upper_sum[1][1][0] += hs1[1][1][0];
  upper_sum[1][1][1] += hs1[1][1][1];
  upper_sum[2][0][0] += hs1[2][0][0];
  upper_sum[2][0][1] += hs1[2][0][1];
  upper_sum[2][1][0] += hs1[2][1][0];
  upper_sum[2][1][1] += hs1[2][1][1];

  /*   
   * The bottom components of be may be reconstructed using the formula
   *  ( b2r + i b2i )  =  ( {a2r - a1r} + i{a2i - a1i} )  =  ( - b1r - i b1i )
   *  ( b3r + i b3i )     ( {a3r + a0r} + i{a3i + a0i} )     (   b0r + i b0i ) 
   */
  lower_sum[0][0][0] -= hs1[0][1][0];
  lower_sum[0][0][1] -= hs1[0][1][1];
  lower_sum[0][1][0] += hs1[0][0][0];
  lower_sum[0][1][1] += hs1[0][0][1];
  
  lower_sum[1][0][0] -= hs1[1][1][0];
  lower_sum[1][0][1] -= hs1[1][1][1];
  lower_sum[1][1][0] += hs1[1][0][0];
  lower_sum[1][1][1] += hs1[1][0][1];
  
  lower_sum[2][0][0] -= hs1[2][1][0];
  lower_sum[2][0][1] -= hs1[2][1][1];
  lower_sum[2][1][0] += hs1[2][0][0];
  lower_sum[2][1][1] += hs1[2][0][1];

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )
   */
  /* Upper part is just a sum, but must deswizzle into 
     normal spinor format. I have rearranged this so that
     the spinor is stored sequentially */
  spinor[0][0][0] = upper_sum[0][0][0] + hs2[0][0][0];
  spinor[0][0][1] = upper_sum[0][0][1] + hs2[0][0][1];
  spinor[0][1][0] = upper_sum[1][0][0] + hs2[1][0][0];
  spinor[0][1][1] = upper_sum[1][0][1] + hs2[1][0][1];
  spinor[0][2][0] = upper_sum[2][0][0] + hs2[2][0][0];
  spinor[0][2][1] = upper_sum[2][0][1] + hs2[2][0][1];

  spinor[1][0][0] = upper_sum[0][1][0] + hs2[0][1][0];
  spinor[1][0][1] = upper_sum[0][1][1] + hs2[0][1][1];
  spinor[1][1][0] = upper_sum[1][1][0] + hs2[1][1][0];
  spinor[1][1][1] = upper_sum[1][1][1] + hs2[1][1][1];
  spinor[1][2][0] = upper_sum[2][1][0] + hs2[2][1][0];
  spinor[1][2][1] = upper_sum[2][1][1] + hs2[2][1][1];


  /*
   * The bottom components of be may be reconstructed using the formula
   *      ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
   *      ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
   */
  

  spinor[2][0][0] = lower_sum[0][0][0] - hs2[0][0][1];
  spinor[2][0][1] = lower_sum[0][0][1] + hs2[0][0][0];

  spinor[2][1][0] = lower_sum[1][0][0] - hs2[1][0][1];
  spinor[2][1][1] = lower_sum[1][0][1] + hs2[1][0][0];

  spinor[2][2][0] = lower_sum[2][0][0] - hs2[2][0][1];
  spinor[2][2][1] = lower_sum[2][0][1] + hs2[2][0][0];

  spinor[3][0][0] = lower_sum[0][1][0] + hs2[0][1][1];
  spinor[3][0][1] = lower_sum[0][1][1] - hs2[0][1][0];
  

  spinor[3][1][0] = lower_sum[1][1][0] + hs2[1][1][1];
  spinor[3][1][1] = lower_sum[1][1][1] - hs2[1][1][0];
  

  spinor[3][2][0] = lower_sum[2][1][0] + hs2[2][1][1];
  spinor[3][2][1] = lower_sum[2][1][1] - hs2[2][1][0];

}


  }
}

#endif
