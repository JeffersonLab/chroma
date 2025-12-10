#ifndef CPP_DSLASH_PARSCALAR_MVV_RECONS_32BIT_C_H
#define CPP_DSLASH_PARSCALAR_MVV_RECONS_32BIT_C_H

#include <cpp_dslash_types.h>

using namespace CPlusPlusWilsonDslash::Dslash32BitTypes;

namespace CPlusPlusWilsonDslash { 

  namespace  DslashParscalar32Bit { 

 
  
inline
void mvv_recons_gamma0_plus(HalfSpinor src, 
			    GaugeMatrix u,
			    HalfSpinor upper_sum, HalfSpinor lower_sum)
{
   su3_mult(upper_sum, u, src);
   /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
    *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
    *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
    *                              ( i  0  0  1)  ( a3 )  
    *
    *
    * Bottom component reconstruction is: 
    *
    *      ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
    *      ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
    */
   
   // col=0
   lower_sum[0][0][0] = -upper_sum[0][1][1];
   lower_sum[0][0][1] = upper_sum[0][1][0];
   lower_sum[0][1][0] = -upper_sum[0][0][1];
   lower_sum[0][1][1] = upper_sum[0][0][0];
   
   // col=1
   lower_sum[1][0][0] = -upper_sum[1][1][1];
   lower_sum[1][0][1] = upper_sum[1][1][0];
   lower_sum[1][1][0] = -upper_sum[1][0][1];
   lower_sum[1][1][1] = upper_sum[1][0][0];
   
   // col=2
   lower_sum[2][0][0] = -upper_sum[2][1][1];
   lower_sum[2][0][1] = upper_sum[2][1][0];
   lower_sum[2][1][0] = -upper_sum[2][0][1];
   lower_sum[2][1][1] = upper_sum[2][0][0];
}
 
inline
void mvv_recons_gamma1_plus_add(HalfSpinor src, 
				GaugeMatrix u,
				HalfSpinor upper_sum, 
				HalfSpinor lower_sum)
{
  HalfSpinor proj_matvec;

  su3_mult(proj_matvec, u, src);

  /* Reconstruction */
  // Top half
  upper_sum[0][0][0] += proj_matvec[0][0][0];
  upper_sum[0][0][1] += proj_matvec[0][0][1];
  upper_sum[0][1][0] += proj_matvec[0][1][0];
  upper_sum[0][1][1] += proj_matvec[0][1][1];
  
  upper_sum[1][0][0] += proj_matvec[1][0][0];
  upper_sum[1][0][1] += proj_matvec[1][0][1];
  upper_sum[1][1][0] += proj_matvec[1][1][0];
  upper_sum[1][1][1] += proj_matvec[1][1][1];
  
  upper_sum[2][0][0] += proj_matvec[2][0][0];
  upper_sum[2][0][1] += proj_matvec[2][0][1];
  upper_sum[2][1][0] += proj_matvec[2][1][0];
  upper_sum[2][1][1] += proj_matvec[2][1][1];

  /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
   *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
   *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
   *
   * The bottom components of be may be reconstructed using the formula
   *      ( b2r + i b2i )  =  ( {a2r - a1r} + i{a2i - a1i} )  =  ( - b1r - i b1i )
   *      ( b3r + i b3i )     ( {a3r + a0r} + i{a3i + a0i} )     (   b0r + i b0i ) 
   */
  lower_sum[0][0][0] -= proj_matvec[0][1][0];
  lower_sum[0][0][1] -= proj_matvec[0][1][1];
  lower_sum[0][1][0] += proj_matvec[0][0][0];
  lower_sum[0][1][1] += proj_matvec[0][0][1];
  
  lower_sum[1][0][0] -= proj_matvec[1][1][0];
  lower_sum[1][0][1] -= proj_matvec[1][1][1];
  lower_sum[1][1][0] += proj_matvec[1][0][0];
  lower_sum[1][1][1] += proj_matvec[1][0][1];
  
  lower_sum[2][0][0] -= proj_matvec[2][1][0];
  lower_sum[2][0][1] -= proj_matvec[2][1][1];
  lower_sum[2][1][0] += proj_matvec[2][0][0];
  lower_sum[2][1][1] += proj_matvec[2][0][1];
}

inline
void mvv_recons_gamma2_plus_add(  HalfSpinor src, 
				  GaugeMatrix u,
				HalfSpinor upper_sum, 
				HalfSpinor lower_sum)
{
  HalfSpinor proj_matvec;
  su3_mult(proj_matvec, u, src);

  
  // Top half
  upper_sum[0][0][0] += proj_matvec[0][0][0];
  upper_sum[0][0][1] += proj_matvec[0][0][1];
  upper_sum[0][1][0] += proj_matvec[0][1][0];
  upper_sum[0][1][1] += proj_matvec[0][1][1];
  
  upper_sum[1][0][0] += proj_matvec[1][0][0];
  upper_sum[1][0][1] += proj_matvec[1][0][1];
  upper_sum[1][1][0] += proj_matvec[1][1][0];
  upper_sum[1][1][1] += proj_matvec[1][1][1];
  
  upper_sum[2][0][0] += proj_matvec[2][0][0];
  upper_sum[2][0][1] += proj_matvec[2][0][1];
  upper_sum[2][1][0] += proj_matvec[2][1][0];
  upper_sum[2][1][1] += proj_matvec[2][1][1];

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )
   * The bottom components of be may be reconstructed using the formula
   *   ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
   *   ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
   */
   
  lower_sum[0][0][0] -= proj_matvec[0][0][1];
  lower_sum[0][0][1] += proj_matvec[0][0][0];
  lower_sum[0][1][0] += proj_matvec[0][1][1];
  lower_sum[0][1][1] -= proj_matvec[0][1][0];
  
  lower_sum[1][0][0] -= proj_matvec[1][0][1];
  lower_sum[1][0][1] += proj_matvec[1][0][0];
  lower_sum[1][1][0] += proj_matvec[1][1][1];
  lower_sum[1][1][1] -= proj_matvec[1][1][0];
  
  lower_sum[2][0][0] -= proj_matvec[2][0][1];
  lower_sum[2][0][1] += proj_matvec[2][0][0];
  lower_sum[2][1][0] += proj_matvec[2][1][1];
  lower_sum[2][1][1] -= proj_matvec[2][1][0];
}

inline
void mvv_recons_gamma2_plus_add_store(  HalfSpinor src, 
				        GaugeMatrix u,
				        HalfSpinor upper_sum, 
				        HalfSpinor lower_sum,
				FourSpinor dst)
{
  HalfSpinor proj_matvec;
  su3_mult(proj_matvec, u, src);

  
  // Top half -- NB dst is not in natural order... 
  // This is because I am imitating the SSE where we
  // didn't yet undeswizzle at this point. This is just 
  // a straight memcpy for now.
  dst[0][0][0] = upper_sum[0][0][0] + proj_matvec[0][0][0];
  dst[0][0][1] = upper_sum[0][0][1] + proj_matvec[0][0][1];
  dst[0][1][0] = upper_sum[0][1][0] + proj_matvec[0][1][0];
  dst[0][1][1] = upper_sum[0][1][1] + proj_matvec[0][1][1];
  dst[0][2][0] = upper_sum[1][0][0] + proj_matvec[1][0][0];
  dst[0][2][1] = upper_sum[1][0][1] + proj_matvec[1][0][1];
  dst[1][0][0] = upper_sum[1][1][0] + proj_matvec[1][1][0];
  dst[1][0][1] = upper_sum[1][1][1] + proj_matvec[1][1][1];
  dst[1][1][0] = upper_sum[2][0][0] + proj_matvec[2][0][0];
  dst[1][1][1] = upper_sum[2][0][1] + proj_matvec[2][0][1];
  dst[1][2][0] = upper_sum[2][1][0] + proj_matvec[2][1][0];
  dst[1][2][1] = upper_sum[2][1][1] + proj_matvec[2][1][1];

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )
   * The bottom components of be may be reconstructed using the formula
   *   ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
   *   ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
   */
   
  dst[2][0][0] = lower_sum[0][0][0] - proj_matvec[0][0][1];
  dst[2][0][1] = lower_sum[0][0][1] + proj_matvec[0][0][0];
  dst[2][1][0] = lower_sum[0][1][0] + proj_matvec[0][1][1];
  dst[2][1][1] = lower_sum[0][1][1] - proj_matvec[0][1][0];
  
  dst[2][2][0] = lower_sum[1][0][0] - proj_matvec[1][0][1];
  dst[2][2][1] = lower_sum[1][0][1] + proj_matvec[1][0][0];
  dst[3][0][0] = lower_sum[1][1][0] + proj_matvec[1][1][1];
  dst[3][0][1] = lower_sum[1][1][1] - proj_matvec[1][1][0];
  
  dst[3][1][0] = lower_sum[2][0][0] - proj_matvec[2][0][1];
  dst[3][1][1] = lower_sum[2][0][1] + proj_matvec[2][0][0];
  dst[3][2][0] = lower_sum[2][1][0] + proj_matvec[2][1][1];
  dst[3][2][1] = lower_sum[2][1][1] - proj_matvec[2][1][0];


}

inline
void mvv_recons_gamma3_plus_add_store(  HalfSpinor src, 
			      GaugeMatrix u,
			      HalfSpinor upper_sum, 
			      HalfSpinor lower_sum,
			    FourSpinor dst)
{

  HalfSpinor proj_matvec;
  su3_mult(proj_matvec, u, src);

  // Top half -- NB dst is not in natural order... 
  // This is because I am imitating the SSE where we
  // didn't yet undeswizzle at this point. This is just 
  // a straight memcpy for now.
  dst[0][0][0] = upper_sum[0][0][0] + proj_matvec[0][0][0];
  dst[0][0][1] = upper_sum[0][0][1] + proj_matvec[0][0][1];
  dst[0][1][0] = upper_sum[0][1][0] + proj_matvec[0][1][0];
  dst[0][1][1] = upper_sum[0][1][1] + proj_matvec[0][1][1];
  dst[0][2][0] = upper_sum[1][0][0] + proj_matvec[1][0][0];
  dst[0][2][1] = upper_sum[1][0][1] + proj_matvec[1][0][1];
  dst[1][0][0] = upper_sum[1][1][0] + proj_matvec[1][1][0];
  dst[1][0][1] = upper_sum[1][1][1] + proj_matvec[1][1][1];
  dst[1][1][0] = upper_sum[2][0][0] + proj_matvec[2][0][0];
  dst[1][1][1] = upper_sum[2][0][1] + proj_matvec[2][0][1];
  dst[1][2][0] = upper_sum[2][1][0] + proj_matvec[2][1][0];
  dst[1][2][1] = upper_sum[2][1][1] + proj_matvec[2][1][1];

  /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
   *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
   *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
   
   *      ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
   *      ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
   */
  
  dst[2][0][0] = lower_sum[0][0][0] - proj_matvec[0][0][0];
  dst[2][0][1] = lower_sum[0][0][1] - proj_matvec[0][0][1];
  dst[2][1][0] = lower_sum[0][1][0] - proj_matvec[0][1][0];
  dst[2][1][1] =  lower_sum[0][1][1] - proj_matvec[0][1][1];

  dst[2][2][0] = lower_sum[1][0][0] - proj_matvec[1][0][0];
  dst[2][2][1] = lower_sum[1][0][1] - proj_matvec[1][0][1];
  dst[3][0][0] = lower_sum[1][1][0] - proj_matvec[1][1][0];
  dst[3][0][1] = lower_sum[1][1][1] - proj_matvec[1][1][1];

  dst[3][1][0] = lower_sum[2][0][0] - proj_matvec[2][0][0];
  dst[3][1][1] = lower_sum[2][0][1] - proj_matvec[2][0][1];
  dst[3][2][0] = lower_sum[2][1][0] - proj_matvec[2][1][0];
  dst[3][2][1] = lower_sum[2][1][1] - proj_matvec[2][1][1];
}



inline
void mvv_recons_gamma0_minus(  HalfSpinor src, 
			      GaugeMatrix u,
			    HalfSpinor upper_sum, HalfSpinor lower_sum)
{
  su3_mult(upper_sum, u, src);

  /*                              ( 1  0  0 +i)  ( a0 )    ( a0 + i a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1 +i  0)  ( a1 )  = ( a1 + i a2 )
   *                    0         ( 0 -i  1  0)  ( a2 )    ( a2 - i a1 )
   *                              (-i  0  0  1)  ( a3 )    ( a3 - i a0 )
   *
   *      ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
   *      ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
   */
  lower_sum[0][0][0] = upper_sum[0][1][1];
  lower_sum[0][0][1] = -upper_sum[0][1][0];
  lower_sum[0][1][0] = upper_sum[0][0][1];
  lower_sum[0][1][1] = -upper_sum[0][0][0];
  
  
  lower_sum[1][0][0] = upper_sum[1][1][1];
  lower_sum[1][0][1] = -upper_sum[1][1][0];
  lower_sum[1][1][0] = upper_sum[1][0][1];
  lower_sum[1][1][1] = -upper_sum[1][0][0];
  
  lower_sum[2][0][0] = upper_sum[2][1][1];
  lower_sum[2][0][1] = -upper_sum[2][1][0];
  lower_sum[2][1][0] = upper_sum[2][0][1];
  lower_sum[2][1][1] = -upper_sum[2][0][0];
  
}

inline
void mvv_recons_gamma1_minus_add(  HalfSpinor src, 
				  GaugeMatrix u,
				HalfSpinor upper_sum, 
				HalfSpinor lower_sum)
{
  HalfSpinor proj_matvec;
  su3_mult(proj_matvec, u, src);

    // Top half
  upper_sum[0][0][0] += proj_matvec[0][0][0];
  upper_sum[0][0][1] += proj_matvec[0][0][1];
  upper_sum[0][1][0] += proj_matvec[0][1][0];
  upper_sum[0][1][1] += proj_matvec[0][1][1];
  
  upper_sum[1][0][0] += proj_matvec[1][0][0];
  upper_sum[1][0][1] += proj_matvec[1][0][1];
  upper_sum[1][1][0] += proj_matvec[1][1][0];
  upper_sum[1][1][1] += proj_matvec[1][1][1];
  
  upper_sum[2][0][0] += proj_matvec[2][0][0];
  upper_sum[2][0][1] += proj_matvec[2][0][1];
  upper_sum[2][1][0] += proj_matvec[2][1][0];
  upper_sum[2][1][1] += proj_matvec[2][1][1];
  

  /*                              ( 1  0  0 -1)  ( a0 )    ( a0 - a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  1  0)  ( a1 )  = ( a1 + a2 )
   *                    1         ( 0  1  1  0)  ( a2 )    ( a2 + a1 )
   *                              (-1  0  0  1)  ( a3 )    ( a3 - a0 )
   * ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
   * ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
   */

    lower_sum[0][0][0] += proj_matvec[0][1][0];
    lower_sum[0][0][1] += proj_matvec[0][1][1];
    lower_sum[0][1][0] -= proj_matvec[0][0][0];
    lower_sum[0][1][1] -= proj_matvec[0][0][1];
    
    lower_sum[1][0][0] += proj_matvec[1][1][0];
    lower_sum[1][0][1] += proj_matvec[1][1][1];
    lower_sum[1][1][0] -= proj_matvec[1][0][0];
    lower_sum[1][1][1] -= proj_matvec[1][0][1];
    

    lower_sum[2][0][0] += proj_matvec[2][1][0];
    lower_sum[2][0][1] += proj_matvec[2][1][1];
    lower_sum[2][1][0] -= proj_matvec[2][0][0];
    lower_sum[2][1][1] -= proj_matvec[2][0][1];
}

inline
void mvv_recons_gamma2_minus_add(  HalfSpinor src, 
				   GaugeMatrix u,
				   HalfSpinor upper_sum, 
				   HalfSpinor lower_sum)
{

  HalfSpinor proj_matvec;
  
  su3_mult(proj_matvec, u, src);

  // Top half
  upper_sum[0][0][0] += proj_matvec[0][0][0];
  upper_sum[0][0][1] += proj_matvec[0][0][1];
  upper_sum[0][1][0] += proj_matvec[0][1][0];
  upper_sum[0][1][1] += proj_matvec[0][1][1];
  
  upper_sum[1][0][0] += proj_matvec[1][0][0];
  upper_sum[1][0][1] += proj_matvec[1][0][1];
  upper_sum[1][1][0] += proj_matvec[1][1][0];
  upper_sum[1][1][1] += proj_matvec[1][1][1];
  
  upper_sum[2][0][0] += proj_matvec[2][0][0];
  upper_sum[2][0][1] += proj_matvec[2][0][1];
  upper_sum[2][1][0] += proj_matvec[2][1][0];
  upper_sum[2][1][1] += proj_matvec[2][1][1];
  
  
  /*                              ( 1  0  i  0)  ( a0 )    ( a0 + i a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0 -i)  ( a1 )  = ( a1 - i a3 )
   *                    2         (-i  0  1  0)  ( a2 )    ( a2 - i a0 )
   *                              ( 0  i  0  1)  ( a3 )    ( a3 + i a1 )
   *
   *      ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *      ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
   */
  
  lower_sum[0][0][0] += proj_matvec[0][0][1];
  lower_sum[0][0][1] -= proj_matvec[0][0][0];
  lower_sum[0][1][0] -= proj_matvec[0][1][1];
  lower_sum[0][1][1] += proj_matvec[0][1][0];
  
  lower_sum[1][0][0] += proj_matvec[1][0][1];
  lower_sum[1][0][1] -= proj_matvec[1][0][0];
  lower_sum[1][1][0] -= proj_matvec[1][1][1];
  lower_sum[1][1][1] += proj_matvec[1][1][0];
  
  lower_sum[2][0][0] += proj_matvec[2][0][1];
  lower_sum[2][0][1] -= proj_matvec[2][0][0];
  lower_sum[2][1][0] -= proj_matvec[2][1][1];
  lower_sum[2][1][1] += proj_matvec[2][1][0];


}

inline
void mvv_recons_gamma2_minus_add_store(  HalfSpinor src, 
				         GaugeMatrix u,
				         HalfSpinor upper_sum, 
				         HalfSpinor lower_sum,
					 FourSpinor dst)
{
  HalfSpinor proj_matvec;
  
  su3_mult(proj_matvec, u, src);

  /* Note that dst is a 4 spinor in natural  ordering
     but we are filling it unnaturally with essentially
     a contiguous memcpy. This emulates the sse
     avoidance of deswizzling */

  dst[0][0][0] = upper_sum[0][0][0] + proj_matvec[0][0][0];
  dst[0][0][1] = upper_sum[0][0][1] + proj_matvec[0][0][1];
  dst[0][1][0] = upper_sum[0][1][0] + proj_matvec[0][1][0];
  dst[0][1][1] = upper_sum[0][1][1] + proj_matvec[0][1][1];
  dst[0][2][0] = upper_sum[1][0][0] + proj_matvec[1][0][0];
  dst[0][2][1] = upper_sum[1][0][1] + proj_matvec[1][0][1];
  dst[1][0][0] = upper_sum[1][1][0] + proj_matvec[1][1][0];
  dst[1][0][1] = upper_sum[1][1][1] + proj_matvec[1][1][1];
  dst[1][1][0] = upper_sum[2][0][0] + proj_matvec[2][0][0];
  dst[1][1][1] = upper_sum[2][0][1] + proj_matvec[2][0][1];
  dst[1][2][0] = upper_sum[2][1][0] + proj_matvec[2][1][0];
  dst[1][2][1] = upper_sum[2][1][1] + proj_matvec[2][1][1];
  
  
  /*                              ( 1  0  i  0)  ( a0 )    ( a0 + i a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0 -i)  ( a1 )  = ( a1 - i a3 )
   *                    2         (-i  0  1  0)  ( a2 )    ( a2 - i a0 )
   *                              ( 0  i  0  1)  ( a3 )    ( a3 + i a1 )
   *
   *      ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *      ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
   */
  
  dst[2][0][0] = lower_sum[0][0][0] + proj_matvec[0][0][1];
  dst[2][0][1] = lower_sum[0][0][1] - proj_matvec[0][0][0];
  dst[2][1][0] = lower_sum[0][1][0] - proj_matvec[0][1][1];
  dst[2][1][1] = lower_sum[0][1][1] + proj_matvec[0][1][0];
  
  dst[2][2][0] = lower_sum[1][0][0] + proj_matvec[1][0][1];
  dst[2][2][1] = lower_sum[1][0][1] - proj_matvec[1][0][0];
  dst[3][0][0] = lower_sum[1][1][0] - proj_matvec[1][1][1];
  dst[3][0][1] = lower_sum[1][1][1] + proj_matvec[1][1][0];
  
  dst[3][1][0] = lower_sum[2][0][0] + proj_matvec[2][0][1];
  dst[3][1][1] = lower_sum[2][0][1] - proj_matvec[2][0][0];
  dst[3][2][0] = lower_sum[2][1][0] - proj_matvec[2][1][1];
  dst[3][2][1] = lower_sum[2][1][1] + proj_matvec[2][1][0];

}

inline
void mvv_recons_gamma3_minus_add_store(  HalfSpinor src, 
					 GaugeMatrix u,
					 HalfSpinor upper_sum, 
					 HalfSpinor lower_sum,
					 FourSpinor dst)
{
  HalfSpinor proj_matvec;

   su3_mult(proj_matvec, u, src);

   // Top half
   dst[0][0][0] = upper_sum[0][0][0] + proj_matvec[0][0][0];
   dst[0][0][1] = upper_sum[0][0][1] + proj_matvec[0][0][1];
   dst[0][1][0] = upper_sum[0][1][0] + proj_matvec[0][1][0];
   dst[0][1][1] = upper_sum[0][1][1] + proj_matvec[0][1][1];
   
   dst[0][2][0] = upper_sum[1][0][0] + proj_matvec[1][0][0];
   dst[0][2][1] = upper_sum[1][0][1] + proj_matvec[1][0][1];
   dst[1][0][0] = upper_sum[1][1][0] + proj_matvec[1][1][0];
   dst[1][0][1] = upper_sum[1][1][1] + proj_matvec[1][1][1];
   
   dst[1][1][0] = upper_sum[2][0][0] + proj_matvec[2][0][0];
   dst[1][1][1] = upper_sum[2][0][1] + proj_matvec[2][0][1];
   dst[1][2][0] = upper_sum[2][1][0] + proj_matvec[2][1][0];
   dst[1][2][1] = upper_sum[2][1][1] + proj_matvec[2][1][1];

    /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
     *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
     *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
     *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
     *
     * ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
     * ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
     */
    dst[2][0][0] = lower_sum[0][0][0]+ proj_matvec[0][0][0];
    dst[2][0][1] = lower_sum[0][0][1]+ proj_matvec[0][0][1];
    dst[2][1][0] = lower_sum[0][1][0]+ proj_matvec[0][1][0];
    dst[2][1][1] = lower_sum[0][1][1]+ proj_matvec[0][1][1];

    dst[2][2][0] = lower_sum[1][0][0]+ proj_matvec[1][0][0];
    dst[2][2][1] = lower_sum[1][0][1]+ proj_matvec[1][0][1];
    dst[3][0][0] = lower_sum[1][1][0]+ proj_matvec[1][1][0];
    dst[3][0][1] = lower_sum[1][1][1]+ proj_matvec[1][1][1];

    dst[3][1][0] = lower_sum[2][0][0]+ proj_matvec[2][0][0];
    dst[3][1][1] = lower_sum[2][0][1]+ proj_matvec[2][0][1];
    dst[3][2][0] = lower_sum[2][1][0]+ proj_matvec[2][1][0];
    dst[3][2][1] = lower_sum[2][1][1]+ proj_matvec[2][1][1];
} 
 


  }
}

#endif
