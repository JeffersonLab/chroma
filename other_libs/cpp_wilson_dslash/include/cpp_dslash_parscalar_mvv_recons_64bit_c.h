#ifndef CPP_DSLASH_PARSCALAR_MVV_RECONS_64BIT_H
#define CPP_DSLASH_PARSCALAR_MVV_RECONS_64BIT_H

#include <cpp_dslash_types.h>

using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;

namespace CPlusPlusWilsonDslash { 
  namespace  DslashParscalar64Bit { 
  
inline void mvv_recons_gamma0_plus( HalfSpinor src, 
			     GaugeMatrix u,
			     FourSpinor dst)
{
     su3_mult(dst, u, src);
    /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
     *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
     *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )
     */
     /*
      *      ( b2r + i b2i )  =  ( {a2r - a1i} + i{a2i + a1r} )  =  ( - b1i + i b1r )
      *      ( b3r + i b3i )     ( {a3r - a0i} + i{a3i + a0r} )     ( - b0i + i b0r ) 
      */

     dst[2][0][0] = -dst[1][0][1];
     dst[2][0][1] = dst[1][0][0];
     dst[3][0][0] = -dst[0][0][1];
     dst[3][0][1] = dst[0][0][0];

     dst[2][1][0] = -dst[1][1][1];
     dst[2][1][1] = dst[1][1][0];
     dst[3][1][0] = -dst[0][1][1];
     dst[3][1][1] = dst[0][1][0];

     dst[2][2][0] = -dst[1][2][1];
     dst[2][2][1] = dst[1][2][0];
     dst[3][2][0] = -dst[0][2][1];
     dst[3][2][1] = dst[0][2][0];


}
inline void mvv_recons_gamma1_plus_add( HalfSpinor src, 
				 GaugeMatrix u,
				FourSpinor dst)
{

  HalfSpinor proj_matvec;

  su3_mult(proj_matvec, u,src);


  dst[0][0][0] += proj_matvec[0][0][0];
  dst[0][0][1] += proj_matvec[0][0][1];
  dst[0][1][0] += proj_matvec[0][1][0];
  dst[0][1][1] += proj_matvec[0][1][1];  
  dst[0][2][0] += proj_matvec[0][2][0];
  dst[0][2][1] += proj_matvec[0][2][1];
  
  dst[1][0][0] += proj_matvec[1][0][0];
  dst[1][0][1] += proj_matvec[1][0][1];
  dst[1][1][0] += proj_matvec[1][1][0];
  dst[1][1][1] += proj_matvec[1][1][1];
  dst[1][2][0] += proj_matvec[1][2][0];
  dst[1][2][1] += proj_matvec[1][2][1];
    
  /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
   *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
   *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
   
   * Therefore the top components are
   
   *      ( b0r + i b0i )  =  ( {a0r + a3r} + i{a0i + a3i} )
   *      ( b1r + i b1i )     ( {a1r - a2r} + i{a1i - a2i} )
   */
  /*
   * The bottom components of be may be reconstructed using the formula
   *      ( b2r + i b2i )  =  ( {a2r - a1r} + i{a2i - a1i} )  =  ( - b1r - i b1i )
   *      ( b3r + i b3i )     ( {a3r + a0r} + i{a3i + a0i} )     (   b0r + i b0i ) 
   */
  
  dst[2][0][0] -= proj_matvec[1][0][0];
  dst[2][0][1] -= proj_matvec[1][0][1];
  dst[3][0][0] += proj_matvec[0][0][0];
  dst[3][0][1] += proj_matvec[0][0][1];
  
  dst[2][1][0] -= proj_matvec[1][1][0];
  dst[2][1][1] -= proj_matvec[1][1][1];
  dst[3][1][0] += proj_matvec[0][1][0];
  dst[3][1][1] += proj_matvec[0][1][1];
  
  dst[2][2][0] -= proj_matvec[1][2][0];
  dst[2][2][1] -= proj_matvec[1][2][1];
  dst[3][2][0] += proj_matvec[0][2][0];
  dst[3][2][1] += proj_matvec[0][2][1];

}

inline void mvv_recons_gamma2_plus_add( HalfSpinor src, 
				 GaugeMatrix u,
				FourSpinor dst)
{

  HalfSpinor proj_matvec;
    /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
     *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
     *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
     *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
     
     */

    su3_mult(proj_matvec, u, src);

    dst[0][0][0] += proj_matvec[0][0][0];
    dst[0][0][1] += proj_matvec[0][0][1];
    dst[0][1][0] += proj_matvec[0][1][0];
    dst[0][1][1] += proj_matvec[0][1][1];  
    dst[0][2][0] += proj_matvec[0][2][0];
    dst[0][2][1] += proj_matvec[0][2][1];

    dst[1][0][0] += proj_matvec[1][0][0];
    dst[1][0][1] += proj_matvec[1][0][1];
    dst[1][1][0] += proj_matvec[1][1][0];
    dst[1][1][1] += proj_matvec[1][1][1];
    dst[1][2][0] += proj_matvec[1][2][0];
    dst[1][2][1] += proj_matvec[1][2][1];

    /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
     *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
     *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
     *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
     
     */

   /*
     * The bottom components of be may be reconstructed using the formula
     *      ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
     *      ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
     */
 
 
     dst[2][0][0] -= proj_matvec[0][0][1];
     dst[2][0][1] += proj_matvec[0][0][0];
     dst[3][0][0] += proj_matvec[1][0][1];
     dst[3][0][1] -= proj_matvec[1][0][0];

     dst[2][1][0] -= proj_matvec[0][1][1];
     dst[2][1][1] += proj_matvec[0][1][0];
     dst[3][1][0] += proj_matvec[1][1][1];
     dst[3][1][1] -= proj_matvec[1][1][0];

     dst[2][2][0] -= proj_matvec[0][2][1];
     dst[2][2][1] += proj_matvec[0][2][0];
     dst[3][2][0] += proj_matvec[1][2][1];
     dst[3][2][1] -= proj_matvec[1][2][0];

  
}

inline void mvv_recons_gamma2_plus_add_store( HalfSpinor src, 
				       GaugeMatrix u,
				       FourSpinor sum,
				       FourSpinor dst)
{
  HalfSpinor proj_matvec;
  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )
   
   * Therefore the top components are
     
   *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
   *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
   
   */
  su3_mult(proj_matvec, u, src);

  dst[0][0][0] = sum[0][0][0] + proj_matvec[0][0][0];
  dst[0][0][1] = sum[0][0][1] + proj_matvec[0][0][1];
  dst[0][1][0] = sum[0][1][0] + proj_matvec[0][1][0];
  dst[0][1][1] = sum[0][1][1] + proj_matvec[0][1][1];  
  dst[0][2][0] = sum[0][2][0] + proj_matvec[0][2][0];
  dst[0][2][1] = sum[0][2][1] + proj_matvec[0][2][1];

  dst[1][0][0] = sum[1][0][0] + proj_matvec[1][0][0];
  dst[1][0][1] = sum[1][0][1] + proj_matvec[1][0][1];
  dst[1][1][0] = sum[1][1][0] + proj_matvec[1][1][0];
  dst[1][1][1] = sum[1][1][1] + proj_matvec[1][1][1];
  dst[1][2][0] = sum[1][2][0] + proj_matvec[1][2][0];
  dst[1][2][1] = sum[1][2][1] + proj_matvec[1][2][1];
   
   /*
     * The bottom components of be may be reconstructed using the formula
     *      ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
     *      ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
     */
 
  dst[2][0][0] = sum[2][0][0] - proj_matvec[0][0][1];
  dst[2][0][1] = sum[2][0][1] + proj_matvec[0][0][0];
  dst[3][0][0] = sum[3][0][0] + proj_matvec[1][0][1];
  dst[3][0][1] = sum[3][0][1] - proj_matvec[1][0][0];
  
  dst[2][1][0] = sum[2][1][0] - proj_matvec[0][1][1];
  dst[2][1][1] = sum[2][1][1] + proj_matvec[0][1][0];
  dst[3][1][0] = sum[3][1][0] + proj_matvec[1][1][1];
  dst[3][1][1] = sum[3][1][1] - proj_matvec[1][1][0];
  
  dst[2][2][0] = sum[2][2][0] - proj_matvec[0][2][1];
  dst[2][2][1] = sum[2][2][1] + proj_matvec[0][2][0];
  dst[3][2][0] = sum[3][2][0] + proj_matvec[1][2][1];
  dst[3][2][1] = sum[3][2][1] -proj_matvec[1][2][0];
 
}




inline void mvv_recons_gamma3_plus_add_store( HalfSpinor src, 
			     GaugeMatrix u,
			     FourSpinor sum,
			    FourSpinor dst)
{
  HalfSpinor proj_matvec;

    /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
     *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
     *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r - a2r} + i{a0i - a2i} )
     *      ( b1r + i b1i )     ( {a1r - a3r} + i{a1i - a3i} )
     */

  su3_mult(proj_matvec, u, src);

  dst[0][0][0] = sum[0][0][0] + proj_matvec[0][0][0];
  dst[0][0][1] = sum[0][0][1] + proj_matvec[0][0][1];
  dst[0][1][0] = sum[0][1][0] + proj_matvec[0][1][0];
  dst[0][1][1] = sum[0][1][1] + proj_matvec[0][1][1];  
  dst[0][2][0] = sum[0][2][0] + proj_matvec[0][2][0];
  dst[0][2][1] = sum[0][2][1] + proj_matvec[0][2][1];

  dst[1][0][0] = sum[1][0][0] + proj_matvec[1][0][0];
  dst[1][0][1] = sum[1][0][1] + proj_matvec[1][0][1];
  dst[1][1][0] = sum[1][1][0] + proj_matvec[1][1][0];
  dst[1][1][1] = sum[1][1][1] + proj_matvec[1][1][1];
  dst[1][2][0] = sum[1][2][0] + proj_matvec[1][2][0];
  dst[1][2][1] = sum[1][2][1] + proj_matvec[1][2][1];

   
  /*
   *      ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
   *      ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
   */
  
  dst[2][0][0] = sum[2][0][0] - proj_matvec[0][0][0];
  dst[2][0][1] = sum[2][0][1] - proj_matvec[0][0][1];
  dst[3][0][0] = sum[3][0][0] - proj_matvec[1][0][0];
  dst[3][0][1] = sum[3][0][1] - proj_matvec[1][0][1];
  
  dst[2][1][0] = sum[2][1][0] - proj_matvec[0][1][0];
  dst[2][1][1] = sum[2][1][1] - proj_matvec[0][1][1];
  dst[3][1][0] = sum[3][1][0] - proj_matvec[1][1][0];
  dst[3][1][1] = sum[3][1][1] - proj_matvec[1][1][1];
  
  dst[2][2][0] = sum[2][2][0] - proj_matvec[0][2][0];
  dst[2][2][1] = sum[2][2][1] - proj_matvec[0][2][1];
  dst[3][2][0] = sum[3][2][0] - proj_matvec[1][2][0];
  dst[3][2][1] = sum[3][2][1] - proj_matvec[1][2][1];


}



inline void mvv_recons_gamma0_minus( HalfSpinor src, 
			     GaugeMatrix u,
			    FourSpinor dst)
{
  su3_mult(dst, u, src);
  

  /*                              ( 1  0  0 +i)  ( a0 )    ( a0 + i a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1 +i  0)  ( a1 )  = ( a1 + i a2 )
   *                    0         ( 0 -i  1  0)  ( a2 )    ( a2 - i a1 )
   *                              (-i  0  0  1)  ( a3 )    ( a3 - i a0 )
   */
  
  /* Bottom is:
   *      ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
   *      ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
   */
   
  dst[2][0][0] = dst[1][0][1];
  dst[2][0][1] = -dst[1][0][0];
  dst[3][0][0] = dst[0][0][1];
  dst[3][0][1] = -dst[0][0][0];
  
  dst[2][1][0] = dst[1][1][1];
  dst[2][1][1] = -dst[1][1][0];
  dst[3][1][0] = dst[0][1][1];
  dst[3][1][1] = -dst[0][1][0];
  
  dst[2][2][0] = dst[1][2][1];
  dst[2][2][1] = -dst[1][2][0];
  dst[3][2][0] = dst[0][2][1];
  dst[3][2][1] = -dst[0][2][0];


}

inline void mvv_recons_gamma1_minus_add( HalfSpinor src, 
				 GaugeMatrix u,
				FourSpinor dst)
{
  HalfSpinor proj_matvec;
  su3_mult(proj_matvec, u, src);
  /*                              ( 1  0  0 -1)  ( a0 )    ( a0 - a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  1  0)  ( a1 )  = ( a1 + a2 )
   *                    1         ( 0  1  1  0)  ( a2 )    ( a2 + a1 )
   *                              (-1  0  0  1)  ( a3 )    ( a3 - a0 )
   */

    dst[0][0][0] += proj_matvec[0][0][0];
    dst[0][0][1] += proj_matvec[0][0][1];
    dst[0][1][0] += proj_matvec[0][1][0];
    dst[0][1][1] += proj_matvec[0][1][1];  
    dst[0][2][0] += proj_matvec[0][2][0];
    dst[0][2][1] += proj_matvec[0][2][1];

    dst[1][0][0] += proj_matvec[1][0][0];
    dst[1][0][1] += proj_matvec[1][0][1];
    dst[1][1][0] += proj_matvec[1][1][0];
    dst[1][1][1] += proj_matvec[1][1][1];
    dst[1][2][0] += proj_matvec[1][2][0];
    dst[1][2][1] += proj_matvec[1][2][1];
   
    /*
     * ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
     * ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
     */
 
    
    dst[2][0][0] += proj_matvec[1][0][0];
    dst[2][0][1] += proj_matvec[1][0][1];
    dst[3][0][0] -= proj_matvec[0][0][0];
    dst[3][0][1] -= proj_matvec[0][0][1];
    
    dst[2][1][0] += proj_matvec[1][1][0];
    dst[2][1][1] += proj_matvec[1][1][1];
    dst[3][1][0] -= proj_matvec[0][1][0];
    dst[3][1][1] -= proj_matvec[0][1][1];
    
    dst[2][2][0] += proj_matvec[1][2][0];
    dst[2][2][1] += proj_matvec[1][2][1];
    dst[3][2][0] -= proj_matvec[0][2][0];
    dst[3][2][1] -= proj_matvec[0][2][1];
     
}


inline void mvv_recons_gamma2_minus_add( HalfSpinor src, 
				 GaugeMatrix u,
				FourSpinor dst)
{
  HalfSpinor proj_matvec;
  su3_mult(proj_matvec, u, src);


    /*                              ( 1  0  i  0)  ( a0 )    ( a0 + i a2 )
     *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0 -i)  ( a1 )  = ( a1 - i a3 )
     *                    2         (-i  0  1  0)  ( a2 )    ( a2 - i a0 )
     *                              ( 0  i  0  1)  ( a3 )    ( a3 + i a1 )
     */
    /*
     * Therefore the top components are
     *      ( b0r + i b0i )  =  ( {a0r - a2i} + i{a0i + a2r} )
     *      ( b1r + i b1i )     ( {a1r + a3i} + i{a1i - a3r} )
     */

    dst[0][0][0] += proj_matvec[0][0][0];
    dst[0][0][1] += proj_matvec[0][0][1];
    dst[0][1][0] += proj_matvec[0][1][0];
    dst[0][1][1] += proj_matvec[0][1][1];  
    dst[0][2][0] += proj_matvec[0][2][0];
    dst[0][2][1] += proj_matvec[0][2][1];

    dst[1][0][0] += proj_matvec[1][0][0];
    dst[1][0][1] += proj_matvec[1][0][1];
    dst[1][1][0] += proj_matvec[1][1][0];
    dst[1][1][1] += proj_matvec[1][1][1];
    dst[1][2][0] += proj_matvec[1][2][0];
    dst[1][2][1] += proj_matvec[1][2][1];
   
    /*
     *      ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
     *      ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
     */

     dst[2][0][0] += proj_matvec[0][0][1];
     dst[2][0][1] -= proj_matvec[0][0][0];
     dst[3][0][0] -= proj_matvec[1][0][1];
     dst[3][0][1] += proj_matvec[1][0][0];

     dst[2][1][0] += proj_matvec[0][1][1];
     dst[2][1][1] -= proj_matvec[0][1][0];
     dst[3][1][0] -= proj_matvec[1][1][1];
     dst[3][1][1] += proj_matvec[1][1][0];

     dst[2][2][0] += proj_matvec[0][2][1];
     dst[2][2][1] -= proj_matvec[0][2][0];
     dst[3][2][0] -= proj_matvec[1][2][1];
     dst[3][2][1] += proj_matvec[1][2][0];

}


inline void mvv_recons_gamma2_minus_add_store( HalfSpinor src, 
					GaugeMatrix u,
					FourSpinor sum,
					FourSpinor dst)
{
  HalfSpinor proj_matvec;
  su3_mult(proj_matvec, u, src);
  
  /*                              ( 1  0  i  0)  ( a0 )    ( a0 + i a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0 -i)  ( a1 )  = ( a1 - i a3 )
   *                    2         (-i  0  1  0)  ( a2 )    ( a2 - i a0 )
   *                              ( 0  i  0  1)  ( a3 )    ( a3 + i a1 )
   */
  /*
   * Therefore the top components are
   *      ( b0r + i b0i )  =  ( {a0r - a2i} + i{a0i + a2r} )
   *      ( b1r + i b1i )     ( {a1r + a3i} + i{a1i - a3r} )
   */
  
  dst[0][0][0] = sum[0][0][0] + proj_matvec[0][0][0];
  dst[0][0][1] = sum[0][0][1] + proj_matvec[0][0][1];
  dst[0][1][0] = sum[0][1][0] + proj_matvec[0][1][0];
  dst[0][1][1] = sum[0][1][1] + proj_matvec[0][1][1];  
  dst[0][2][0] = sum[0][2][0] + proj_matvec[0][2][0];
  dst[0][2][1] = sum[0][2][1] + proj_matvec[0][2][1];

  dst[1][0][0] = sum[1][0][0] + proj_matvec[1][0][0];
  dst[1][0][1] = sum[1][0][1] + proj_matvec[1][0][1];
  dst[1][1][0] = sum[1][1][0] + proj_matvec[1][1][0];
  dst[1][1][1] = sum[1][1][1] + proj_matvec[1][1][1];
  dst[1][2][0] = sum[1][2][0] + proj_matvec[1][2][0];
  dst[1][2][1] = sum[1][2][1] + proj_matvec[1][2][1];
   
  /*
   *   ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *   ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
   */

  dst[2][0][0] = sum[2][0][0] + proj_matvec[0][0][1];
  dst[2][0][1] = sum[2][0][1] - proj_matvec[0][0][0];
  dst[3][0][0] = sum[3][0][0] - proj_matvec[1][0][1];
  dst[3][0][1] = sum[3][0][1] + proj_matvec[1][0][0];
  
  dst[2][1][0] = sum[2][1][0] + proj_matvec[0][1][1];
  dst[2][1][1] = sum[2][1][1] - proj_matvec[0][1][0];
  dst[3][1][0] = sum[3][1][0] - proj_matvec[1][1][1];
  dst[3][1][1] = sum[3][1][1] + proj_matvec[1][1][0];
  
  dst[2][2][0] = sum[2][2][0] + proj_matvec[0][2][1];
  dst[2][2][1] = sum[2][2][1] - proj_matvec[0][2][0];
  dst[3][2][0] = sum[3][2][0] - proj_matvec[1][2][1];
  dst[3][2][1] = sum[3][2][1] + proj_matvec[1][2][0];

}




inline void mvv_recons_gamma3_minus_add_store( HalfSpinor src, 
			     GaugeMatrix u,
			     FourSpinor sum,
			    FourSpinor dst)
{
  HalfSpinor proj_matvec;
  su3_mult(proj_matvec, u, src);

  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
   *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
   *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
   
   * Therefore the top components are
   
   *      ( b0r + i b0i )  =  ( {a0r + a2r} + i{a0i + a2i} )
   *      ( b1r + i b1i )     ( {a1r + a3r} + i{a1i + a3i} )
   */

  dst[0][0][0] = sum[0][0][0] + proj_matvec[0][0][0];
  dst[0][0][1] = sum[0][0][1] + proj_matvec[0][0][1];
  dst[0][1][0] = sum[0][1][0] + proj_matvec[0][1][0];
  dst[0][1][1] = sum[0][1][1] + proj_matvec[0][1][1];  
  dst[0][2][0] = sum[0][2][0] + proj_matvec[0][2][0];
  dst[0][2][1] = sum[0][2][1] + proj_matvec[0][2][1];

  dst[1][0][0] = sum[1][0][0] + proj_matvec[1][0][0];
  dst[1][0][1] = sum[1][0][1] + proj_matvec[1][0][1];
  dst[1][1][0] = sum[1][1][0] + proj_matvec[1][1][0];
  dst[1][1][1] = sum[1][1][1] + proj_matvec[1][1][1];
  dst[1][2][0] = sum[1][2][0] + proj_matvec[1][2][0];
  dst[1][2][1] = sum[1][2][1] + proj_matvec[1][2][1];
 
    /*
     *      ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
     *      ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
     */
     dst[2][0][0] = sum[2][0][0] + proj_matvec[0][0][0];
     dst[2][0][1] = sum[2][0][1] + proj_matvec[0][0][1];
     dst[3][0][0] = sum[3][0][0] + proj_matvec[1][0][0];
     dst[3][0][1] = sum[3][0][1] + proj_matvec[1][0][1];

     dst[2][1][0] = sum[2][1][0] + proj_matvec[0][1][0];
     dst[2][1][1] = sum[2][1][1] + proj_matvec[0][1][1];
     dst[3][1][0] = sum[3][1][0] + proj_matvec[1][1][0];
     dst[3][1][1] = sum[3][1][1] + proj_matvec[1][1][1];

     dst[2][2][0] = sum[2][2][0] + proj_matvec[0][2][0];
     dst[2][2][1] = sum[2][2][1] + proj_matvec[0][2][1];
     dst[3][2][0] = sum[3][2][0] + proj_matvec[1][2][0];
     dst[3][2][1] = sum[3][2][1] + proj_matvec[1][2][1];

}


  }
}

#endif
