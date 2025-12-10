#ifndef CPP_DSLASH_SCALAR_32BIT_C_H
#define CPP_DSLASH_SCALAR_32BIT_C_H


#warning "using C stuff"

#ifndef ALIGN
#define ALIGN __attribute__ ((aligned(16)))
#endif

#include <cpp_dslash_types.h>
#include <cpp_dslash_matvec32bit_c.h>

namespace CPlusPlusWilsonDslash {

  namespace DslashScalar32Bit {

       /* Thread dispatch function for D */
    void DPsiPlus(size_t lo, size_t hi, int id, const void *ptr);
    
    /* Thread dispatch function for D^\dagger */
    void DPsiMinus(size_t lo, size_t hi, int id, const void *ptr);

       /* Thread dispatch function for D */
    void DPsiPlus3D(size_t lo, size_t hi, int id, const void *ptr);
    
    /* Thread dispatch function for D^\dagger */
    void DPsiMinus3D(size_t lo, size_t hi, int id, const void *ptr);

    using namespace Dslash32BitTypes;


    inline
    void dslash_plus_dir0_forward(FourSpinor spinor_in,
				  GaugeMatrix u,
				  HalfSpinor  upper_sum,
				  HalfSpinor  lower_sum   )
      
    {
      
      HalfSpinor proj_result;
       
      /* 1 - gamma_0 projection */


      /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
       *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
       *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
       *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )
       
       * Therefore the top components are
       
       *      ( b0r + i b0i )  =  ( {a0r + a3i} + i{a0i - a3r} )
       *      ( b1r + i b1i )     ( {a1r + a2i} + i{a1i - a2r} )
       */

      // Color = 0
     proj_result[0][0][0] = spinor_in[0][0][0] + spinor_in[3][0][1];
     proj_result[0][0][1] = spinor_in[0][0][1] - spinor_in[3][0][0];
     proj_result[0][1][0] = spinor_in[1][0][0] + spinor_in[2][0][1];
     proj_result[0][1][1] = spinor_in[1][0][1] - spinor_in[2][0][0];

     proj_result[1][0][0] = spinor_in[0][1][0] + spinor_in[3][1][1];
     proj_result[1][0][1] = spinor_in[0][1][1] - spinor_in[3][1][0];
     proj_result[1][1][0] = spinor_in[1][1][0] + spinor_in[2][1][1];
     proj_result[1][1][1] = spinor_in[1][1][1] - spinor_in[2][1][0];
    
     proj_result[2][0][0] = spinor_in[0][2][0] + spinor_in[3][2][1];
     proj_result[2][0][1] = spinor_in[0][2][1] - spinor_in[3][2][0];
     proj_result[2][1][0] = spinor_in[1][2][0] + spinor_in[2][2][1];
     proj_result[2][1][1] = spinor_in[1][2][1] - spinor_in[2][2][0];
     
     // su3_multiply
     su3_mult(upper_sum, u, proj_result);
     
     /*
      *      ( b2r + i b2i )  =  ( {a2r - a1i} + i{a2i + a1r} )  =  ( - b1i + i b1r )
      *      ( b3r + i b3i )     ( {a3r - a0i} + i{a3i + a0r} )     ( - b0i + i b0r ) 
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
  void dslash_plus_dir0_backward_add( FourSpinor spinor_in,
				      GaugeMatrix u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )

  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;


    /*                              ( 1  0  0 +i)  ( a0 )    ( a0 + i a3 )
     *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1 +i  0)  ( a1 )  = ( a1 + i a2 )
     *                    0         ( 0 -i  1  0)  ( a2 )    ( a2 - i a1 )
     *                              (-i  0  0  1)  ( a3 )    ( a3 - i a0 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r - a3i} + i{a0i + a3r} )
     *      ( b1r + i b1i )     ( {a1r - a2i} + i{a1i + a2r} )
     
     * The bottom components of be may be reconstructed using the formula
     
     *      ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
     *      ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
     */
    // Color = 0

    proj_result[0][0][0] = spinor_in[0][0][0] - spinor_in[3][0][1];
    proj_result[0][0][1] = spinor_in[0][0][1] + spinor_in[3][0][0];
    proj_result[0][1][0] = spinor_in[1][0][0] - spinor_in[2][0][1];
    proj_result[0][1][1] = spinor_in[1][0][1] + spinor_in[2][0][0];

    // Color = 1
    proj_result[1][0][0] = spinor_in[0][1][0] - spinor_in[3][1][1];
    proj_result[1][0][1] = spinor_in[0][1][1] + spinor_in[3][1][0];
    proj_result[1][1][0] = spinor_in[1][1][0] - spinor_in[2][1][1];
    proj_result[1][1][1] = spinor_in[1][1][1] + spinor_in[2][1][0];

    // Color = 2
    proj_result[2][0][0] = spinor_in[0][2][0] - spinor_in[3][2][1];
    proj_result[2][0][1] = spinor_in[0][2][1] + spinor_in[3][2][0];
    proj_result[2][1][0] = spinor_in[1][2][0] - spinor_in[2][2][1];
    proj_result[2][1][1] = spinor_in[1][2][1] + spinor_in[2][2][0];

    su3_adj_mult(proj_matvec, u, proj_result);

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
    
    /*
     *      ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
     *      ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
     */
    lower_sum[0][0][0] += proj_matvec[0][1][1];
    lower_sum[0][0][1] -= proj_matvec[0][1][0];
    lower_sum[0][1][0] += proj_matvec[0][0][1];
    lower_sum[0][1][1] -= proj_matvec[0][0][0];
    
    
    lower_sum[1][0][0] += proj_matvec[1][1][1];
    lower_sum[1][0][1] -= proj_matvec[1][1][0];
    lower_sum[1][1][0] += proj_matvec[1][0][1];
    lower_sum[1][1][1] -= proj_matvec[1][0][0];
    
    
    lower_sum[2][0][0] += proj_matvec[2][1][1];
    lower_sum[2][0][1] -= proj_matvec[2][1][0];
    lower_sum[2][1][0] += proj_matvec[2][0][1];
    lower_sum[2][1][1] -= proj_matvec[2][0][0];
    


  }

  inline
  void dslash_plus_dir1_forward_add( FourSpinor  spinor_in,
				     GaugeMatrix  u,
				     HalfSpinor  upper_sum,
				     HalfSpinor  lower_sum   )
  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;


  /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
   *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
   *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
	 
   * Therefore the top components are
      
   *      ( b0r + i b0i )  =  ( {a0r + a3r} + i{a0i + a3i} )
   *      ( b1r + i b1i )     ( {a1r - a2r} + i{a1i - a2i} )
   */
    // color =0
    proj_result[0][0][0] = spinor_in[0][0][0] + spinor_in[3][0][0];
    proj_result[0][0][1] = spinor_in[0][0][1] + spinor_in[3][0][1];
    proj_result[0][1][0] = spinor_in[1][0][0] - spinor_in[2][0][0];
    proj_result[0][1][1] = spinor_in[1][0][1] - spinor_in[2][0][1];

    // color = 1
    proj_result[1][0][0] = spinor_in[0][1][0] + spinor_in[3][1][0];
    proj_result[1][0][1] = spinor_in[0][1][1] + spinor_in[3][1][1];
    proj_result[1][1][0] = spinor_in[1][1][0] - spinor_in[2][1][0];
    proj_result[1][1][1] = spinor_in[1][1][1] - spinor_in[2][1][1];

    // color = 2
    proj_result[2][0][0] = spinor_in[0][2][0] + spinor_in[3][2][0];
    proj_result[2][0][1] = spinor_in[0][2][1] + spinor_in[3][2][1];
    proj_result[2][1][0] = spinor_in[1][2][0] - spinor_in[2][2][0];
    proj_result[2][1][1] = spinor_in[1][2][1] - spinor_in[2][2][1];


    su3_mult(proj_matvec, u, proj_result);

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
    
    /*
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
  void dslash_plus_dir1_backward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )
  {

    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

    // Color = 0

  /*                              ( 1  0  0 -1)  ( a0 )    ( a0 - a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  1  0)  ( a1 )  = ( a1 + a2 )
   *                    1         ( 0  1  1  0)  ( a2 )    ( a2 + a1 )
   *                              (-1  0  0  1)  ( a3 )    ( a3 - a0 )

   * Therefore the top components are

   *      ( b0r + i b0i )  =  ( {a0r - a3r} + i{a0i - a3i} )
   *      ( b1r + i b1i )     ( {a1r + a2r} + i{a1i + a2i} )
   */

    proj_result[0][0][0] = spinor_in[0][0][0] - spinor_in[3][0][0];
    proj_result[0][0][1] = spinor_in[0][0][1] - spinor_in[3][0][1];
    proj_result[0][1][0] = spinor_in[1][0][0] + spinor_in[2][0][0];
    proj_result[0][1][1] = spinor_in[1][0][1] + spinor_in[2][0][1];


    proj_result[1][0][0] = spinor_in[0][1][0] - spinor_in[3][1][0];
    proj_result[1][0][1] = spinor_in[0][1][1] - spinor_in[3][1][1];
    proj_result[1][1][0] = spinor_in[1][1][0] + spinor_in[2][1][0];
    proj_result[1][1][1] = spinor_in[1][1][1] + spinor_in[2][1][1];


    proj_result[2][0][0] = spinor_in[0][2][0] - spinor_in[3][2][0];
    proj_result[2][0][1] = spinor_in[0][2][1] - spinor_in[3][2][1];
    proj_result[2][1][0] = spinor_in[1][2][0] + spinor_in[2][2][0];
    proj_result[2][1][1] = spinor_in[1][2][1] + spinor_in[2][2][1];


    su3_adj_mult(proj_matvec, u, proj_result);

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
    
    /*
     *      ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
     *      ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
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
  void dslash_plus_dir2_forward_add( FourSpinor  spinor_in,
				     GaugeMatrix  u,
				     HalfSpinor  upper_sum,
				     HalfSpinor  lower_sum   )
  {

    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )

   * Therefore the top components are
      
   *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
   *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
      
   */
    proj_result[0][0][0] = spinor_in[0][0][0] + spinor_in[2][0][1];
    proj_result[0][0][1] = spinor_in[0][0][1] - spinor_in[2][0][0];
    proj_result[0][1][0] = spinor_in[1][0][0] - spinor_in[3][0][1];
    proj_result[0][1][1] = spinor_in[1][0][1] + spinor_in[3][0][0];


    proj_result[1][0][0] = spinor_in[0][1][0] + spinor_in[2][1][1];
    proj_result[1][0][1] = spinor_in[0][1][1] - spinor_in[2][1][0];
    proj_result[1][1][0] = spinor_in[1][1][0] - spinor_in[3][1][1];
    proj_result[1][1][1] = spinor_in[1][1][1] + spinor_in[3][1][0];


    proj_result[2][0][0] = spinor_in[0][2][0] + spinor_in[2][2][1];
    proj_result[2][0][1] = spinor_in[0][2][1] - spinor_in[2][2][0];
    proj_result[2][1][0] = spinor_in[1][2][0] - spinor_in[3][2][1];
    proj_result[2][1][1] = spinor_in[1][2][1] + spinor_in[3][2][0];


    su3_mult(proj_matvec, u, proj_result);


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

    /*
     * The bottom components of be may be reconstructed using the formula
     *      ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
     *      ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
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
  void dslash_plus_dir2_backward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )
  {

    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

    
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

    proj_result[0][0][0] = spinor_in[0][0][0] - spinor_in[2][0][1];
    proj_result[0][0][1] = spinor_in[0][0][1] + spinor_in[2][0][0];
    proj_result[0][1][0] = spinor_in[1][0][0] + spinor_in[3][0][1];
    proj_result[0][1][1] = spinor_in[1][0][1] - spinor_in[3][0][0];

    proj_result[1][0][0] = spinor_in[0][1][0] - spinor_in[2][1][1];
    proj_result[1][0][1] = spinor_in[0][1][1] + spinor_in[2][1][0];
    proj_result[1][1][0] = spinor_in[1][1][0] + spinor_in[3][1][1];
    proj_result[1][1][1] = spinor_in[1][1][1] - spinor_in[3][1][0];

    proj_result[2][0][0] = spinor_in[0][2][0] - spinor_in[2][2][1];
    proj_result[2][0][1] = spinor_in[0][2][1] + spinor_in[2][2][0];
    proj_result[2][1][0] = spinor_in[1][2][0] + spinor_in[3][2][1];
    proj_result[2][1][1] = spinor_in[1][2][1] - spinor_in[3][2][0];
    
    su3_adj_mult(proj_matvec, u, proj_result);

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

    /*
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
  void dslash_plus_dir2_backward_add_store( FourSpinor  spinor_in,
					    GaugeMatrix  u,
					    HalfSpinor  upper_sum,
					    HalfSpinor  lower_sum,
					    FourSpinor spinor_out)
  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

    
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

    proj_result[0][0][0] = spinor_in[0][0][0] - spinor_in[2][0][1];
    proj_result[0][0][1] = spinor_in[0][0][1] + spinor_in[2][0][0];
    proj_result[0][1][0] = spinor_in[1][0][0] + spinor_in[3][0][1];
    proj_result[0][1][1] = spinor_in[1][0][1] - spinor_in[3][0][0];

    proj_result[1][0][0] = spinor_in[0][1][0] - spinor_in[2][1][1];
    proj_result[1][0][1] = spinor_in[0][1][1] + spinor_in[2][1][0];
    proj_result[1][1][0] = spinor_in[1][1][0] + spinor_in[3][1][1];
    proj_result[1][1][1] = spinor_in[1][1][1] - spinor_in[3][1][0];

    proj_result[2][0][0] = spinor_in[0][2][0] - spinor_in[2][2][1];
    proj_result[2][0][1] = spinor_in[0][2][1] + spinor_in[2][2][0];
    proj_result[2][1][0] = spinor_in[1][2][0] + spinor_in[3][2][1];
    proj_result[2][1][1] = spinor_in[1][2][1] - spinor_in[3][2][0];

    su3_adj_mult(proj_matvec, u, proj_result);

    // Top half
    spinor_out[0][0][0] = upper_sum[0][0][0] + proj_matvec[0][0][0];
    spinor_out[0][0][1] = upper_sum[0][0][1] + proj_matvec[0][0][1];
    spinor_out[0][1][0] = upper_sum[1][0][0] + proj_matvec[1][0][0];
    spinor_out[0][1][1] = upper_sum[1][0][1] + proj_matvec[1][0][1];
    spinor_out[0][2][0] = upper_sum[2][0][0] + proj_matvec[2][0][0];
    spinor_out[0][2][1] = upper_sum[2][0][1] + proj_matvec[2][0][1];

    spinor_out[1][0][0] = upper_sum[0][1][0] + proj_matvec[0][1][0];
    spinor_out[1][0][1] = upper_sum[0][1][1] + proj_matvec[0][1][1];
    spinor_out[1][1][0] = upper_sum[1][1][0] + proj_matvec[1][1][0];
    spinor_out[1][1][1] = upper_sum[1][1][1] + proj_matvec[1][1][1];
    spinor_out[1][2][0] = upper_sum[2][1][0] + proj_matvec[2][1][0];
    spinor_out[1][2][1] = upper_sum[2][1][1] + proj_matvec[2][1][1];
    
    /*
     *      ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
     *      ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
     */

    // Arrange for linear storage order 
    spinor_out[2][0][0] = lower_sum[0][0][0] + proj_matvec[0][0][1];
    spinor_out[2][0][1] = lower_sum[0][0][1] - proj_matvec[0][0][0];
    spinor_out[2][1][0] = lower_sum[1][0][0] + proj_matvec[1][0][1];
    spinor_out[2][1][1] = lower_sum[1][0][1] - proj_matvec[1][0][0];
    spinor_out[2][2][0] = lower_sum[2][0][0] + proj_matvec[2][0][1];
    spinor_out[2][2][1] = lower_sum[2][0][1] - proj_matvec[2][0][0];

    spinor_out[3][0][0] = lower_sum[0][1][0] - proj_matvec[0][1][1];
    spinor_out[3][0][1] = lower_sum[0][1][1] + proj_matvec[0][1][0];
    spinor_out[3][1][0] = lower_sum[1][1][0] - proj_matvec[1][1][1];
    spinor_out[3][1][1] = lower_sum[1][1][1] + proj_matvec[1][1][0];
    spinor_out[3][2][0] = lower_sum[2][1][0] - proj_matvec[2][1][1];
    spinor_out[3][2][1] = lower_sum[2][1][1] + proj_matvec[2][1][0];

    
    
  }


  inline
  void dslash_plus_dir3_forward_add( FourSpinor  spinor_in,
				     GaugeMatrix  u,
				     HalfSpinor  upper_sum,
				     HalfSpinor  lower_sum)
  {

    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

    /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
     *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
     *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r - a2r} + i{a0i - a2i} )
     *      ( b1r + i b1i )     ( {a1r - a3r} + i{a1i - a3i} )
     */

    proj_result[0][0][0] = spinor_in[0][0][0] - spinor_in[2][0][0];
    proj_result[0][0][1] = spinor_in[0][0][1] - spinor_in[2][0][1];
    proj_result[0][1][0] = spinor_in[1][0][0] - spinor_in[3][0][0];
    proj_result[0][1][1] = spinor_in[1][0][1] - spinor_in[3][0][1];

    proj_result[1][0][0] = spinor_in[0][1][0] - spinor_in[2][1][0];
    proj_result[1][0][1] = spinor_in[0][1][1] - spinor_in[2][1][1];
    proj_result[1][1][0] = spinor_in[1][1][0] - spinor_in[3][1][0];
    proj_result[1][1][1] = spinor_in[1][1][1] - spinor_in[3][1][1];

    proj_result[2][0][0] = spinor_in[0][2][0] - spinor_in[2][2][0];
    proj_result[2][0][1] = spinor_in[0][2][1] - spinor_in[2][2][1];
    proj_result[2][1][0] = spinor_in[1][2][0] - spinor_in[3][2][0];
    proj_result[2][1][1] = spinor_in[1][2][1] - spinor_in[3][2][1];


    su3_mult(proj_matvec, u, proj_result);

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


    /*
     *      ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
     *      ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
     */

    lower_sum[0][0][0] -= proj_matvec[0][0][0];
    lower_sum[0][0][1] -= proj_matvec[0][0][1];
    lower_sum[0][1][0] -= proj_matvec[0][1][0];
    lower_sum[0][1][1] -= proj_matvec[0][1][1];

    lower_sum[1][0][0] -= proj_matvec[1][0][0];
    lower_sum[1][0][1] -= proj_matvec[1][0][1];
    lower_sum[1][1][0] -= proj_matvec[1][1][0];
    lower_sum[1][1][1] -= proj_matvec[1][1][1];

    lower_sum[2][0][0] -= proj_matvec[2][0][0];
    lower_sum[2][0][1] -= proj_matvec[2][0][1];
    lower_sum[2][1][0] -= proj_matvec[2][1][0];
    lower_sum[2][1][1] -= proj_matvec[2][1][1];

  }

  inline
  void dslash_plus_dir3_backward_add_store( FourSpinor  spinor_in,
					    GaugeMatrix  u,
					    HalfSpinor  upper_sum,
					    HalfSpinor  lower_sum,
					    FourSpinor spinor_out)

  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

    /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
     *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
     *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
     *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r + a2r} + i{a0i + a2i} )
     *      ( b1r + i b1i )     ( {a1r + a3r} + i{a1i + a3i} )
     */

    proj_result[0][0][0] = spinor_in[0][0][0] + spinor_in[2][0][0];
    proj_result[0][0][1] = spinor_in[0][0][1] + spinor_in[2][0][1];
    proj_result[0][1][0] = spinor_in[1][0][0] + spinor_in[3][0][0];
    proj_result[0][1][1] = spinor_in[1][0][1] + spinor_in[3][0][1];

    proj_result[1][0][0] = spinor_in[0][1][0] + spinor_in[2][1][0];
    proj_result[1][0][1] = spinor_in[0][1][1] + spinor_in[2][1][1];
    proj_result[1][1][0] = spinor_in[1][1][0] + spinor_in[3][1][0];
    proj_result[1][1][1] = spinor_in[1][1][1] + spinor_in[3][1][1];

    proj_result[2][0][0] = spinor_in[0][2][0] + spinor_in[2][2][0];
    proj_result[2][0][1] = spinor_in[0][2][1] + spinor_in[2][2][1];
    proj_result[2][1][0] = spinor_in[1][2][0] + spinor_in[3][2][0];
    proj_result[2][1][1] = spinor_in[1][2][1] + spinor_in[3][2][1];


    su3_adj_mult(proj_matvec, u, proj_result);

    // Top half 
    // rearrange for linear storage
    spinor_out[0][0][0] = upper_sum[0][0][0] + proj_matvec[0][0][0];
    spinor_out[0][0][1] = upper_sum[0][0][1] + proj_matvec[0][0][1];
    spinor_out[0][1][0] = upper_sum[1][0][0] + proj_matvec[1][0][0];
    spinor_out[0][1][1] = upper_sum[1][0][1] + proj_matvec[1][0][1];
    spinor_out[0][2][0] = upper_sum[2][0][0] + proj_matvec[2][0][0];
    spinor_out[0][2][1] = upper_sum[2][0][1] + proj_matvec[2][0][1];

    spinor_out[1][0][0] = upper_sum[0][1][0] + proj_matvec[0][1][0];
    spinor_out[1][0][1] = upper_sum[0][1][1] + proj_matvec[0][1][1];
    spinor_out[1][1][0] = upper_sum[1][1][0] + proj_matvec[1][1][0];
    spinor_out[1][1][1] = upper_sum[1][1][1] + proj_matvec[1][1][1];
    spinor_out[1][2][0] = upper_sum[2][1][0] + proj_matvec[2][1][0];
    spinor_out[1][2][1] = upper_sum[2][1][1] + proj_matvec[2][1][1];
 
    /*
     *      ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
     *      ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
     */
    
    spinor_out[2][0][0] = lower_sum[0][0][0]+ proj_matvec[0][0][0];
    spinor_out[2][0][1] = lower_sum[0][0][1]+ proj_matvec[0][0][1];
    spinor_out[2][1][0] = lower_sum[1][0][0]+ proj_matvec[1][0][0];
    spinor_out[2][1][1] = lower_sum[1][0][1]+ proj_matvec[1][0][1];
    spinor_out[2][2][0] = lower_sum[2][0][0]+ proj_matvec[2][0][0];
    spinor_out[2][2][1] = lower_sum[2][0][1]+ proj_matvec[2][0][1];

    spinor_out[3][0][0] = lower_sum[0][1][0]+ proj_matvec[0][1][0];
    spinor_out[3][0][1] = lower_sum[0][1][1]+ proj_matvec[0][1][1];
    spinor_out[3][1][0] = lower_sum[1][1][0]+ proj_matvec[1][1][0];
    spinor_out[3][1][1] = lower_sum[1][1][1]+ proj_matvec[1][1][1];
    spinor_out[3][2][0] = lower_sum[2][1][0]+ proj_matvec[2][1][0];
    spinor_out[3][2][1] = lower_sum[2][1][1]+ proj_matvec[2][1][1];
  }


  inline
  void dslash_minus_dir0_forward( FourSpinor spinor_in,
				  GaugeMatrix u,
				  HalfSpinor  upper_sum,
				  HalfSpinor  lower_sum   )
       
  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

    
  /*                              ( 1  0  0 +i)  ( a0 )    ( a0 + i a3 )
     *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1 +i  0)  ( a1 )  = ( a1 + i a2 )
     *                    0         ( 0 -i  1  0)  ( a2 )    ( a2 - i a1 )
     *                              (-i  0  0  1)  ( a3 )    ( a3 - i a0 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r - a3i} + i{a0i + a3r} )
     *      ( b1r + i b1i )     ( {a1r - a2i} + i{a1i + a2r} )
     
     * The bottom components of be may be reconstructed using the formula
     
     *      ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
     *      ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
     */
    // Color = 0

    proj_result[0][0][0] = spinor_in[0][0][0] - spinor_in[3][0][1];
    proj_result[0][0][1] = spinor_in[0][0][1] + spinor_in[3][0][0];
    proj_result[0][1][0] = spinor_in[1][0][0] - spinor_in[2][0][1];
    proj_result[0][1][1] = spinor_in[1][0][1] + spinor_in[2][0][0];

    // Color = 1
    proj_result[1][0][0] = spinor_in[0][1][0] - spinor_in[3][1][1];
    proj_result[1][0][1] = spinor_in[0][1][1] + spinor_in[3][1][0];
    proj_result[1][1][0] = spinor_in[1][1][0] - spinor_in[2][1][1];
    proj_result[1][1][1] = spinor_in[1][1][1] + spinor_in[2][1][0];

    // Color = 2
    proj_result[2][0][0] = spinor_in[0][2][0] - spinor_in[3][2][1];
    proj_result[2][0][1] = spinor_in[0][2][1] + spinor_in[3][2][0];
    proj_result[2][1][0] = spinor_in[1][2][0] - spinor_in[2][2][1];
    proj_result[2][1][1] = spinor_in[1][2][1] + spinor_in[2][2][0];

    su3_mult(upper_sum, u, proj_result);

    /*
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
  void dslash_minus_dir0_backward_add( FourSpinor spinor_in,
				       GaugeMatrix u,
				       HalfSpinor  upper_sum,
				       HalfSpinor  lower_sum   )

  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;
      /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
       *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
       *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
       *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )
       
       * Therefore the top components are
       
       *      ( b0r + i b0i )  =  ( {a0r + a3i} + i{a0i - a3r} )
       *      ( b1r + i b1i )     ( {a1r + a2i} + i{a1i - a2r} )
       */

    // Color = 0

    proj_result[0][0][0] = spinor_in[0][0][0] + spinor_in[3][0][1];
    proj_result[0][0][1] = spinor_in[0][0][1] - spinor_in[3][0][0];
    proj_result[0][1][0] = spinor_in[1][0][0] + spinor_in[2][0][1];
    proj_result[0][1][1] = spinor_in[1][0][1] - spinor_in[2][0][0];
    
    // Color = 1
    proj_result[1][0][0] = spinor_in[0][1][0] + spinor_in[3][1][1];
    proj_result[1][0][1] = spinor_in[0][1][1] - spinor_in[3][1][0];
    proj_result[1][1][0] = spinor_in[1][1][0] + spinor_in[2][1][1];
    proj_result[1][1][1] = spinor_in[1][1][1] - spinor_in[2][1][0];
    
    // Color = 2
    proj_result[2][0][0] = spinor_in[0][2][0] + spinor_in[3][2][1];
    proj_result[2][0][1] = spinor_in[0][2][1] - spinor_in[3][2][0];
    proj_result[2][1][0] = spinor_in[1][2][0] + spinor_in[2][2][1];
    proj_result[2][1][1] = spinor_in[1][2][1] - spinor_in[2][2][0];

    su3_adj_mult(proj_matvec, u, proj_result);

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
    
    /*
     *      ( b2r + i b2i )  =  ( {a2r - a1i} + i{a2i + a1r} )  =  ( - b1i + i b1r )
     *      ( b3r + i b3i )     ( {a3r - a0i} + i{a3i + a0r} )     ( - b0i + i b0r ) 
     */
    // col=0
    lower_sum[0][0][0] -= proj_matvec[0][1][1];
    lower_sum[0][0][1] += proj_matvec[0][1][0];
    lower_sum[0][1][0] -= proj_matvec[0][0][1];
    lower_sum[0][1][1] += proj_matvec[0][0][0];

    // col=1
    lower_sum[1][0][0] -= proj_matvec[1][1][1];
    lower_sum[1][0][1] += proj_matvec[1][1][0];
    lower_sum[1][1][0] -= proj_matvec[1][0][1];
    lower_sum[1][1][1] += proj_matvec[1][0][0];

    // col=2
    lower_sum[2][0][0] -= proj_matvec[2][1][1];
    lower_sum[2][0][1] += proj_matvec[2][1][0];
    lower_sum[2][1][0] -= proj_matvec[2][0][1];
    lower_sum[2][1][1] += proj_matvec[2][0][0];
  }


  inline
  void dslash_minus_dir1_forward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )
  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

      /*                              ( 1  0  0 -1)  ( a0 )    ( a0 - a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  1  0)  ( a1 )  = ( a1 + a2 )
   *                    1         ( 0  1  1  0)  ( a2 )    ( a2 + a1 )
   *                              (-1  0  0  1)  ( a3 )    ( a3 - a0 )

   * Therefore the top components are

   *      ( b0r + i b0i )  =  ( {a0r - a3r} + i{a0i - a3i} )
   *      ( b1r + i b1i )     ( {a1r + a2r} + i{a1i + a2i} )
   */

    proj_result[0][0][0] = spinor_in[0][0][0] - spinor_in[3][0][0];
    proj_result[0][0][1] = spinor_in[0][0][1] - spinor_in[3][0][1];
    proj_result[0][1][0] = spinor_in[1][0][0] + spinor_in[2][0][0];
    proj_result[0][1][1] = spinor_in[1][0][1] + spinor_in[2][0][1];


    proj_result[1][0][0] = spinor_in[0][1][0] - spinor_in[3][1][0];
    proj_result[1][0][1] = spinor_in[0][1][1] - spinor_in[3][1][1];
    proj_result[1][1][0] = spinor_in[1][1][0] + spinor_in[2][1][0];
    proj_result[1][1][1] = spinor_in[1][1][1] + spinor_in[2][1][1];


    proj_result[2][0][0] = spinor_in[0][2][0] - spinor_in[3][2][0];
    proj_result[2][0][1] = spinor_in[0][2][1] - spinor_in[3][2][1];
    proj_result[2][1][0] = spinor_in[1][2][0] + spinor_in[2][2][0];
    proj_result[2][1][1] = spinor_in[1][2][1] + spinor_in[2][2][1];

    su3_mult(proj_matvec, u, proj_result);

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
    
    /*
     *      ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
     *      ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
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
  void dslash_minus_dir1_backward_add( FourSpinor  spinor_in,
				       GaugeMatrix  u,
				       HalfSpinor  upper_sum,
				       HalfSpinor  lower_sum   )
  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;
 

  /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
   *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
   *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
	 
   * Therefore the top components are
      
   *      ( b0r + i b0i )  =  ( {a0r + a3r} + i{a0i + a3i} )
   *      ( b1r + i b1i )     ( {a1r - a2r} + i{a1i - a2i} )
   */
    // color =0
    proj_result[0][0][0] = spinor_in[0][0][0] + spinor_in[3][0][0];
    proj_result[0][0][1] = spinor_in[0][0][1] + spinor_in[3][0][1];
    proj_result[0][1][0] = spinor_in[1][0][0] - spinor_in[2][0][0];
    proj_result[0][1][1] = spinor_in[1][0][1] - spinor_in[2][0][1];

    // color = 1
    proj_result[1][0][0] = spinor_in[0][1][0] + spinor_in[3][1][0];
    proj_result[1][0][1] = spinor_in[0][1][1] + spinor_in[3][1][1];
    proj_result[1][1][0] = spinor_in[1][1][0] - spinor_in[2][1][0];
    proj_result[1][1][1] = spinor_in[1][1][1] - spinor_in[2][1][1];

    // color = 2
    proj_result[2][0][0] = spinor_in[0][2][0] + spinor_in[3][2][0];
    proj_result[2][0][1] = spinor_in[0][2][1] + spinor_in[3][2][1];
    proj_result[2][1][0] = spinor_in[1][2][0] - spinor_in[2][2][0];
    proj_result[2][1][1] = spinor_in[1][2][1] - spinor_in[2][2][1];


    su3_adj_mult(proj_matvec, u, proj_result);

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
    
    /*
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
  void dslash_minus_dir2_forward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum   )
  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

    
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

    proj_result[0][0][0] = spinor_in[0][0][0] - spinor_in[2][0][1];
    proj_result[0][0][1] = spinor_in[0][0][1] + spinor_in[2][0][0];
    proj_result[0][1][0] = spinor_in[1][0][0] + spinor_in[3][0][1];
    proj_result[0][1][1] = spinor_in[1][0][1] - spinor_in[3][0][0];

    proj_result[1][0][0] = spinor_in[0][1][0] - spinor_in[2][1][1];
    proj_result[1][0][1] = spinor_in[0][1][1] + spinor_in[2][1][0];
    proj_result[1][1][0] = spinor_in[1][1][0] + spinor_in[3][1][1];
    proj_result[1][1][1] = spinor_in[1][1][1] - spinor_in[3][1][0];

    proj_result[2][0][0] = spinor_in[0][2][0] - spinor_in[2][2][1];
    proj_result[2][0][1] = spinor_in[0][2][1] + spinor_in[2][2][0];
    proj_result[2][1][0] = spinor_in[1][2][0] + spinor_in[3][2][1];
    proj_result[2][1][1] = spinor_in[1][2][1] - spinor_in[3][2][0];
    
    su3_mult(proj_matvec, u, proj_result);

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

    /*
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
  void dslash_minus_dir2_backward_add( FourSpinor  spinor_in,
				       GaugeMatrix  u,
				       HalfSpinor  upper_sum,
				       HalfSpinor  lower_sum   )
  {
    
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )

   * Therefore the top components are
      
   *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
   *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
      
   */
    proj_result[0][0][0] = spinor_in[0][0][0] + spinor_in[2][0][1];
    proj_result[0][0][1] = spinor_in[0][0][1] - spinor_in[2][0][0];
    proj_result[0][1][0] = spinor_in[1][0][0] - spinor_in[3][0][1];
    proj_result[0][1][1] = spinor_in[1][0][1] + spinor_in[3][0][0];


    proj_result[1][0][0] = spinor_in[0][1][0] + spinor_in[2][1][1];
    proj_result[1][0][1] = spinor_in[0][1][1] - spinor_in[2][1][0];
    proj_result[1][1][0] = spinor_in[1][1][0] - spinor_in[3][1][1];
    proj_result[1][1][1] = spinor_in[1][1][1] + spinor_in[3][1][0];


    proj_result[2][0][0] = spinor_in[0][2][0] + spinor_in[2][2][1];
    proj_result[2][0][1] = spinor_in[0][2][1] - spinor_in[2][2][0];
    proj_result[2][1][0] = spinor_in[1][2][0] - spinor_in[3][2][1];
    proj_result[2][1][1] = spinor_in[1][2][1] + spinor_in[3][2][0];


    su3_adj_mult(proj_matvec, u, proj_result);


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

    /*
     * The bottom components of be may be reconstructed using the formula
     *      ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
     *      ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
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
  void dslash_minus_dir2_backward_add_store( FourSpinor  spinor_in,
					     GaugeMatrix  u,
					     HalfSpinor  upper_sum,
					     HalfSpinor  lower_sum,
					     FourSpinor spinor_out)
  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )

   * Therefore the top components are
      
   *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
   *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
      
   */
    proj_result[0][0][0] = spinor_in[0][0][0] + spinor_in[2][0][1];
    proj_result[0][0][1] = spinor_in[0][0][1] - spinor_in[2][0][0];
    proj_result[0][1][0] = spinor_in[1][0][0] - spinor_in[3][0][1];
    proj_result[0][1][1] = spinor_in[1][0][1] + spinor_in[3][0][0];


    proj_result[1][0][0] = spinor_in[0][1][0] + spinor_in[2][1][1];
    proj_result[1][0][1] = spinor_in[0][1][1] - spinor_in[2][1][0];
    proj_result[1][1][0] = spinor_in[1][1][0] - spinor_in[3][1][1];
    proj_result[1][1][1] = spinor_in[1][1][1] + spinor_in[3][1][0];


    proj_result[2][0][0] = spinor_in[0][2][0] + spinor_in[2][2][1];
    proj_result[2][0][1] = spinor_in[0][2][1] - spinor_in[2][2][0];
    proj_result[2][1][0] = spinor_in[1][2][0] - spinor_in[3][2][1];
    proj_result[2][1][1] = spinor_in[1][2][1] + spinor_in[3][2][0];

    su3_adj_mult(proj_matvec, u, proj_result);


    // Top half
    // rearrange for linear storage order 
    spinor_out[0][0][0] = upper_sum[0][0][0] + proj_matvec[0][0][0];
    spinor_out[0][0][1] = upper_sum[0][0][1] + proj_matvec[0][0][1];
    spinor_out[0][1][0] = upper_sum[1][0][0] + proj_matvec[1][0][0];
    spinor_out[0][1][1] = upper_sum[1][0][1] + proj_matvec[1][0][1];
    spinor_out[0][2][0] = upper_sum[2][0][0] + proj_matvec[2][0][0];
    spinor_out[0][2][1] = upper_sum[2][0][1] + proj_matvec[2][0][1];

    spinor_out[1][0][0] = upper_sum[0][1][0] + proj_matvec[0][1][0];
    spinor_out[1][0][1] = upper_sum[0][1][1] + proj_matvec[0][1][1];
    spinor_out[1][1][0] = upper_sum[1][1][0] + proj_matvec[1][1][0];
    spinor_out[1][1][1] = upper_sum[1][1][1] + proj_matvec[1][1][1];
    spinor_out[1][2][0] = upper_sum[2][1][0] + proj_matvec[2][1][0];
    spinor_out[1][2][1] = upper_sum[2][1][1] + proj_matvec[2][1][1];

    /*
     * The bottom components of be may be reconstructed using the formula
     *      ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
     *      ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
     */
    

    spinor_out[2][0][0] = lower_sum[0][0][0] - proj_matvec[0][0][1];
    spinor_out[2][0][1] = lower_sum[0][0][1] + proj_matvec[0][0][0];
    spinor_out[2][1][0] = lower_sum[1][0][0] - proj_matvec[1][0][1];
    spinor_out[2][1][1] = lower_sum[1][0][1] + proj_matvec[1][0][0];
    spinor_out[2][2][0] = lower_sum[2][0][0] - proj_matvec[2][0][1];
    spinor_out[2][2][1] = lower_sum[2][0][1] + proj_matvec[2][0][0];

    spinor_out[3][0][0] = lower_sum[0][1][0] + proj_matvec[0][1][1];
    spinor_out[3][0][1] = lower_sum[0][1][1] - proj_matvec[0][1][0];
    spinor_out[3][1][0] = lower_sum[1][1][0] + proj_matvec[1][1][1];
    spinor_out[3][1][1] = lower_sum[1][1][1] - proj_matvec[1][1][0];
    spinor_out[3][2][0] = lower_sum[2][1][0] + proj_matvec[2][1][1];
    spinor_out[3][2][1] = lower_sum[2][1][1] - proj_matvec[2][1][0];

  }

  inline
  void dslash_minus_dir3_forward_add( FourSpinor  spinor_in,
				      GaugeMatrix  u,
				      HalfSpinor  upper_sum,
				      HalfSpinor  lower_sum)
  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;
  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
     *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
     *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
     *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r + a2r} + i{a0i + a2i} )
     *      ( b1r + i b1i )     ( {a1r + a3r} + i{a1i + a3i} )
     */

    proj_result[0][0][0] = spinor_in[0][0][0] + spinor_in[2][0][0];
    proj_result[0][0][1] = spinor_in[0][0][1] + spinor_in[2][0][1];
    proj_result[0][1][0] = spinor_in[1][0][0] + spinor_in[3][0][0];
    proj_result[0][1][1] = spinor_in[1][0][1] + spinor_in[3][0][1];

    proj_result[1][0][0] = spinor_in[0][1][0] + spinor_in[2][1][0];
    proj_result[1][0][1] = spinor_in[0][1][1] + spinor_in[2][1][1];
    proj_result[1][1][0] = spinor_in[1][1][0] + spinor_in[3][1][0];
    proj_result[1][1][1] = spinor_in[1][1][1] + spinor_in[3][1][1];

    proj_result[2][0][0] = spinor_in[0][2][0] + spinor_in[2][2][0];
    proj_result[2][0][1] = spinor_in[0][2][1] + spinor_in[2][2][1];
    proj_result[2][1][0] = spinor_in[1][2][0] + spinor_in[3][2][0];
    proj_result[2][1][1] = spinor_in[1][2][1] + spinor_in[3][2][1];


    su3_mult(proj_matvec, u, proj_result);

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

    /*
     *      ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
     *      ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
     */
     lower_sum[0][0][0]+= proj_matvec[0][0][0];
     lower_sum[0][0][1]+= proj_matvec[0][0][1];
     lower_sum[0][1][0]+= proj_matvec[0][1][0];
     lower_sum[0][1][1]+= proj_matvec[0][1][1];

     lower_sum[1][0][0]+= proj_matvec[1][0][0];
     lower_sum[1][0][1]+= proj_matvec[1][0][1];
     lower_sum[1][1][0]+= proj_matvec[1][1][0];
     lower_sum[1][1][1]+= proj_matvec[1][1][1];
 
     lower_sum[2][0][0]+= proj_matvec[2][0][0];
     lower_sum[2][0][1]+= proj_matvec[2][0][1];
     lower_sum[2][1][0]+= proj_matvec[2][1][0];
     lower_sum[2][1][1]+= proj_matvec[2][1][1];

  }

  inline
  void dslash_minus_dir3_backward_add_store( FourSpinor  spinor_in,
					     GaugeMatrix  u,
					     HalfSpinor  upper_sum,
					     HalfSpinor  lower_sum,
					     FourSpinor spinor_out)
       
  {
    HalfSpinor proj_result;
    HalfSpinor proj_matvec;
    /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
     *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
     *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r - a2r} + i{a0i - a2i} )
     *      ( b1r + i b1i )     ( {a1r - a3r} + i{a1i - a3i} )
     */

    proj_result[0][0][0] = spinor_in[0][0][0] - spinor_in[2][0][0];
    proj_result[0][0][1] = spinor_in[0][0][1] - spinor_in[2][0][1];
    proj_result[0][1][0] = spinor_in[1][0][0] - spinor_in[3][0][0];
    proj_result[0][1][1] = spinor_in[1][0][1] - spinor_in[3][0][1];

    proj_result[1][0][0] = spinor_in[0][1][0] - spinor_in[2][1][0];
    proj_result[1][0][1] = spinor_in[0][1][1] - spinor_in[2][1][1];
    proj_result[1][1][0] = spinor_in[1][1][0] - spinor_in[3][1][0];
    proj_result[1][1][1] = spinor_in[1][1][1] - spinor_in[3][1][1];

    proj_result[2][0][0] = spinor_in[0][2][0] - spinor_in[2][2][0];
    proj_result[2][0][1] = spinor_in[0][2][1] - spinor_in[2][2][1];
    proj_result[2][1][0] = spinor_in[1][2][0] - spinor_in[3][2][0];
    proj_result[2][1][1] = spinor_in[1][2][1] - spinor_in[3][2][1];


    su3_adj_mult(proj_matvec, u, proj_result);
    
    // Top half
    // rearrange for linear storage order 
    spinor_out[0][0][0] = upper_sum[0][0][0] + proj_matvec[0][0][0];
    spinor_out[0][0][1] = upper_sum[0][0][1] + proj_matvec[0][0][1];
    spinor_out[0][1][0] = upper_sum[1][0][0] + proj_matvec[1][0][0];
    spinor_out[0][1][1] = upper_sum[1][0][1] + proj_matvec[1][0][1];
    spinor_out[0][2][0] = upper_sum[2][0][0] + proj_matvec[2][0][0];
    spinor_out[0][2][1] = upper_sum[2][0][1] + proj_matvec[2][0][1];

    spinor_out[1][0][0] = upper_sum[0][1][0] + proj_matvec[0][1][0];
    spinor_out[1][0][1] = upper_sum[0][1][1] + proj_matvec[0][1][1];
    spinor_out[1][1][0] = upper_sum[1][1][0] + proj_matvec[1][1][0];
    spinor_out[1][1][1] = upper_sum[1][1][1] + proj_matvec[1][1][1];
    spinor_out[1][2][0] = upper_sum[2][1][0] + proj_matvec[2][1][0];
    spinor_out[1][2][1] = upper_sum[2][1][1] + proj_matvec[2][1][1];

      /*
     *      ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
     *      ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
     */


    spinor_out[2][0][0] = lower_sum[0][0][0] - proj_matvec[0][0][0];
    spinor_out[2][0][1] = lower_sum[0][0][1] - proj_matvec[0][0][1];
    spinor_out[2][1][0] = lower_sum[1][0][0] - proj_matvec[1][0][0];
    spinor_out[2][1][1] = lower_sum[1][0][1] - proj_matvec[1][0][1];
    spinor_out[2][2][0] = lower_sum[2][0][0] - proj_matvec[2][0][0];
    spinor_out[2][2][1] = lower_sum[2][0][1] - proj_matvec[2][0][1];

    spinor_out[3][0][0] = lower_sum[0][1][0] - proj_matvec[0][1][0];
    spinor_out[3][0][1] = lower_sum[0][1][1] - proj_matvec[0][1][1];
    spinor_out[3][1][0] = lower_sum[1][1][0] - proj_matvec[1][1][0];
    spinor_out[3][1][1] = lower_sum[1][1][1] - proj_matvec[1][1][1];
    spinor_out[3][2][0] = lower_sum[2][1][0] - proj_matvec[2][1][0];
    spinor_out[3][2][1] = lower_sum[2][1][1] - proj_matvec[2][1][1];

  
  }

  } // Namespace 32 bit
   
} // Namespace 
#endif

    



