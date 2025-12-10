#ifndef CPP_DSLASH_PARSCALAR_DECOMP_64BIT_C_H
#define CPP_DSLASH_PARSCALAR_DECOMP_64BIT_C_H

#warning "Using CStuff"
#include <cpp_dslash_types.h>

using namespace CPlusPlusWilsonDslash::Dslash64BitTypes;

namespace CPlusPlusWilsonDslash { 
  namespace  DslashParscalar64Bit { 

inline void decomp_gamma0_minus( FourSpinor src, HalfSpinor dst) 
{
    /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
     *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
     *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r + a3i} + i{a0i - a3r} )
     *      ( b1r + i b1i )     ( {a1r + a2i} + i{a1i - a2r} )
     */

     dst[0][0][0] = src[0][0][0] + src[3][0][1];
     dst[0][0][1] = src[0][0][1] - src[3][0][0];
     dst[1][0][0] = src[1][0][0] + src[2][0][1];
     dst[1][0][1] = src[1][0][1] - src[2][0][0];

     dst[0][1][0] = src[0][1][0] + src[3][1][1];
     dst[0][1][1] = src[0][1][1] - src[3][1][0];
     dst[1][1][0] = src[1][1][0] + src[2][1][1];
     dst[1][1][1] = src[1][1][1] - src[2][1][0];

     dst[0][2][0] = src[0][2][0] + src[3][2][1];
     dst[0][2][1] = src[0][2][1] - src[3][2][0];
     dst[1][2][0] = src[1][2][0] + src[2][2][1];
     dst[1][2][1] = src[1][2][1] - src[2][2][0];
}


inline void decomp_gamma1_minus( FourSpinor src, HalfSpinor dst)
{
  /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
   *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
   *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
	 
   * Therefore the top components are
      
   *      ( b0r + i b0i )  =  ( {a0r + a3r} + i{a0i + a3i} )
   *      ( b1r + i b1i )     ( {a1r - a2r} + i{a1i - a2i} )
   */
    dst[0][0][0] = src[0][0][0] + src[3][0][0];
    dst[0][0][1] = src[0][0][1] + src[3][0][1];
    dst[1][0][0] = src[1][0][0] - src[2][0][0];
    dst[1][0][1] = src[1][0][1] - src[2][0][1];

    dst[0][1][0] = src[0][1][0] + src[3][1][0];
    dst[0][1][1] = src[0][1][1] + src[3][1][1];
    dst[1][1][0] = src[1][1][0] - src[2][1][0];
    dst[1][1][1] = src[1][1][1] - src[2][1][1];

    dst[0][2][0] = src[0][2][0] + src[3][2][0];
    dst[0][2][1] = src[0][2][1] + src[3][2][1];
    dst[1][2][0] = src[1][2][0] - src[2][2][0];
    dst[1][2][1] = src[1][2][1] - src[2][2][1];
     
}

inline void decomp_gamma2_minus( FourSpinor src, HalfSpinor dst) 
{

    /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
     *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
     *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
     *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
     
     */
    dst[0][0][0] = src[0][0][0] + src[2][0][1];
    dst[0][0][1] = src[0][0][1] - src[2][0][0];
    dst[1][0][0] = src[1][0][0] - src[3][0][1];
    dst[1][0][1] = src[1][0][1] + src[3][0][0];

    dst[0][1][0] = src[0][1][0] + src[2][1][1];
    dst[0][1][1] = src[0][1][1] - src[2][1][0];
    dst[1][1][0] = src[1][1][0] - src[3][1][1];
    dst[1][1][1] = src[1][1][1] + src[3][1][0];

    dst[0][2][0] = src[0][2][0] + src[2][2][1];
    dst[0][2][1] = src[0][2][1] - src[2][2][0];
    dst[1][2][0] = src[1][2][0] - src[3][2][1];
    dst[1][2][1] = src[1][2][1] + src[3][2][0];

}

inline void decomp_gamma3_minus( FourSpinor src, HalfSpinor dst) 

{
    /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
     *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
     *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
     *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
     
     * Therefore the top components are
     
     *      ( b0r + i b0i )  =  ( {a0r - a2r} + i{a0i - a2i} )
     *      ( b1r + i b1i )     ( {a1r - a3r} + i{a1i - a3i} )
     */

    dst[0][0][0] = src[0][0][0] - src[2][0][0];
    dst[0][0][1] = src[0][0][1] - src[2][0][1];
    dst[1][0][0] = src[1][0][0] - src[3][0][0];
    dst[1][0][1] = src[1][0][1] - src[3][0][1];

    dst[0][1][0] = src[0][1][0] - src[2][1][0];
    dst[0][1][1] = src[0][1][1] - src[2][1][1];
    dst[1][1][0] = src[1][1][0] - src[3][1][0];
    dst[1][1][1] = src[1][1][1] - src[3][1][1];

    dst[0][2][0] = src[0][2][0] - src[2][2][0];
    dst[0][2][1] = src[0][2][1] - src[2][2][1];
    dst[1][2][0] = src[1][2][0] - src[3][2][0];
    dst[1][2][1] = src[1][2][1] - src[3][2][1];

}



inline void decomp_gamma0_plus( FourSpinor src, HalfSpinor dst) 
{
  /*                              ( 1  0  0 +i)  ( a0 )    ( a0 + i a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1 +i  0)  ( a1 )  = ( a1 + i a2 )
   *                    0         ( 0 -i  1  0)  ( a2 )    ( a2 - i a1 )
   *                              (-i  0  0  1)  ( a3 )    ( a3 - i a0 )
   
   * Therefore the top components are
   
   *      ( b0r + i b0i )  =  ( {a0r - a3i} + i{a0i + a3r} )
   *      ( b1r + i b1i )     ( {a1r - a2i} + i{a1i + a2r} )
   */
  
  dst[0][0][0] = src[0][0][0] - src[3][0][1];
  dst[0][0][1] = src[0][0][1] + src[3][0][0];
  dst[1][0][0] = src[1][0][0] - src[2][0][1];
  dst[1][0][1] = src[1][0][1] + src[2][0][0];

  dst[0][1][0] = src[0][1][0] - src[3][1][1];
  dst[0][1][1] = src[0][1][1] + src[3][1][0];
  dst[1][1][0] = src[1][1][0] - src[2][1][1];
  dst[1][1][1] = src[1][1][1] + src[2][1][0];
  
  dst[0][2][0] = src[0][2][0] - src[3][2][1];
  dst[0][2][1] = src[0][2][1] + src[3][2][0];
  dst[1][2][0] = src[1][2][0] - src[2][2][1];
  dst[1][2][1] = src[1][2][1] + src[2][2][0];
  

}

inline void decomp_gamma1_plus( FourSpinor src, HalfSpinor dst) 
{

  
  /*                              ( 1  0  0 -1)  ( a0 )    ( a0 - a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  1  0)  ( a1 )  = ( a1 + a2 )
   *                    1         ( 0  1  1  0)  ( a2 )    ( a2 + a1 )
   *                              (-1  0  0  1)  ( a3 )    ( a3 - a0 )

   * Therefore the top components are
   
   *      ( b0r + i b0i )  =  ( {a0r - a3r} + i{a0i - a3i} )
   *      ( b1r + i b1i )     ( {a1r + a2r} + i{a1i + a2i} )
   */
  
  dst[0][0][0] = src[0][0][0] - src[3][0][0];
  dst[0][0][1] = src[0][0][1] - src[3][0][1];
  dst[1][0][0] = src[1][0][0] + src[2][0][0];
  dst[1][0][1] = src[1][0][1] + src[2][0][1];
  
  dst[0][1][0] = src[0][1][0] - src[3][1][0];
  dst[0][1][1] = src[0][1][1] - src[3][1][1];
  dst[1][1][0] = src[1][1][0] + src[2][1][0];
  dst[1][1][1] = src[1][1][1] + src[2][1][1];
  
  dst[0][2][0] = src[0][2][0] - src[3][2][0];
  dst[0][2][1] = src[0][2][1] - src[3][2][1];
  dst[1][2][0] = src[1][2][0] + src[2][2][0];
  dst[1][2][1] = src[1][2][1] + src[2][2][1];
  
}

inline void decomp_gamma2_plus( FourSpinor src, HalfSpinor dst) 
{
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
  
  dst[0][0][0] = src[0][0][0] - src[2][0][1];
  dst[0][0][1] = src[0][0][1] + src[2][0][0];
  dst[1][0][0] = src[1][0][0] + src[3][0][1];
  dst[1][0][1] = src[1][0][1] - src[3][0][0];
  
  dst[0][1][0] = src[0][1][0] - src[2][1][1];
  dst[0][1][1] = src[0][1][1] + src[2][1][0];
  dst[1][1][0] = src[1][1][0] + src[3][1][1];
  dst[1][1][1] = src[1][1][1] - src[3][1][0];
  
  dst[0][2][0] = src[0][2][0] - src[2][2][1];
  dst[0][2][1] = src[0][2][1] + src[2][2][0];
  dst[1][2][0] = src[1][2][0] + src[3][2][1];
  dst[1][2][1] = src[1][2][1] - src[3][2][0];
  
}

inline void decomp_gamma3_plus( FourSpinor src, HalfSpinor dst) 
{
  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
   *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
   *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
   
   * Therefore the top components are
   
   *      ( b0r + i b0i )  =  ( {a0r + a2r} + i{a0i + a2i} )
   *      ( b1r + i b1i )     ( {a1r + a3r} + i{a1i + a3i} )
   */
  
  dst[0][0][0] = src[0][0][0] + src[2][0][0];
  dst[0][0][1] = src[0][0][1] + src[2][0][1];
  dst[1][0][0] = src[1][0][0] + src[3][0][0];
  dst[1][0][1] = src[1][0][1] + src[3][0][1];
  
  dst[0][1][0] = src[0][1][0] + src[2][1][0];
  dst[0][1][1] = src[0][1][1] + src[2][1][1];
  dst[1][1][0] = src[1][1][0] + src[3][1][0];
  dst[1][1][1] = src[1][1][1] + src[3][1][1];
  
  dst[0][2][0] = src[0][2][0] + src[2][2][0];
  dst[0][2][1] = src[0][2][1] + src[2][2][1];
  dst[1][2][0] = src[1][2][0] + src[3][2][0];
  dst[1][2][1] = src[1][2][1] + src[3][2][1];
  

}

  }
}

#endif
