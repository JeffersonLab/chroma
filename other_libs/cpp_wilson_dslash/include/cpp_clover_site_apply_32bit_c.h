#ifndef CPP_CLOVER_SITE_APPLY_32_BIT_C_H
#define CPP_CLOVER_SITE_APPLY_32_BIT_C_H

#include <cpp_clover_types.h>
#include <cpp_dslash_types.h>


namespace CPlusPlusClover { 

  // Use 32bit types
  using namespace Clover32BitTypes;
using namespace Dslash32BitTypes;

  namespace CPlusPlusClover32Bit { 


    inline 
    void cloverSiteApply(FourSpinor result, const CloverTerm clov, const FourSpinor src) 
    {

      // Unrolled version 3. 
      // Took unrolled version 2 and wrote out in real() and imag() 
      // parts. Rearranged so that all the reals follow each other
      //  in the output so that we can write linearly


      result[0][0][0]  = clov[0].diag[0]  * src[0][0][0];
      result[0][0][0] += clov[0].off_diag[0][0]  * src[0][1][0];
      result[0][0][0] += clov[0].off_diag[0][1]  * src[0][1][1];
      result[0][0][0] += clov[0].off_diag[1][0]  * src[0][2][0];
      result[0][0][0] += clov[0].off_diag[1][1]  * src[0][2][1];
      result[0][0][0] += clov[0].off_diag[3][0]  * src[1][0][0];
      result[0][0][0] += clov[0].off_diag[3][1]  * src[1][0][1];
      result[0][0][0] += clov[0].off_diag[6][0]  * src[1][1][0];
      result[0][0][0] += clov[0].off_diag[6][1]  * src[1][1][1];
      result[0][0][0] += clov[0].off_diag[10][0] * src[1][2][0];
      result[0][0][0] += clov[0].off_diag[10][1] * src[1][2][1];

      result[0][0][1]  = clov[0].diag[0]  * src[0][0][1];
      result[0][0][1] += clov[0].off_diag[0][0]  * src[0][1][1];
      result[0][0][1] -= clov[0].off_diag[0][1]  * src[0][1][0];
      result[0][0][1] += clov[0].off_diag[3][0]  * src[1][0][1];
      result[0][0][1] -= clov[0].off_diag[3][1]  * src[1][0][0];
      result[0][0][1] += clov[0].off_diag[1][0]  * src[0][2][1];
      result[0][0][1] -= clov[0].off_diag[1][1]  * src[0][2][0];
      result[0][0][1] += clov[0].off_diag[6][0]  * src[1][1][1];
      result[0][0][1] -= clov[0].off_diag[6][1]  * src[1][1][0];
      result[0][0][1] += clov[0].off_diag[10][0] * src[1][2][1];
      result[0][0][1] -= clov[0].off_diag[10][1] * src[1][2][0];


      result[0][1][0]  = clov[0].diag[ 1] * src[0][1][0];
      result[0][1][0] += clov[0].off_diag[ 0][0] * src[0][0][0];
      result[0][1][0] -= clov[0].off_diag[ 0][1] * src[0][0][1];
      result[0][1][0] += clov[0].off_diag[ 2][0] * src[0][2][0];
      result[0][1][0] += clov[0].off_diag[ 2][1] * src[0][2][1];
      result[0][1][0] += clov[0].off_diag[ 4][0] * src[1][0][0];
      result[0][1][0] += clov[0].off_diag[ 4][1] * src[1][0][1];
      result[0][1][0] += clov[0].off_diag[ 7][0] * src[1][1][0];
      result[0][1][0] += clov[0].off_diag[ 7][1] * src[1][1][1];
      result[0][1][0] += clov[0].off_diag[11][0] * src[1][2][0];
      result[0][1][0] += clov[0].off_diag[11][1] * src[1][2][1];


      result[0][1][1]  = clov[0].diag[ 1] * src[0][1][1];
      result[0][1][1] += clov[0].off_diag[ 0][0] * src[0][0][1];
      result[0][1][1] += clov[0].off_diag[ 0][1] * src[0][0][0];
      result[0][1][1] += clov[0].off_diag[ 2][0] * src[0][2][1];
      result[0][1][1] -= clov[0].off_diag[ 2][1] * src[0][2][0];
      result[0][1][1] += clov[0].off_diag[ 4][0] * src[1][0][1];
      result[0][1][1] -= clov[0].off_diag[ 4][1] * src[1][0][0];
      result[0][1][1] += clov[0].off_diag[ 7][0] * src[1][1][1];
      result[0][1][1] -= clov[0].off_diag[ 7][1] * src[1][1][0];
      result[0][1][1] += clov[0].off_diag[11][0] * src[1][2][1];
      result[0][1][1] -= clov[0].off_diag[11][1] * src[1][2][0];


      result[0][2][0] = clov[0].diag[ 2]  * src[0][2][0];
      result[0][2][0] += clov[0].off_diag[ 1][0] * src[0][0][0];
      result[0][2][0] -= clov[0].off_diag[ 1][1] * src[0][0][1];
      result[0][2][0] += clov[0].off_diag[ 2][0] * src[0][1][0];
      result[0][2][0] -= clov[0].off_diag[ 2][1] * src[0][1][1];
      result[0][2][0] += clov[0].off_diag[5][0]  * src[1][0][0];
      result[0][2][0] += clov[0].off_diag[5][1]  * src[1][0][1];
      result[0][2][0] += clov[0].off_diag[8][0]  * src[1][1][0];
      result[0][2][0] += clov[0].off_diag[8][1]  * src[1][1][1];
      result[0][2][0] += clov[0].off_diag[12][0] * src[1][2][0];
      result[0][2][0] += clov[0].off_diag[12][1] * src[1][2][1];


      result[0][2][1] = clov[0].diag[ 2]  * src[0][2][1];
      result[0][2][1] += clov[0].off_diag[ 1][0] * src[0][0][1];
      result[0][2][1] += clov[0].off_diag[ 1][1] * src[0][0][0];
      result[0][2][1] += clov[0].off_diag[ 2][0] * src[0][1][1];
      result[0][2][1] += clov[0].off_diag[ 2][1] * src[0][1][0];
      result[0][2][1] += clov[0].off_diag[5][0]  * src[1][0][1];
      result[0][2][1] -= clov[0].off_diag[5][1]  * src[1][0][0];
      result[0][2][1] += clov[0].off_diag[8][0]  * src[1][1][1];
      result[0][2][1] -= clov[0].off_diag[8][1]  * src[1][1][0];
      result[0][2][1] += clov[0].off_diag[12][0] * src[1][2][1];
      result[0][2][1] -= clov[0].off_diag[12][1] * src[1][2][0];


      result[1][0][0]  = clov[0].diag[ 3] * src[1][0][0];
      result[1][0][0] += clov[0].off_diag[ 3][0] * src[0][0][0];
      result[1][0][0] -= clov[0].off_diag[ 3][1] * src[0][0][1];
      result[1][0][0] += clov[0].off_diag[ 4][0] * src[0][1][0];
      result[1][0][0] -= clov[0].off_diag[ 4][1] * src[0][1][1];
      result[1][0][0] += clov[0].off_diag[ 5][0] * src[0][2][0];
      result[1][0][0] -= clov[0].off_diag[ 5][1] * src[0][2][1];
      result[1][0][0] += clov[0].off_diag[ 9][0] * src[1][1][0];
      result[1][0][0] += clov[0].off_diag[ 9][1] * src[1][1][1];
      result[1][0][0] += clov[0].off_diag[13][0] * src[1][2][0];
      result[1][0][0] += clov[0].off_diag[13][1] * src[1][2][1];


      result[1][0][1]  = clov[0].diag[ 3] * src[1][0][1];
      result[1][0][1] += clov[0].off_diag[ 3][0] * src[0][0][1];
      result[1][0][1] += clov[0].off_diag[ 3][1] * src[0][0][0];
      result[1][0][1] += clov[0].off_diag[ 4][0] * src[0][1][1];
      result[1][0][1] += clov[0].off_diag[ 4][1] * src[0][1][0];
      result[1][0][1] += clov[0].off_diag[ 5][0] * src[0][2][1];
      result[1][0][1] += clov[0].off_diag[ 5][1] * src[0][2][0];
      result[1][0][1] += clov[0].off_diag[ 9][0] * src[1][1][1];
      result[1][0][1] -= clov[0].off_diag[ 9][1] * src[1][1][0];
      result[1][0][1] += clov[0].off_diag[13][0] * src[1][2][1];
      result[1][0][1] -= clov[0].off_diag[13][1] * src[1][2][0];


      result[1][1][0]  = clov[0].diag[ 4] * src[1][1][0];
      result[1][1][0] += clov[0].off_diag[ 6][0] * src[0][0][0];
      result[1][1][0] -= clov[0].off_diag[ 6][1] * src[0][0][1];
      result[1][1][0] += clov[0].off_diag[ 7][0] * src[0][1][0];
      result[1][1][0] -= clov[0].off_diag[ 7][1] * src[0][1][1];
      result[1][1][0] += clov[0].off_diag[ 8][0] * src[0][2][0];
      result[1][1][0] -= clov[0].off_diag[ 8][1] * src[0][2][1];
      result[1][1][0] += clov[0].off_diag[ 9][0] * src[1][0][0];
      result[1][1][0] -= clov[0].off_diag[ 9][1] * src[1][0][1];
      result[1][1][0] += clov[0].off_diag[14][0] * src[1][2][0];
      result[1][1][0] += clov[0].off_diag[14][1] * src[1][2][1];


      result[1][1][1]  = clov[0].diag[ 4] * src[1][1][1];
      result[1][1][1] += clov[0].off_diag[ 6][0] * src[0][0][1];
      result[1][1][1] += clov[0].off_diag[ 6][1] * src[0][0][0];
      result[1][1][1] += clov[0].off_diag[ 7][0] * src[0][1][1];
      result[1][1][1] += clov[0].off_diag[ 7][1] * src[0][1][0];
      result[1][1][1] += clov[0].off_diag[ 8][0] * src[0][2][1];
      result[1][1][1] += clov[0].off_diag[ 8][1] * src[0][2][0];
      result[1][1][1] += clov[0].off_diag[ 9][0] * src[1][0][1];
      result[1][1][1] += clov[0].off_diag[ 9][1] * src[1][0][0];
      result[1][1][1] += clov[0].off_diag[14][0] * src[1][2][1];
      result[1][1][1] -= clov[0].off_diag[14][1] * src[1][2][0];


      result[1][2][0]  = clov[0].diag[ 5] * src[1][2][0];
      result[1][2][0] += clov[0].off_diag[10][0] * src[0][0][0];
      result[1][2][0] -= clov[0].off_diag[10][1] * src[0][0][1];
      result[1][2][0] += clov[0].off_diag[11][0] * src[0][1][0];
      result[1][2][0] -= clov[0].off_diag[11][1] * src[0][1][1];
      result[1][2][0] += clov[0].off_diag[12][0] * src[0][2][0];
      result[1][2][0] -= clov[0].off_diag[12][1] * src[0][2][1];
      result[1][2][0] += clov[0].off_diag[13][0] * src[1][0][0];
      result[1][2][0] -= clov[0].off_diag[13][1] * src[1][0][1];
      result[1][2][0] += clov[0].off_diag[14][0] * src[1][1][0];
      result[1][2][0] -= clov[0].off_diag[14][1] * src[1][1][1];


      result[1][2][1]  = clov[0].diag[ 5] * src[1][2][1];
      result[1][2][1] += clov[0].off_diag[10][0] * src[0][0][1];
      result[1][2][1] += clov[0].off_diag[10][1] * src[0][0][0];
      result[1][2][1] += clov[0].off_diag[11][0] * src[0][1][1];
      result[1][2][1] += clov[0].off_diag[11][1] * src[0][1][0];
      result[1][2][1] += clov[0].off_diag[12][0] * src[0][2][1];
      result[1][2][1] += clov[0].off_diag[12][1] * src[0][2][0];
      result[1][2][1] += clov[0].off_diag[13][0] * src[1][0][1];
      result[1][2][1] += clov[0].off_diag[13][1] * src[1][0][0];
      result[1][2][1] += clov[0].off_diag[14][0] * src[1][1][1];
      result[1][2][1] += clov[0].off_diag[14][1] * src[1][1][0];


      result[2][0][0]  = clov[1].diag[0]  * src[2][0][0];
      result[2][0][0] += clov[1].off_diag[0][0]  * src[2][1][0];
      result[2][0][0] += clov[1].off_diag[0][1]  * src[2][1][1];
      result[2][0][0] += clov[1].off_diag[1][0]  * src[2][2][0];
      result[2][0][0] += clov[1].off_diag[1][1]  * src[2][2][1];
      result[2][0][0] += clov[1].off_diag[3][0]  * src[3][0][0];
      result[2][0][0] += clov[1].off_diag[3][1]  * src[3][0][1];
      result[2][0][0] += clov[1].off_diag[6][0]  * src[3][1][0];
      result[2][0][0] += clov[1].off_diag[6][1]  * src[3][1][1];
      result[2][0][0] += clov[1].off_diag[10][0] * src[3][2][0];
      result[2][0][0] += clov[1].off_diag[10][1] * src[3][2][1];

      result[2][0][1]  = clov[1].diag[0]  * src[2][0][1];
      result[2][0][1] += clov[1].off_diag[0][0]  * src[2][1][1];
      result[2][0][1] -= clov[1].off_diag[0][1]  * src[2][1][0];
      result[2][0][1] += clov[1].off_diag[1][0]  * src[2][2][1];
      result[2][0][1] -= clov[1].off_diag[1][1]  * src[2][2][0];
      result[2][0][1] += clov[1].off_diag[3][0]  * src[3][0][1];
      result[2][0][1] -= clov[1].off_diag[3][1]  * src[3][0][0];
      result[2][0][1] += clov[1].off_diag[6][0]  * src[3][1][1];
      result[2][0][1] -= clov[1].off_diag[6][1]  * src[3][1][0];
      result[2][0][1] += clov[1].off_diag[10][0] * src[3][2][1];
      result[2][0][1] -= clov[1].off_diag[10][1] * src[3][2][0];


      result[2][1][0]  = clov[1].diag[ 1]  * src[2][1][0];
      result[2][1][0] += clov[1].off_diag[ 0][0]  * src[2][0][0];
      result[2][1][0] -= clov[1].off_diag[ 0][1]  * src[2][0][1];
      result[2][1][0] += clov[1].off_diag[ 2][0]  * src[2][2][0];
      result[2][1][0] += clov[1].off_diag[ 2][1]  * src[2][2][1];
      result[2][1][0] += clov[1].off_diag[ 4][0]  * src[3][0][0];
      result[2][1][0] += clov[1].off_diag[ 4][1]  * src[3][0][1];
      result[2][1][0] += clov[1].off_diag[ 7][0]  * src[3][1][0];
      result[2][1][0] += clov[1].off_diag[ 7][1]  * src[3][1][1];
      result[2][1][0] += clov[1].off_diag[11][0]  * src[3][2][0];
      result[2][1][0] += clov[1].off_diag[11][1]  * src[3][2][1];

      result[2][1][1]  = clov[1].diag[ 1]  * src[2][1][1];
      result[2][1][1] += clov[1].off_diag[ 0][0]  * src[2][0][1];
      result[2][1][1] += clov[1].off_diag[ 0][1]  * src[2][0][0];
      result[2][1][1] += clov[1].off_diag[ 2][0]  * src[2][2][1];
      result[2][1][1] -= clov[1].off_diag[ 2][1]  * src[2][2][0];
      result[2][1][1] += clov[1].off_diag[ 4][0]  * src[3][0][1];
      result[2][1][1] -= clov[1].off_diag[ 4][1]  * src[3][0][0];
      result[2][1][1] += clov[1].off_diag[ 7][0]  * src[3][1][1];
      result[2][1][1] -= clov[1].off_diag[ 7][1]  * src[3][1][0];
      result[2][1][1] += clov[1].off_diag[11][0]  * src[3][2][1];
      result[2][1][1] -= clov[1].off_diag[11][1]  * src[3][2][0];


      result[2][2][0]  = clov[1].diag[ 2]  * src[2][2][0];
      result[2][2][0] += clov[1].off_diag[ 1][0]  * src[2][0][0];
      result[2][2][0] -= clov[1].off_diag[ 1][1]  * src[2][0][1];
      result[2][2][0] += clov[1].off_diag[ 2][0]  * src[2][1][0];
      result[2][2][0] -= clov[1].off_diag[ 2][1]  * src[2][1][1];
      result[2][2][0] += clov[1].off_diag[5][0]   * src[3][0][0];
      result[2][2][0] += clov[1].off_diag[5][1]   * src[3][0][1];
      result[2][2][0] += clov[1].off_diag[8][0]   * src[3][1][0];
      result[2][2][0] += clov[1].off_diag[8][1]   * src[3][1][1];
      result[2][2][0] += clov[1].off_diag[12][0]  * src[3][2][0];
      result[2][2][0] += clov[1].off_diag[12][1]  * src[3][2][1];

      result[2][2][1] = clov[1].diag[ 2]   * src[2][2][1];
      result[2][2][1] += clov[1].off_diag[ 1][0]  * src[2][0][1];
      result[2][2][1] += clov[1].off_diag[ 1][1]  * src[2][0][0];
      result[2][2][1] += clov[1].off_diag[ 2][0]  * src[2][1][1];
      result[2][2][1] += clov[1].off_diag[ 2][1]  * src[2][1][0];
      result[2][2][1] += clov[1].off_diag[5][0]   * src[3][0][1];
      result[2][2][1] -= clov[1].off_diag[5][1]   * src[3][0][0];
      result[2][2][1] += clov[1].off_diag[8][0]   * src[3][1][1];
      result[2][2][1] -= clov[1].off_diag[8][1]   * src[3][1][0];
      result[2][2][1] += clov[1].off_diag[12][0]  * src[3][2][1];
      result[2][2][1] -= clov[1].off_diag[12][1]  * src[3][2][0];


      result[3][0][0]  = clov[1].diag[ 3]  * src[3][0][0];
      result[3][0][0] += clov[1].off_diag[ 3][0]  * src[2][0][0];
      result[3][0][0] -= clov[1].off_diag[ 3][1]  * src[2][0][1];
      result[3][0][0] += clov[1].off_diag[ 4][0]  * src[2][1][0];
      result[3][0][0] -= clov[1].off_diag[ 4][1]  * src[2][1][1];
      result[3][0][0] += clov[1].off_diag[ 5][0]  * src[2][2][0];
      result[3][0][0] -= clov[1].off_diag[ 5][1]  * src[2][2][1];
      result[3][0][0] += clov[1].off_diag[ 9][0]  * src[3][1][0];
      result[3][0][0] += clov[1].off_diag[ 9][1]  * src[3][1][1];
      result[3][0][0] += clov[1].off_diag[13][0]  * src[3][2][0];
      result[3][0][0] += clov[1].off_diag[13][1]  * src[3][2][1];

      result[3][0][1]  = clov[1].diag[ 3]  * src[3][0][1];
      result[3][0][1] += clov[1].off_diag[ 3][0]  * src[2][0][1];
      result[3][0][1] += clov[1].off_diag[ 3][1]  * src[2][0][0];
      result[3][0][1] += clov[1].off_diag[ 4][0]  * src[2][1][1];
      result[3][0][1] += clov[1].off_diag[ 4][1]  * src[2][1][0];
      result[3][0][1] += clov[1].off_diag[ 5][0]  * src[2][2][1];
      result[3][0][1] += clov[1].off_diag[ 5][1]  * src[2][2][0];
      result[3][0][1] += clov[1].off_diag[ 9][0]  * src[3][1][1];
      result[3][0][1] -= clov[1].off_diag[ 9][1]  * src[3][1][0];
      result[3][0][1] += clov[1].off_diag[13][0]  * src[3][2][1];
      result[3][0][1] -= clov[1].off_diag[13][1]  * src[3][2][0];


      result[3][1][0]  = clov[1].diag[ 4]  * src[3][1][0];
      result[3][1][0] += clov[1].off_diag[ 6][0]  * src[2][0][0];
      result[3][1][0] -= clov[1].off_diag[ 6][1]  * src[2][0][1];
      result[3][1][0] += clov[1].off_diag[ 7][0]  * src[2][1][0];
      result[3][1][0] -= clov[1].off_diag[ 7][1]  * src[2][1][1];
      result[3][1][0] += clov[1].off_diag[ 8][0]  * src[2][2][0];
      result[3][1][0] -= clov[1].off_diag[ 8][1]  * src[2][2][1];
      result[3][1][0] += clov[1].off_diag[ 9][0]  * src[3][0][0];
      result[3][1][0] -= clov[1].off_diag[ 9][1]  * src[3][0][1];
      result[3][1][0] += clov[1].off_diag[14][0]  * src[3][2][0];
      result[3][1][0] += clov[1].off_diag[14][1]  * src[3][2][1];

      result[3][1][1]  = clov[1].diag[ 4]  * src[3][1][1];
      result[3][1][1] += clov[1].off_diag[ 6][0]  * src[2][0][1];
      result[3][1][1] += clov[1].off_diag[ 6][1]  * src[2][0][0];
      result[3][1][1] += clov[1].off_diag[ 7][0]  * src[2][1][1];
      result[3][1][1] += clov[1].off_diag[ 7][1]  * src[2][1][0];
      result[3][1][1] += clov[1].off_diag[ 8][0]  * src[2][2][1];
      result[3][1][1] += clov[1].off_diag[ 8][1]  * src[2][2][0];
      result[3][1][1] += clov[1].off_diag[ 9][0]  * src[3][0][1];
      result[3][1][1] += clov[1].off_diag[ 9][1]  * src[3][0][0];
      result[3][1][1] += clov[1].off_diag[14][0]  * src[3][2][1];
      result[3][1][1] -= clov[1].off_diag[14][1]  * src[3][2][0];


      result[3][2][0]  = clov[1].diag[ 5]  * src[3][2][0];
      result[3][2][0] += clov[1].off_diag[10][0]  * src[2][0][0];
      result[3][2][0] -= clov[1].off_diag[10][1]  * src[2][0][1];
      result[3][2][0] += clov[1].off_diag[11][0]  * src[2][1][0];
      result[3][2][0] -= clov[1].off_diag[11][1]  * src[2][1][1];
      result[3][2][0] += clov[1].off_diag[12][0]  * src[2][2][0];
      result[3][2][0] -= clov[1].off_diag[12][1]  * src[2][2][1];
      result[3][2][0] += clov[1].off_diag[13][0]  * src[3][0][0];
      result[3][2][0] -= clov[1].off_diag[13][1]  * src[3][0][1];
      result[3][2][0] += clov[1].off_diag[14][0]  * src[3][1][0];
      result[3][2][0] -= clov[1].off_diag[14][1]  * src[3][1][1];

      result[3][2][1]  = clov[1].diag[ 5]  * src[3][2][1];
      result[3][2][1] += clov[1].off_diag[10][0]  * src[2][0][1];
      result[3][2][1] += clov[1].off_diag[10][1]  * src[2][0][0];
      result[3][2][1] += clov[1].off_diag[11][0]  * src[2][1][1];
      result[3][2][1] += clov[1].off_diag[11][1]  * src[2][1][0];
      result[3][2][1] += clov[1].off_diag[12][0]  * src[2][2][1];
      result[3][2][1] += clov[1].off_diag[12][1]  * src[2][2][0];
      result[3][2][1] += clov[1].off_diag[13][0]  * src[3][0][1];
      result[3][2][1] += clov[1].off_diag[13][1]  * src[3][0][0];
      result[3][2][1] += clov[1].off_diag[14][0]  * src[3][1][1];
      result[3][2][1] += clov[1].off_diag[14][1]  * src[3][1][0];

    }
  }
}

#endif
