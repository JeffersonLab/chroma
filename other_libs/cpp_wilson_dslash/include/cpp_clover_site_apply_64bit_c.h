#ifndef CPP_CLOVER_SITE_APPLY_64_BIT_C_H
#define CPP_CLOVER_SITE_APPLY_64_BIT_C_H

#include <cpp_clover_types.h>
#include <cpp_dslash_types.h>


#undef CMADD
#undef CONJMADD
#undef SMUL

// Complex Multiply add: r += x*y
// Real Part += x.re*y.re - x.im*y.im
// Imag Part += x.re*y.im + x.im*y.re

#define CMADD(r, x, y) { \
    r[0] += x[0]*y[0];  \
    r[0] -= x[1]*y[1];  \
    r[1] += x[0]*y[1];  \
    r[1] += x[1]*y[0];  \
    }

// Complex Conjugate Multiply add: r += conj(x)*y
// conj x is x.re, -x.im
		       // 
// Real Part += x.re*y.re - (-x.im)*y.im = x.re*y.re + x.im*y.im
// Imag Part += x.re*y.im + (-x.im)*y.re = x.re*y.im - x.im*y.re

#define CONJMADD(r,x,y) { \
    r[0] += x[0]*y[0];  \
    r[0] += x[1]*y[1];  \
    r[1] += x[0]*y[1]; \
    r[1] -= x[1]*y[0]; \
  }

// Scalar mult: r = x*y, x is real and y is complex
#define SMUL(r,x,y) { \
    r[0] = x*y[0]; \
    r[1] = x*y[1]; \
  }

namespace CPlusPlusClover { 

  // Use 64bit types
  using namespace Clover64BitTypes;
  using namespace Dslash64BitTypes;

  namespace CPlusPlusClover64Bit { 

    inline 
    void cloverSiteApply(FourSpinor result, const CloverTerm clov, const FourSpinor src) 
    {

      //    result[0][0] =clov[0].diag[ 0]  * src[0][0]
      // +   conj(clov[0].off_diag[ 0]) * src[ 1]
      // +   conj(clov[0].off_diag[ 1]) * src[0][2]
      // +   conj(clov[0].off_diag[ 3]) * src[1][0]
      // +   conj(clov[0].off_diag[ 6]) * src[1][1]
      // +   conj(clov[0].off_diag[10]) * src[1][2];
      SMUL(result[0][0],     clov[0].diag[0],      src[0][0]);
      CONJMADD(result[0][0], clov[0].off_diag[ 0], src[0][1]);
      CONJMADD(result[0][0], clov[0].off_diag[ 1], src[0][2]);
      CONJMADD(result[0][0], clov[0].off_diag[ 3], src[1][0]);
      CONJMADD(result[0][0], clov[0].off_diag[ 6], src[1][1]);
      CONJMADD(result[0][0], clov[0].off_diag[10], src[1][2]);

      // result[0][1] =clov[0].diag[ 1]  * src[0][1]
      // +        clov[0].off_diag[ 0]  * src[0][0]
      // +   conj(clov[0].off_diag[ 2]) * src[0][2]
      // +   conj(clov[0].off_diag[ 4]) * src[1][0]
      // +   conj(clov[0].off_diag[ 7]) * src[1][1]
      // +   conj(clov[0].off_diag[11]) * src[1][2];
      SMUL(result[0][1],     clov[0].diag[1],      src[0][1]);
      CMADD(result[0][1],    clov[0].off_diag[ 0], src[0][0]);
      CONJMADD(result[0][1], clov[0].off_diag[ 2], src[0][2]);
      CONJMADD(result[0][1], clov[0].off_diag[ 4], src[1][0]);
      CONJMADD(result[0][1], clov[0].off_diag[ 7], src[1][1]);
      CONJMADD(result[0][1], clov[0].off_diag[11], src[1][2]);
    
      // result[0][2] =clov[0].diag[ 2]  * src[0][2]
      // +        clov[0].off_diag[ 1]  * src[0][0]
      // +        clov[0].off_diag[ 2]  * src[0][1]
      // +   conj(clov[0].off_diag[ 5]) * src[1][0]
      // +   conj(clov[0].off_diag[ 8]) * src[1][1]
      // +   conj(clov[0].off_diag[12]) * src[1][2];
      SMUL(result[0][2],     clov[0].diag[2],      src[0][2]);
      CMADD(result[0][2],    clov[0].off_diag[ 1], src[0][0]);
      CMADD(result[0][2],    clov[0].off_diag[ 2], src[0][1]);
      CONJMADD(result[0][2], clov[0].off_diag[ 5], src[1][0]);
      CONJMADD(result[0][2], clov[0].off_diag[ 8], src[1][1]);
      CONJMADD(result[0][2], clov[0].off_diag[12], src[1][2]);
    
      // result[1][0] =clov[0].diag[ 3]  * src[1][0]
      // +        clov[0].off_diag[ 3]  * src[0][0]
      // +        clov[0].off_diag[ 4]  * src[0][1]
      // +        clov[0].off_diag[ 5]  * src[0][2]
      // +   conj(clov[0].off_diag[ 9]) * src[1][1]
      // +   conj(clov[0].off_diag[13]) * src[1][2];
      SMUL(result[1][0],     clov[0].diag[3],      src[1][0]);
      CMADD(result[1][0],    clov[0].off_diag[ 3], src[0][0]);
      CMADD(result[1][0],    clov[0].off_diag[ 4], src[0][1]);
      CMADD(result[1][0],    clov[0].off_diag[ 5], src[0][2]);
      CONJMADD(result[1][0], clov[0].off_diag[ 9], src[1][1]);
      CONJMADD(result[1][0], clov[0].off_diag[13], src[1][2]);
    
      // result[1][1] =clov[0].diag[ 4]  * src[1][1]
      // +        clov[0].off_diag[ 6]  * src[0][0]
      // +        clov[0].off_diag[ 7]  * src[0][1]
      // +        clov[0].off_diag[ 8]  * src[0][2]
      // +        clov[0].off_diag[ 9]  * src[1][0]
      // +   conj(clov[0].off_diag[14]) * src[1][2];
      SMUL(result[1][1],     clov[0].diag[4],      src[1][1]);
      CMADD(result[1][1],    clov[0].off_diag[ 6], src[0][0]);
      CMADD(result[1][1],    clov[0].off_diag[ 7], src[0][1]);
      CMADD(result[1][1],    clov[0].off_diag[ 8], src[0][2]);
      CMADD(result[1][1],    clov[0].off_diag[ 9], src[1][0]);
      CONJMADD(result[1][1], clov[0].off_diag[14], src[1][2]);
    
      //    result[1][2] =clov[0].diag[ 5]  * src[1][2]
      // +        clov[0].off_diag[10]  * src[0][0]
      // +        clov[0].off_diag[11]  * src[0][1]
      // +        clov[0].off_diag[12]  * src[0][2]
      // +        clov[0].off_diag[13]  * src[1][0]
      // +        clov[0].off_diag[14]  * src[1][1];
      SMUL(result[1][2],     clov[0].diag[5],      src[1][2]);
      CMADD(result[1][2],    clov[0].off_diag[10], src[0][0]);
      CMADD(result[1][2],    clov[0].off_diag[11], src[0][1]);
      CMADD(result[1][2],    clov[0].off_diag[12], src[0][2]);
      CMADD(result[1][2],    clov[0].off_diag[13], src[1][0]);
      CONJMADD(result[1][2], clov[0].off_diag[14], src[1][1]);
    
 
      // result[2][0] =clov[1].diag[ 0]  * src[2][0]
      // +   conj(clov[1].off_diag[ 0]) * src[2][1]
      // +   conj(clov[1].off_diag[ 1]) * src[2][2]
      // +   conj(clov[1].off_diag[ 3]) * src[3][0]
      // +   conj(clov[1].off_diag[ 6]) * src[3][1]
      // +   conj(clov[1].off_diag[10]) * src[3][2];
      SMUL(result[2][0],     clov[1].diag[0],      src[2][0]);
      CONJMADD(result[2][0], clov[1].off_diag[0],  src[2][1]);
      CONJMADD(result[2][0], clov[1].off_diag[1],  src[2][2]);
      CONJMADD(result[2][0], clov[1].off_diag[3],  src[3][0]);
      CONJMADD(result[2][0], clov[1].off_diag[6],  src[3][1]);
      CONJMADD(result[2][0], clov[1].off_diag[10], src[3][2]);

    
      // result[2][1] =clov[1].diag[ 1]  * src[2][1]
      // +        clov[1].off_diag[ 0]  * src[2][0]
      // +   conj(clov[1].off_diag[ 2]) * src[2][2]
      // +   conj(clov[1].off_diag[ 4]) * src[3][0]
      // +   conj(clov[1].off_diag[ 7]) * src[3][1]
      // +   conj(clov[1].off_diag[11]) * src[3][2];
      SMUL(result[2][1],     clov[1].diag[1],      src[2][1]);
      CMADD(result[2][1],    clov[1].off_diag[0],  src[2][0]);
      CONJMADD(result[2][1], clov[1].off_diag[2],  src[2][2]);
      CONJMADD(result[2][1], clov[1].off_diag[4],  src[3][0]);
      CONJMADD(result[2][1], clov[1].off_diag[7],  src[3][1]);
      CONJMADD(result[2][1], clov[1].off_diag[11], src[3][2]);

    
      //  result[2][2] =clov[1].diag[ 2]  * src[2][2]
      // +        clov[1].off_diag[ 1]  * src[2][0]
      // +        clov[1].off_diag[ 2]  * src[2][1]
      // +   conj(clov[1].off_diag[ 5]) * src[3][0]
      // +   conj(clov[1].off_diag[ 8]) * src[3][1]
      // +   conj(clov[1].off_diag[12]) * src[3][2];
      SMUL(result[2][2],     clov[1].diag[2],      src[2][2]);
      CMADD(result[2][2],    clov[1].off_diag[1],  src[2][0]);
      CMADD(result[2][2],    clov[1].off_diag[2],  src[2][1]);
      CONJMADD(result[2][2], clov[1].off_diag[5],  src[3][0]);
      CONJMADD(result[2][2], clov[1].off_diag[8],  src[3][1]);
      CONJMADD(result[2][2], clov[1].off_diag[12], src[3][2]);

    
      //    result[3][0] =clov[1].diag[ 3]  * src[3][0]
      // +        clov[1].off_diag[ 3]  * src[2][0]
      // +        clov[1].off_diag[ 4]  * src[2][1]
      // +        clov[1].off_diag[ 5]  * src[2][2]
      // +   conj(clov[1].off_diag[ 9]) * src[3][1]
      // +   conj(clov[1].off_diag[13]) * src[3][2];
      SMUL(result[3][0],     clov[1].diag[3],      src[3][0]);
      CMADD(result[3][0],    clov[1].off_diag[3],  src[2][0]);
      CMADD(result[3][0],    clov[1].off_diag[4],  src[2][1]);
      CMADD(result[3][0],    clov[1].off_diag[5],  src[2][2]);
      CONJMADD(result[3][0], clov[1].off_diag[9],  src[3][1]);
      CONJMADD(result[3][0], clov[1].off_diag[13], src[3][2]);

      // result[3][1] =clov[1].diag[ 4]  * src[3][1]
      // +        clov[1].off_diag[ 6]  * src[2][0]
      // +        clov[1].off_diag[ 7]  * src[2][1]
      // +        clov[1].off_diag[ 8]  * src[2][2]
      // +        clov[1].off_diag[ 9]  * src[3][0]
      // +   conj(clov[1].off_diag[14]) * src[3][2];
      SMUL(result[3][1],     clov[1].diag[4],      src[3][1]);
      CMADD(result[3][1],    clov[1].off_diag[6],  src[2][0]);
      CMADD(result[3][1],    clov[1].off_diag[7],  src[2][1]);
      CMADD(result[3][1],    clov[1].off_diag[8],  src[2][2]);
      CMADD(result[3][1],    clov[1].off_diag[9],  src[3][0]);
      CONJMADD(result[3][1], clov[1].off_diag[14], src[3][2]);

    
      // result[3][2] =clov[1].diag[ 5]  * src[3][2]
      // +        clov[1].off_diag[10]  * src[2][0]
      // +        clov[1].off_diag[11]  * src[2][1]
      // +        clov[1].off_diag[12]  * src[2][2]
      // +        clov[1].off_diag[13]  * src[3][0]
      // +        clov[1].off_diag[14]  * src[3][1];
      SMUL(result[3][2],     clov[1].diag[5],      src[3][2]);
      CMADD(result[3][2],    clov[1].off_diag[10],  src[2][0]);
      CMADD(result[3][2],    clov[1].off_diag[11],  src[2][1]);
      CMADD(result[3][2],    clov[1].off_diag[12],  src[2][2]);
      CMADD(result[3][2],    clov[1].off_diag[13],  src[3][0]);
      CMADD(result[3][2],    clov[1].off_diag[14],  src[3][1]);

    }

  } // Namespace
}  // Namespace

#undef CMADD
#undef CONJMADD
#undef SMUL

#endif
