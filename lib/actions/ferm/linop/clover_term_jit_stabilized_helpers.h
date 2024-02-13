//! Exponential of triangular form matrix for spinors

#ifndef __clover_term_jit_stabilized_helpers_h__
#define __clover_term_jit_stabilized_helpers_h__

#include "state.h"

namespace QDP{
    template<typename T> struct PCompJIT;
    template<typename T> struct PTriDiaJIT;
    template<typename T> struct PTriOffJIT;
    template<typename T> struct PCompREG;
    template<typename T> struct PTriDiaREG;
    template<typename T> struct PTriOffREG;
}


namespace Chroma{
  
  template<typename T>
  void tri6_mult(QDP::PCompREG<QDP::PTriDiaREG<QDP::RScalarREG<QDP::WordREG<T>>>>& ADia, 
                 QDP::PCompREG<QDP::PTriOffREG<QDP::RComplexREG<QDP::WordREG<T>>>>& AOff,
                 QDP::PCompREG<QDP::PTriDiaREG<QDP::RScalarREG<QDP::WordREG<T>>>>& BDia, 
                 QDP::PCompREG<QDP::PTriOffREG<QDP::RComplexREG<QDP::WordREG<T>>>>& BOff,
                 QDP::PCompREG<QDP::PTriDiaREG<QDP::RScalarREG<QDP::WordREG<T>>>>& CDia, 
                 QDP::PCompREG<QDP::PTriOffREG<QDP::RComplexREG<QDP::WordREG<T>>>>& COff,
                 int i){
    // Diagonal Terms:

    ADia.elem(i).elem(0) = BDia.elem(i).elem(0) * CDia.elem(i).elem(0)
                         + real(BOff.elem(i).elem(0) ) * real(COff.elem(i).elem(0) ) + imag(BOff.elem(i).elem(0) ) * imag(COff.elem(i).elem(0) ) 
                         + real(BOff.elem(i).elem(1) ) * real(COff.elem(i).elem(1) ) + imag(BOff.elem(i).elem(1) ) * imag(COff.elem(i).elem(1) ) 
                         + real(BOff.elem(i).elem(3) ) * real(COff.elem(i).elem(3) ) + imag(BOff.elem(i).elem(3) ) * imag(COff.elem(i).elem(3) ) 
                         + real(BOff.elem(i).elem(6) ) * real(COff.elem(i).elem(6) ) + imag(BOff.elem(i).elem(6) ) * imag(COff.elem(i).elem(6) ) 
                         + real(BOff.elem(i).elem(10)) * real(COff.elem(i).elem(10)) + imag(BOff.elem(i).elem(10)) * imag(COff.elem(i).elem(10));
    
    ADia.elem(i).elem(1) = BDia.elem(i).elem(1) * CDia.elem(i).elem(1) 
                         + real(BOff.elem(i).elem(0) ) * real(COff.elem(i).elem(0) ) + imag(BOff.elem(i).elem(0) ) * imag(COff.elem(i).elem(0) ) 
                         + real(BOff.elem(i).elem(2) ) * real(COff.elem(i).elem(2) ) + imag(BOff.elem(i).elem(2) ) * imag(COff.elem(i).elem(2) ) 
                         + real(BOff.elem(i).elem(4) ) * real(COff.elem(i).elem(4) ) + imag(BOff.elem(i).elem(4) ) * imag(COff.elem(i).elem(4) ) 
                         + real(BOff.elem(i).elem(7) ) * real(COff.elem(i).elem(7) ) + imag(BOff.elem(i).elem(7) ) * imag(COff.elem(i).elem(7) ) 
                         + real(BOff.elem(i).elem(11)) * real(COff.elem(i).elem(11)) + imag(BOff.elem(i).elem(11)) * imag(COff.elem(i).elem(11));

    ADia.elem(i).elem(2) = BDia.elem(i).elem(2) * CDia.elem(i).elem(2) 
                         + real(BOff.elem(i).elem(1) ) * real(COff.elem(i).elem(1) ) + imag(BOff.elem(i).elem(1) ) * imag(COff.elem(i).elem(1) ) 
                         + real(BOff.elem(i).elem(2) ) * real(COff.elem(i).elem(2) ) + imag(BOff.elem(i).elem(2) ) * imag(COff.elem(i).elem(2) ) 
                         + real(BOff.elem(i).elem(5) ) * real(COff.elem(i).elem(5) ) + imag(BOff.elem(i).elem(5) ) * imag(COff.elem(i).elem(5) ) 
                         + real(BOff.elem(i).elem(8) ) * real(COff.elem(i).elem(8) ) + imag(BOff.elem(i).elem(8) ) * imag(COff.elem(i).elem(8) ) 
                         + real(BOff.elem(i).elem(12)) * real(COff.elem(i).elem(12)) + imag(BOff.elem(i).elem(12)) * imag(COff.elem(i).elem(12));

    ADia.elem(i).elem(3) = BDia.elem(i).elem(3) * CDia.elem(i).elem(3) 
                         + real(BOff.elem(i).elem(3) ) * real(COff.elem(i).elem(3) ) + imag(BOff.elem(i).elem(3) ) * imag(COff.elem(i).elem(3) ) 
                         + real(BOff.elem(i).elem(4) ) * real(COff.elem(i).elem(4) ) + imag(BOff.elem(i).elem(4) ) * imag(COff.elem(i).elem(4) ) 
                         + real(BOff.elem(i).elem(5) ) * real(COff.elem(i).elem(5) ) + imag(BOff.elem(i).elem(5) ) * imag(COff.elem(i).elem(5) ) 
                         + real(BOff.elem(i).elem(9) ) * real(COff.elem(i).elem(9) ) + imag(BOff.elem(i).elem(9) ) * imag(COff.elem(i).elem(9) ) 
                         + real(BOff.elem(i).elem(13)) * real(COff.elem(i).elem(13)) + imag(BOff.elem(i).elem(13)) * imag(COff.elem(i).elem(13));

    ADia.elem(i).elem(4) = BDia.elem(i).elem(4) * CDia.elem(i).elem(4) 
                         + real(BOff.elem(i).elem(6) ) * real(COff.elem(i).elem(6) ) + imag(BOff.elem(i).elem(6) ) * imag(COff.elem(i).elem(6) ) 
                         + real(BOff.elem(i).elem(7) ) * real(COff.elem(i).elem(7) ) + imag(BOff.elem(i).elem(7) ) * imag(COff.elem(i).elem(7) ) 
                         + real(BOff.elem(i).elem(8) ) * real(COff.elem(i).elem(8) ) + imag(BOff.elem(i).elem(8) ) * imag(COff.elem(i).elem(8) ) 
                         + real(BOff.elem(i).elem(9) ) * real(COff.elem(i).elem(9) ) + imag(BOff.elem(i).elem(9) ) * imag(COff.elem(i).elem(9) ) 
                         + real(BOff.elem(i).elem(14)) * real(COff.elem(i).elem(14)) + imag(BOff.elem(i).elem(14)) * imag(COff.elem(i).elem(14));

    ADia.elem(i).elem(5) = BDia.elem(i).elem(5) * CDia.elem(i).elem(5) 
                         + real(BOff.elem(i).elem(10)) * real(COff.elem(i).elem(10)) + imag(BOff.elem(i).elem(10)) * imag(COff.elem(i).elem(10)) 
                         + real(BOff.elem(i).elem(11)) * real(COff.elem(i).elem(11)) + imag(BOff.elem(i).elem(11)) * imag(COff.elem(i).elem(11)) 
                         + real(BOff.elem(i).elem(12)) * real(COff.elem(i).elem(12)) + imag(BOff.elem(i).elem(12)) * imag(COff.elem(i).elem(12)) 
                         + real(BOff.elem(i).elem(13)) * real(COff.elem(i).elem(13)) + imag(BOff.elem(i).elem(13)) * imag(COff.elem(i).elem(13)) 
                         + real(BOff.elem(i).elem(14)) * real(COff.elem(i).elem(14)) + imag(BOff.elem(i).elem(14)) * imag(COff.elem(i).elem(14));     

    // Off-diagonal Terms:
    AOff.elem(i).elem(0) = BDia.elem(i).elem(0) *      COff.elem(i).elem(0)
                         + BOff.elem(i).elem(0) *      CDia.elem(i).elem(1)
                         + BOff.elem(i).elem(1) * conj(COff.elem(i).elem(2))
                         + BOff.elem(i).elem(3) * conj(COff.elem(i).elem(4))
                         + BOff.elem(i).elem(6) * conj(COff.elem(i).elem(7))
                         + BOff.elem(i).elem(10)* conj(COff.elem(i).elem(11));
    // A.offd[i][0] = B.diag[i][0]  *      C.offd[i][0]  
    //              + B.offd[i][0]  *      C.diag[i][1]
    //              + B.offd[i][1]  * conj(C.offd[i][2])  
    //              + B.offd[i][3]  * conj(C.offd[i][4])
    //              + B.offd[i][6]  * conj(C.offd[i][7])  
    //              + B.offd[i][10] * conj(C.offd[i][11]); 

    AOff.elem(i).elem(1) = BDia.elem(i).elem(0) *      COff.elem(i).elem(1)
                         + BOff.elem(i).elem(0) *      COff.elem(i).elem(2)
                         + BOff.elem(i).elem(1) *      CDia.elem(i).elem(2)
                         + BOff.elem(i).elem(3) * conj(COff.elem(i).elem(5))
                         + BOff.elem(i).elem(6) * conj(COff.elem(i).elem(8))
                         + BOff.elem(i).elem(10)* conj(COff.elem(i).elem(12));
    // A.offd[i][1] = B.diag[i][0]  *      C.offd[i][1]  
    //              + B.offd[i][0]  *      C.offd[i][2]
    //              + B.offd[i][1]  *      C.diag[i][2]  
    //              + B.offd[i][3]  * conj(C.offd[i][5])
    //              + B.offd[i][6]  * conj(C.offd[i][8]) 
    //              + B.offd[i][10] * conj(C.offd[i][12]); 

    AOff.elem(i).elem(2) = conj(BOff.elem(i).elem(0)) *      COff.elem(i).elem(1)
                         +      BDia.elem(i).elem(1)  *      COff.elem(i).elem(2)
                         +      BOff.elem(i).elem(2)  *      CDia.elem(i).elem(2)
                         +      BOff.elem(i).elem(4)  * conj(COff.elem(i).elem(5))
                         +      BOff.elem(i).elem(7)  * conj(COff.elem(i).elem(8))
                         +      BOff.elem(i).elem(11) * conj(COff.elem(i).elem(12));
    // A.offd[i][2] = conj(B.offd[i][0])  *      C.offd[i][1]  
    //              +      B.diag[i][1]   *      C.offd[i][2]
    //              +      B.offd[i][2]   *      C.diag[i][2]  
    //              +      B.offd[i][4]   * conj(C.offd[i][5])
    //              +      B.offd[i][7]   * conj(C.offd[i][8]) 
    //              +      B.offd[i][11]  * conj(C.offd[i][12]);   

    AOff.elem(i).elem(3) = BDia.elem(i).elem(0) * COff.elem(i).elem(3)
                         + BOff.elem(i).elem(0) * COff.elem(i).elem(4)
                         + BOff.elem(i).elem(1) * COff.elem(i).elem(5)
                         + BOff.elem(i).elem(3) * CDia.elem(i).elem(3)
                         + BOff.elem(i).elem(6) * conj(COff.elem(i).elem(9))
                         + BOff.elem(i).elem(10)* conj(COff.elem(i).elem(13));
    // A.offd[i][3] = B.diag[i][0]  *      C.offd[i][3]  
    //              + B.offd[i][0]  *      C.offd[i][4]
    //              + B.offd[i][1]  *      C.offd[i][5]  
    //              + B.offd[i][3]  *      C.diag[i][3]
    //              + B.offd[i][6]  * conj(C.offd[i][9]) 
    //              + B.offd[i][10] * conj(C.offd[i][13]); 

    AOff.elem(i).elem(4) = conj(BOff.elem(i).elem(0)) *      COff.elem(i).elem(3)
                         +      BDia.elem(i).elem(1)  *      COff.elem(i).elem(4)
                         +      BOff.elem(i).elem(2)  *      COff.elem(i).elem(5)
                         +      BOff.elem(i).elem(4)  *      CDia.elem(i).elem(3)
                         +      BOff.elem(i).elem(7)  * conj(COff.elem(i).elem(9))
                         +      BOff.elem(i).elem(11) * conj(COff.elem(i).elem(13));
    // A.offd[i][4] = conj(B.offd[i][0])  *      C.offd[i][3]  
    //                   + B.diag[i][1]   *      C.offd[i][4]
    //                   + B.offd[i][2]   *      C.offd[i][5]  
    //                   + B.offd[i][4]   *      C.diag[i][3]
    //                   + B.offd[i][7]   * conj(C.offd[i][9]) 
    //                   + B.offd[i][11]  * conj(C.offd[i][13]);  

    AOff.elem(i).elem(5) = conj(BOff.elem(i).elem(1)) *      COff.elem(i).elem(3)
                         + conj(BOff.elem(i).elem(2)) *      COff.elem(i).elem(4)
                         +      BDia.elem(i).elem(2)  *      COff.elem(i).elem(5)
                         +      BOff.elem(i).elem(5)  *      CDia.elem(i).elem(3)
                         +      BOff.elem(i).elem(8)  * conj(COff.elem(i).elem(9))
                         +      BOff.elem(i).elem(12) * conj(COff.elem(i).elem(13));
    // A.offd[i][5]  = conj(B.offd[i][1])  *      C.offd[i][3]  
    //               + conj(B.offd[i][2])  *      C.offd[i][4]
    //               +      B.diag[i][2]   *      C.offd[i][5]  
    //               +      B.offd[i][5]   *      C.diag[i][3]
    //               +      B.offd[i][8]   * conj(C.offd[i][9]) 
    //               +      B.offd[i][12]  * conj(C.offd[i][13]);   

    AOff.elem(i).elem(6) = BDia.elem(i).elem(0)  *      COff.elem(i).elem(6)
                         + BOff.elem(i).elem(0)  *      COff.elem(i).elem(7)
                         + BOff.elem(i).elem(1)  *      COff.elem(i).elem(8)
                         + BOff.elem(i).elem(3)  *      COff.elem(i).elem(9)
                         + BOff.elem(i).elem(6)  *      CDia.elem(i).elem(4)
                         + BOff.elem(i).elem(10) * conj(COff.elem(i).elem(14));
    // A.offd[i][6] = B.diag[i][0]  *      C.offd[i][6]  
    //              + B.offd[i][0]  *      C.offd[i][7]
    //              + B.offd[i][1]  *      C.offd[i][8]  
    //              + B.offd[i][3]  *      C.offd[i][9]
    //              + B.offd[i][6]  *      C.diag[i][4]  
    //              + B.offd[i][10] * conj(C.offd[i][14]); 

    AOff.elem(i).elem(7) = conj(BOff.elem(i).elem(0)) *      COff.elem(i).elem(6)
                         +      BDia.elem(i).elem(1)  *      COff.elem(i).elem(7)
                         +      BOff.elem(i).elem(2)  *      COff.elem(i).elem(8)
                         +      BOff.elem(i).elem(4)  *      COff.elem(i).elem(9)
                         +      BOff.elem(i).elem(7)  *      CDia.elem(i).elem(4)
                         +      BOff.elem(i).elem(11) * conj(COff.elem(i).elem(14));
    // A.offd[i][7] = conj(B.offd[i][0])  *      C.offd[i][6]  
    //              +      B.diag[i][1]   *      C.offd[i][7]
    //              +      B.offd[i][2]   *      C.offd[i][8]  
    //              +      B.offd[i][4]   *      C.offd[i][9]
    //              +      B.offd[i][7]   *      C.diag[i][4]  
    //              +      B.offd[i][11]  * conj(C.offd[i][14]);   

    AOff.elem(i).elem(8) = conj(BOff.elem(i).elem(1)) *      COff.elem(i).elem(6)
                         + conj(BOff.elem(i).elem(2)) *      COff.elem(i).elem(7)
                         +      BDia.elem(i).elem(2)  *      COff.elem(i).elem(8)
                         +      BOff.elem(i).elem(5)  *      COff.elem(i).elem(9)
                         +      BOff.elem(i).elem(8)  *      CDia.elem(i).elem(4)
                         +      BOff.elem(i).elem(12) * conj(COff.elem(i).elem(14));
    // A.offd[i][8] = conj(B.offd[i][1])  *      C.offd[i][6]  
    //              + conj(B.offd[i][2])  *      C.offd[i][7]
    //              +      B.diag[i][2]   *      C.offd[i][8]  
    //              +      B.offd[i][5]   *      C.offd[i][9]
    //              +      B.offd[i][8]   *      C.diag[i][4]  
    //              +      B.offd[i][12]  * conj(C.offd[i][14]);   
 
    AOff.elem(i).elem(9) = conj(BOff.elem(i).elem(3)) *      COff.elem(i).elem(6)
                         + conj(BOff.elem(i).elem(4)) *      COff.elem(i).elem(7)
                         + conj(BOff.elem(i).elem(5)) *      COff.elem(i).elem(8)
                         +      BDia.elem(i).elem(3)  *      COff.elem(i).elem(9)
                         +      BOff.elem(i).elem(9)  *      CDia.elem(i).elem(4)
                         +      BOff.elem(i).elem(13) * conj(COff.elem(i).elem(14));
    // A.offd[i][9] = conj(B.offd[i][3])  *      C.offd[i][6]  
    //              + conj(B.offd[i][4])  *      C.offd[i][7]
    //              + conj(B.offd[i][5])  *      C.offd[i][8]  
    //              +      B.diag[i][3]   *      C.offd[i][9]
    //              +      B.offd[i][9]   *      C.diag[i][4]  
    //              +      B.offd[i][13]  * conj(C.offd[i][14]);   

    AOff.elem(i).elem(10) = BDia.elem(i).elem(0)  * COff.elem(i).elem(10)
                          + BOff.elem(i).elem(0)  * COff.elem(i).elem(11)
                          + BOff.elem(i).elem(1)  * COff.elem(i).elem(12)
                          + BOff.elem(i).elem(3)  * COff.elem(i).elem(13)
                          + BOff.elem(i).elem(6)  * COff.elem(i).elem(14)
                          + BOff.elem(i).elem(10) * CDia.elem(i).elem(5);
    // A.offd[i][10] = B.diag[i][0]  * C.offd[i][10]  
    //               + B.offd[i][0]  * C.offd[i][11]
    //               + B.offd[i][1]  * C.offd[i][12]  
    //               + B.offd[i][3]  * C.offd[i][13]
    //               + B.offd[i][6]  * C.offd[i][14]  
    //               + B.offd[i][10] * C.diag[i][5];      

    AOff.elem(i).elem(11) = conj(BOff.elem(i).elem(0)) * COff.elem(i).elem(10)
                          +      BDia.elem(i).elem(1)  * COff.elem(i).elem(11)
                          +      BOff.elem(i).elem(2)  * COff.elem(i).elem(12)
                          +      BOff.elem(i).elem(4)  * COff.elem(i).elem(13)
                          +      BOff.elem(i).elem(7)  * COff.elem(i).elem(14)
                          +      BOff.elem(i).elem(11) * CDia.elem(i).elem(5);
    // A.offd[i][11] = conj(B.offd[i][0])  * C.offd[i][10]  
    //               +      B.diag[i][1]   * C.offd[i][11]
    //               +      B.offd[i][2]   * C.offd[i][12]  
    //               +      B.offd[i][4]   * C.offd[i][13]
    //               +      B.offd[i][7]   * C.offd[i][14]  
    //               +      B.offd[i][11]  * C.diag[i][5];      

    AOff.elem(i).elem(12) = conj(BOff.elem(i).elem(1)) * COff.elem(i).elem(10)
                          + conj(BOff.elem(i).elem(2)) * COff.elem(i).elem(11)
                          +      BDia.elem(i).elem(2)  * COff.elem(i).elem(12)
                          +      BOff.elem(i).elem(5)  * COff.elem(i).elem(13)
                          +      BOff.elem(i).elem(8)  * COff.elem(i).elem(14)
                          +      BOff.elem(i).elem(12) * CDia.elem(i).elem(5);
    // A.offd[i][12] = conj(B.offd[i][1])  * C.offd[i][10]  
    //               + conj(B.offd[i][2])  * C.offd[i][11]
    //               +      B.diag[i][2]   * C.offd[i][12]  
    //               +      B.offd[i][5]   * C.offd[i][13]
    //               +      B.offd[i][8]   * C.offd[i][14]  
    //               +      B.offd[i][12]  * C.diag[i][5];    

    AOff.elem(i).elem(13) = conj(BOff.elem(i).elem(3)) * COff.elem(i).elem(10)
                          + conj(BOff.elem(i).elem(4)) * COff.elem(i).elem(11)
                          + conj(BOff.elem(i).elem(5)) * COff.elem(i).elem(12)
                          +      BDia.elem(i).elem(3)  * COff.elem(i).elem(13)
                          +      BOff.elem(i).elem(9)  * COff.elem(i).elem(14)
                          +      BOff.elem(i).elem(13) * CDia.elem(i).elem(5);
    // A.offd[i][13] = conj(B.offd[i][3])  * C.offd[i][10]  
    //               + conj(B.offd[i][4])  * C.offd[i][11]
    //               + conj(B.offd[i][5])  * C.offd[i][12]  
    //               +      B.diag[i][3]   * C.offd[i][13]
    //               +      B.offd[i][9]   * C.offd[i][14]  
    //               +      B.offd[i][13]  * C.diag[i][5];   

    AOff.elem(i).elem(14) = conj(BOff.elem(i).elem(6)) * COff.elem(i).elem(10)
                          + conj(BOff.elem(i).elem(7)) * COff.elem(i).elem(11)
                          + conj(BOff.elem(i).elem(8)) * COff.elem(i).elem(12)
                          + conj(BOff.elem(i).elem(9)) * COff.elem(i).elem(13)
                          +      BDia.elem(i).elem(4)  * COff.elem(i).elem(14)
                          +      BOff.elem(i).elem(14) * CDia.elem(i).elem(5);
    // A.offd[i][14] = conj(B.offd[i][6])  * C.offd[i][10]  
    //               + conj(B.offd[i][7])  * C.offd[i][11]
    //               + conj(B.offd[i][8])  * C.offd[i][12]  
    //               + conj(B.offd[i][9])  * C.offd[i][13]
    //               +      B.diag[i][4]   * C.offd[i][14]  
    //               +      B.offd[i][14]  * C.diag[i][5];    
  
  } 

  //! Computes Tr(A)
  template<typename T>
  RScalarREG<WordREG<T> > tri6_trace( QDP::PCompREG<QDP::PTriDiaREG<QDP::RScalarREG<QDP::WordREG<T>>>>& ADia, 
                         QDP::PCompREG<QDP::PTriOffREG<QDP::RComplexREG<QDP::WordREG<T>>>>& AOff,
                        int comp){
    RScalarREG<WordREG<T> > tr = 0.0;
    for (size_t i = 0; i < 2*Nc; i++)
        tr += ADia.elem(comp).elem(i);
    return tr;
  }

  // Computes Tr(A*B) using hermiticity of the matrices
  template<typename T>
  RScalarREG<WordREG<T> >  tri6_trace_mul( QDP::PCompREG<QDP::PTriDiaREG<QDP::RScalarREG<QDP::WordREG<T>>>>& ADia, 
                             QDP::PCompREG<QDP::PTriOffREG<QDP::RComplexREG<QDP::WordREG<T>>>>& AOff,
                             QDP::PCompREG<QDP::PTriDiaREG<QDP::RScalarREG<QDP::WordREG<T>>>>& BDia, 
                             QDP::PCompREG<QDP::PTriOffREG<QDP::RComplexREG<QDP::WordREG<T>>>>& BOff,
                            int comp){
    RScalarREG<WordREG<T> > trd = 0.0, tro = 0.0;
    for (size_t i = 0; i < 2*Nc; i++)
        trd += ADia.elem(comp).elem(i) * BDia.elem(comp).elem(i);

    for (size_t i = 0; i < 2*Nc*Nc-Nc; i++)
        tro += real(AOff.elem(comp).elem(i)) * real(BOff.elem(comp).elem(i))
             + imag(AOff.elem(comp).elem(i)) * imag(BOff.elem(comp).elem(i));
    return trd + tro + tro;
  }   

  template<typename T>
  void exponentiate(QDP::PCompJIT<QDP::PTriDiaJIT<QDP::RScalarJIT<QDP::WordJIT<T>>>>& ADia_Jit, 
                    QDP::PCompJIT<QDP::PTriOffJIT<QDP::RComplexJIT<QDP::WordJIT<T>>>>& AOff_Jit, 
                    int sign){
    const int Niter = 17;
    QDP::PCompREG<QDP::PTriDiaREG<QDP::RScalarREG<QDP::WordREG<T>>>> ADia, A2Dia, A3Dia, tmpDia;
    QDP::PCompREG<QDP::PTriOffREG<QDP::RComplexREG<QDP::WordREG<T>>>> AOff, A2Off, A3Off, tmpOff;

    for(size_t comp = 0; comp < 2; comp++){
        for(size_t i = 0; i < 2*Nc; i++){
            ADia.elem(comp).elem(i) = ADia_Jit.elem(comp).elem(i);
        }
        for(size_t i = 0; i < 2*Nc*Nc-Nc; i++){
            AOff.elem(comp).elem(i) = AOff_Jit.elem(comp).elem(i);
        }
    }

    for(size_t comp = 0; comp < 2; comp++){

      tri6_mult(A2Dia, A2Off, ADia, AOff, ADia, AOff, comp);
      tri6_mult(A3Dia, A3Off, ADia, AOff, A2Dia, A2Off, comp);
      
      RScalarREG<WordREG<T>> tr[7];
      tr[0] = 0;
      tr[1] = 0;
      tr[2] = tri6_trace(A2Dia, A2Off, comp);
      tr[3] = tri6_trace(A3Dia, A3Off, comp);
      tr[4] = tri6_trace_mul(A2Dia, A2Off, A2Dia, A2Off, comp);
      tr[5] = tri6_trace_mul(A2Dia, A2Off, A3Dia, A3Off, comp);
      tr[6] = tri6_trace_mul(A3Dia, A3Off, A3Dia, A3Off, comp);

      RScalarREG<WordREG<T>> p[5];
      p[0] = RScalarREG<WordREG<T> >(-1.0/6.0)*tr[6] + RScalarREG<WordREG<T> >(1.0/18.0)*tr[3]*tr[3] 
           - RScalarREG<WordREG<T> >(1.0/48.0)*tr[2]*tr[2]*tr[2] + RScalarREG<WordREG<T> >(1.0/8.0)*tr[4]*tr[2];
      p[1] = RScalarREG<WordREG<T> >(-1.0/5.0)*tr[5] + RScalarREG<WordREG<T> >(1.0/6.0) *tr[2]*tr[3];
      p[2] = RScalarREG<WordREG<T> >(-1.0/4.0)*tr[4] + RScalarREG<WordREG<T> >(1.0/8.0) *tr[2]*tr[2];
      p[3] = RScalarREG<WordREG<T> >(-1.0/3.0)*tr[3];
      p[4] = RScalarREG<WordREG<T> >(-1.0/2.0)*tr[2];

      RScalarREG<WordREG<T>> c[Niter+1];
      c[0] = 1;
      for (size_t i = 1; i < Niter+1; i++)
          c[i] = c[i - 1] / RScalarREG<WordREG<T> >(i);
      if (sign == 1)
          for (size_t i = 1; i < Niter+1; i += 2)
              c[i] = -c[i];

      RScalarREG<WordREG<T>> q[6];    
      if (Niter > 5){ 
          int ic = Niter - 6;
          for (size_t i = 0; i < 6; i++)
              q[i] = c[ic+1+i];
          
          while (ic >= 0) {
              RScalarREG<WordREG<T> > q5 = q[5];
              q[5] = q[4];
              q[4] = q[3]  - q5 * p[4];
              q[3] = q[2]  - q5 * p[3];
              q[2] = q[1]  - q5 * p[2];
              q[1] = q[0]  - q5 * p[1];
              q[0] = c[ic] - q5 * p[0];
              ic -= 1;
          }
      }
      else {
          QDPIO::cerr << "error" << std::endl;
      }

      //   I*q0 + A*q1 + A^2*q2 + A^3*q3
      for (size_t i = 0; i < 2*Nc; i++)
        ADia.elem(comp).elem(i) = q[0] 
                                + ADia.elem(comp).elem(i)*q[1] 
                                + A2Dia.elem(comp).elem(i)*q[2]
                                + A3Dia.elem(comp).elem(i)*q[3];
        //   A.diag[comp][i] = q[0] + A.diag[comp][i]*q[1] + A2.diag[comp][i]*q[2] + A3.diag[comp][i]*q[3];   
      for (size_t i = 0; i < 2*Nc*Nc-Nc; i++)
        AOff.elem(comp).elem(i) = AOff.elem(comp).elem(i)*q[1] 
                                + A2Off.elem(comp).elem(i)*q[2]
                                + A3Off.elem(comp).elem(i)*q[3];
        //   A.offd[comp][i] = A.offd[comp][i]*q[1] + A2.offd[comp][i]*q[2] + A3.offd[comp][i]*q[3];

      // A^4*q4
      tri6_mult(tmpDia, tmpOff, A2Dia, A2Off, A2Dia, A2Off, comp);
      for (size_t i = 0; i < 2*Nc; i++)
        ADia.elem(comp).elem(i) += tmpDia.elem(comp).elem(i)*q[4];
        //   A.diag[comp][i] += tmp.diag[comp][i]*q[4];   
      for (size_t i = 0; i < 2*Nc*Nc-Nc; i++)
        AOff.elem(comp).elem(i) += tmpOff.elem(comp).elem(i)*q[4];
        //   A.offd[comp][i] += tmp.offd[comp][i]*q[4];

      // A^5*q5
      tri6_mult(tmpDia, tmpOff, A3Dia, A3Off, A2Dia, A2Off, comp);
      for (size_t i = 0; i < 2*Nc; i++)
        ADia.elem(comp).elem(i) += tmpDia.elem(comp).elem(i)*q[5]; 
        //   A.diag[comp][i] += tmp.diag[comp][i]*q[5];   
      for (size_t i = 0; i < 2*Nc*Nc-Nc; i++)
        AOff.elem(comp).elem(i) += tmpOff.elem(comp).elem(i)*q[5]; 
        //  A.offd[comp][i] += tmp.offd[comp][i]*q[5];

        
    }

    for(size_t comp = 0; comp < 2; comp++){
        for(size_t i = 0; i < 2*Nc; i++){
            ADia_Jit.elem(comp).elem(i) = ADia.elem(comp).elem(i);
        }
        for(size_t i = 0; i < 2*Nc*Nc-Nc; i++){
            AOff_Jit.elem(comp).elem(i) = AOff.elem(comp).elem(i);
        }
    }
  }   
}

#endif
