//! Exponential of triangular form matrix for spinors

#ifndef __clover_term_qdp_stabilized_helpers_h__
#define __clover_term_qdp_stabilized_helpers_h__

#include "state.h"
#include "qdp_allocator.h" 

#include <complex>
#include <array>

namespace Chroma{
    
  template<typename R> struct PrimitiveClovTriang;
  
  template<typename R>
  void tri6_mult(PrimitiveClovTriang<R>& A, const PrimitiveClovTriang<R>& B, const PrimitiveClovTriang<R>& C, int i){
    // Diagonal Terms:
    A.diag[i][0] = B.diag[i][0].elem()  * C.diag[i][0].elem()  
              + B.offd[i][0].real()  * C.offd[i][0].real()  + B.offd[i][0].imag()  * C.offd[i][0].imag()
              + B.offd[i][1].real()  * C.offd[i][1].real()  + B.offd[i][1].imag()  * C.offd[i][1].imag() 
              + B.offd[i][3].real()  * C.offd[i][3].real()  + B.offd[i][3].imag()  * C.offd[i][3].imag() 
              + B.offd[i][6].real()  * C.offd[i][6].real()  + B.offd[i][6].imag()  * C.offd[i][6].imag() 
              + B.offd[i][10].real()  * C.offd[i][10].real()  + B.offd[i][10].imag()  * C.offd[i][10].imag(); 

    A.diag[i][1] = B.diag[i][1].elem()  * C.diag[i][1].elem()  
              + B.offd[i][0].real()  * C.offd[i][0].real()  + B.offd[i][0].imag()  * C.offd[i][0].imag()
              + B.offd[i][2].real()  * C.offd[i][2].real()  + B.offd[i][2].imag()  * C.offd[i][2].imag() 
              + B.offd[i][4].real()  * C.offd[i][4].real()  + B.offd[i][4].imag()  * C.offd[i][4].imag() 
              + B.offd[i][7].real()  * C.offd[i][7].real()  + B.offd[i][7].imag()  * C.offd[i][7].imag() 
              + B.offd[i][11].real()  * C.offd[i][11].real()  + B.offd[i][11].imag()  * C.offd[i][11].imag(); 

    A.diag[i][2] = B.diag[i][2].elem()  * C.diag[i][2].elem()  
              + B.offd[i][1].real()  * C.offd[i][1].real()  + B.offd[i][1].imag()  * C.offd[i][1].imag()
              + B.offd[i][2].real()  * C.offd[i][2].real()  + B.offd[i][2].imag()  * C.offd[i][2].imag() 
              + B.offd[i][5].real()  * C.offd[i][5].real()  + B.offd[i][5].imag()  * C.offd[i][5].imag() 
              + B.offd[i][8].real() * C.offd[i][8].real() + B.offd[i][8].imag() * C.offd[i][8].imag() 
              + B.offd[i][12].real() * C.offd[i][12].real() + B.offd[i][12].imag() * C.offd[i][12].imag(); 

    A.diag[i][3] = B.diag[i][3].elem()  * C.diag[i][3].elem()  
              + B.offd[i][3].real()  * C.offd[i][3].real()  + B.offd[i][3].imag()  * C.offd[i][3].imag()
              + B.offd[i][4].real()  * C.offd[i][4].real()  + B.offd[i][4].imag()  * C.offd[i][4].imag() 
              + B.offd[i][5].real()  * C.offd[i][5].real()  + B.offd[i][5].imag()  * C.offd[i][5].imag() 
              + B.offd[i][9].real() * C.offd[i][9].real() + B.offd[i][9].imag() * C.offd[i][9].imag() 
              + B.offd[i][13].real() * C.offd[i][13].real() + B.offd[i][13].imag() * C.offd[i][13].imag(); 

    A.diag[i][4] = B.diag[i][4].elem()  * C.diag[i][4].elem()  
              + B.offd[i][6].real()  * C.offd[i][6].real()  + B.offd[i][6].imag()  * C.offd[i][6].imag()
              + B.offd[i][7].real()  * C.offd[i][7].real()  + B.offd[i][7].imag()  * C.offd[i][7].imag() 
              + B.offd[i][8].real() * C.offd[i][8].real() + B.offd[i][8].imag() * C.offd[i][8].imag() 
              + B.offd[i][9].real() * C.offd[i][9].real() + B.offd[i][9].imag() * C.offd[i][9].imag() 
              + B.offd[i][14].real() * C.offd[i][14].real() + B.offd[i][14].imag() * C.offd[i][14].imag(); 

    A.diag[i][5] = B.diag[i][5].elem()  * C.diag[i][5].elem()  
              + B.offd[i][10].real()  * C.offd[i][10].real()  + B.offd[i][10].imag()  * C.offd[i][10].imag()
              + B.offd[i][11].real()  * C.offd[i][11].real()  + B.offd[i][11].imag()  * C.offd[i][11].imag() 
              + B.offd[i][12].real() * C.offd[i][12].real() + B.offd[i][12].imag() * C.offd[i][12].imag() 
              + B.offd[i][13].real() * C.offd[i][13].real() + B.offd[i][13].imag() * C.offd[i][13].imag() 
              + B.offd[i][14].real() * C.offd[i][14].real() + B.offd[i][14].imag() * C.offd[i][14].imag();                 

    // Off-diagonal Terms:
    A.offd[i][0] = B.diag[i][0]  * C.offd[i][0]  + B.offd[i][0]  * C.diag[i][1]
              + B.offd[i][1]  * conj(C.offd[i][2])  + B.offd[i][3]  * conj(C.offd[i][4])
              + B.offd[i][6]  * conj(C.offd[i][7])  + B.offd[i][10]  * conj(C.offd[i][11]); 

    A.offd[i][1] = B.diag[i][0]  * C.offd[i][1]  + B.offd[i][0]  * C.offd[i][2]
              + B.offd[i][1]  * C.diag[i][2]  + B.offd[i][3]  * conj(C.offd[i][5])
              + B.offd[i][6]  * conj(C.offd[i][8]) + B.offd[i][10]  * conj(C.offd[i][12]); 

    A.offd[i][2] = conj(B.offd[i][0])  * C.offd[i][1]  + B.diag[i][1]  * C.offd[i][2]
              + B.offd[i][2]  * C.diag[i][2]  + B.offd[i][4]  * conj(C.offd[i][5])
              + B.offd[i][7]  * conj(C.offd[i][8]) + B.offd[i][11]  * conj(C.offd[i][12]);   

    A.offd[i][3] = B.diag[i][0]  * C.offd[i][3]  + B.offd[i][0]  * C.offd[i][4]
              + B.offd[i][1]  * C.offd[i][5]  + B.offd[i][3]  * C.diag[i][3]
              + B.offd[i][6]  * conj(C.offd[i][9]) + B.offd[i][10]  * conj(C.offd[i][13]); 

    A.offd[i][4] = conj(B.offd[i][0])  * C.offd[i][3]  + B.diag[i][1]  * C.offd[i][4]
              + B.offd[i][2]  * C.offd[i][5]  + B.offd[i][4]  * C.diag[i][3]
              + B.offd[i][7]  * conj(C.offd[i][9]) + B.offd[i][11]  * conj(C.offd[i][13]);   

    A.offd[i][5]  = conj(B.offd[i][1])  * C.offd[i][3]  + conj(B.offd[i][2])  * C.offd[i][4]
              + B.diag[i][2]  * C.offd[i][5]  + B.offd[i][5]  * C.diag[i][3]
              + B.offd[i][8] * conj(C.offd[i][9]) + B.offd[i][12] * conj(C.offd[i][13]);   

    A.offd[i][6] = B.diag[i][0]  * C.offd[i][6]  + B.offd[i][0]  * C.offd[i][7]
              + B.offd[i][1]  * C.offd[i][8] + B.offd[i][3]  * C.offd[i][9]
              + B.offd[i][6]  * C.diag[i][4]  + B.offd[i][10]  * conj(C.offd[i][14]); 

    A.offd[i][7] = conj(B.offd[i][0])  * C.offd[i][6]  + B.diag[i][1]  * C.offd[i][7]
              + B.offd[i][2]  * C.offd[i][8] + B.offd[i][4]  * C.offd[i][9]
              + B.offd[i][7]  * C.diag[i][4]  + B.offd[i][11]  * conj(C.offd[i][14]);   

    A.offd[i][8] = conj(B.offd[i][1])  * C.offd[i][6]  + conj(B.offd[i][2])  * C.offd[i][7]
              + B.diag[i][2]  * C.offd[i][8] + B.offd[i][5]  * C.offd[i][9]
              + B.offd[i][8] * C.diag[i][4]  + B.offd[i][12] * conj(C.offd[i][14]);   

    A.offd[i][9] = conj(B.offd[i][3])  * C.offd[i][6]  + conj(B.offd[i][4])  * C.offd[i][7]
              + conj(B.offd[i][5])  * C.offd[i][8] + B.diag[i][3]  * C.offd[i][9]
              + B.offd[i][9] * C.diag[i][4]  + B.offd[i][13] * conj(C.offd[i][14]);   

    A.offd[i][10] = B.diag[i][0]  * C.offd[i][10]  + B.offd[i][0]  * C.offd[i][11]
              + B.offd[i][1]  * C.offd[i][12] + B.offd[i][3]  * C.offd[i][13]
              + B.offd[i][6]  * C.offd[i][14] + B.offd[i][10]  * C.diag[i][5];      

    A.offd[i][11] = conj(B.offd[i][0])  * C.offd[i][10]  + B.diag[i][1]  * C.offd[i][11]
              + B.offd[i][2]  * C.offd[i][12] + B.offd[i][4]  * C.offd[i][13]
              + B.offd[i][7]  * C.offd[i][14] + B.offd[i][11]  * C.diag[i][5];      


    A.offd[i][12] = conj(B.offd[i][1])  * C.offd[i][10]  + conj(B.offd[i][2])  * C.offd[i][11]
              + B.diag[i][2]  * C.offd[i][12] + B.offd[i][5]  * C.offd[i][13]
              + B.offd[i][8] * C.offd[i][14] + B.offd[i][12] * C.diag[i][5];    

    A.offd[i][13] = conj(B.offd[i][3])  * C.offd[i][10]  + conj(B.offd[i][4])  * C.offd[i][11]
              + conj(B.offd[i][5])  * C.offd[i][12] + B.diag[i][3]  * C.offd[i][13]
              + B.offd[i][9] * C.offd[i][14] + B.offd[i][13] * C.diag[i][5];   

    A.offd[i][14] = conj(B.offd[i][6])  * C.offd[i][10]  + conj(B.offd[i][7])  * C.offd[i][11]
              + conj(B.offd[i][8]) * C.offd[i][12] + conj(B.offd[i][9]) * C.offd[i][13]
              + B.diag[i][4]  * C.offd[i][14] + B.offd[i][14] * C.diag[i][5];    
  
  } 

  //! Computes Tr(A)
  template<typename R>
  RScalar<R> tri6_trace(const PrimitiveClovTriang<R>& A, int comp){
    RScalar<R> tr = 0.0;
    for (size_t i = 0; i < 2*Nc; i++)
        tr += A.diag[comp][i];
    return tr;
  }

  // Computes Tr(A*B) using hermiticity of the matrices
  template<typename R>
  RScalar<R> tri6_trace_mul(const PrimitiveClovTriang<R>& A, const PrimitiveClovTriang<R>& B, int comp){
    RScalar<R> trd = 0.0, tro = 0.0;
    for (size_t i = 0; i < 2*Nc; i++)
        trd += A.diag[comp][i] * B.diag[comp][i];

    for (size_t i = 0; i < 2*Nc*Nc-Nc; i++)
        tro += RScalar<R>(A.offd[comp][i].real() * B.offd[comp][i].real() + A.offd[comp][i].imag() * B.offd[comp][i].imag());
    return trd + tro + tro;
  }   

  template<typename R>
  void exponentiate(PrimitiveClovTriang<R>& A, int sign){
    const int Niter = 17;
    PrimitiveClovTriang<R> A2, A3, tmp;
    for(size_t comp = 0; comp < 2; comp++){

      tri6_mult(A2, A, A, comp);
      tri6_mult(A3, A, A2, comp);
      
      std::array<RScalar<R>, 7> tr;
      tr[0] = 0;
      tr[1] = 0;
      tr[2] = tri6_trace(A2, comp);
      tr[3] = tri6_trace(A3, comp);
      tr[4] = tri6_trace_mul(A2, A2, comp);
      tr[5] = tri6_trace_mul(A2, A3, comp);
      tr[6] = tri6_trace_mul(A3, A3, comp);

      std::array<RScalar<R>, 5> p;
      p[0] = RScalar<R>(-1.0/6.0)*tr[6] + RScalar<R>(1.0/18.0)*tr[3]*tr[3] 
           - RScalar<R>(1.0/48.0)*tr[2]*tr[2]*tr[2] + RScalar<R>(1.0/8.0)*tr[4]*tr[2];
      p[1] = RScalar<R>(-1.0/5.0)*tr[5] + RScalar<R>(1.0/6.0) *tr[2]*tr[3];
      p[2] = RScalar<R>(-1.0/4.0)*tr[4] + RScalar<R>(1.0/8.0) *tr[2]*tr[2];
      p[3] = RScalar<R>(-1.0/3.0)*tr[3];
      p[4] = RScalar<R>(-1.0/2.0)*tr[2];

      std::array<RScalar<R>, Niter+1> c;
      c[0] = 1;
      for (size_t i = 1; i < Niter+1; i++)
          c[i] = c[i - 1] / RScalar<R>(i);
      if (sign == 1)
          for (size_t i = 1; i < Niter+1; i += 2)
              c[i] = -c[i];

      std::array<RScalar<R>, 6> q;    
      if (Niter > 5){ 
          int ic = Niter - 6;
          for (size_t i = 0; i < 6; i++)
              q[i] = c[ic+1+i];
          
          while (ic >= 0) {
              RScalar<R> q5 = q[5];
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

      // I*q0 + A*q1 + A^2*q2 + A^3*q3
      for (size_t i = 0; i < 2*Nc; i++)
          A.diag[comp][i] = q[0] + A.diag[comp][i]*q[1] + A2.diag[comp][i]*q[2] + A3.diag[comp][i]*q[3];   
      for (size_t i = 0; i < 2*Nc*Nc-Nc; i++)
          A.offd[comp][i] = A.offd[comp][i]*q[1] + A2.offd[comp][i]*q[2] + A3.offd[comp][i]*q[3];

      // A^4*q4
      tri6_mult(tmp, A2, A2, comp);
      for (size_t i = 0; i < 2*Nc; i++)
          A.diag[comp][i] += tmp.diag[comp][i]*q[4];   
      for (size_t i = 0; i < 2*Nc*Nc-Nc; i++)
          A.offd[comp][i] += tmp.offd[comp][i]*q[4];

      // A^5*q5
      tri6_mult(tmp, A3, A2, comp);
      for (size_t i = 0; i < 2*Nc; i++)
          A.diag[comp][i] += tmp.diag[comp][i]*q[5];   
      for (size_t i = 0; i < 2*Nc*Nc-Nc; i++)
          A.offd[comp][i] += tmp.offd[comp][i]*q[5];
    }
  }   
}

#endif