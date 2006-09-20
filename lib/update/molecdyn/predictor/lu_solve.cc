#include "chromabase.h"
#include "update/molecdyn/predictor/lu_solve.h"

namespace Chroma 
{
  // Solve M a = b by LU decomposition with partial pivoting
  //
  void LUSolve( multi1d<DComplex>& a, 
		const multi2d<DComplex>& M, 
		const multi1d<DComplex>& b ) 
  {
    START_CODE();
    
    int Nvec1 = M.size1();
    int Nvec2 = M.size2();
    
    if( Nvec1 != Nvec2 ) { 
      QDPIO::cerr << "Barf" << endl;
      QDP_abort(1);
    }
    
    int Nvec = Nvec1;
    
    // Copy b and M so that we can change elements of the copies
    // during the pivoting and the LU decomposition. We can save
    // memory by just destroying M, but we never really expect
    // it to be really big: At most 20x20 I should imagine.

    multi1d<DComplex> b_local(b);
    multi2d<DComplex> M_local(M);
    
    // -----------------------------------------------------------------
    // LU Decompose M_local, in place (Crone's algorithm?)
    // It's in Numerical Recipes but also a more understandable
    // description can be found at:
    //          http://csep10.phys.utk.edu/guidry/
    //               phys594/lectures/linear_algebra/lanotes/node3.html
    // 
    // OR look in your favourite Matrix Analysis text
    // -----------------------------------------------------------------
    
    // -------------------------------------------------------------
    // Start LU Decomp. Definition. 0-th row of U is 0-th row of M
    //   and L_{i,i} = 1 for all i
    // 
    // So we start with the 1-th (2nd) row
    // ------------------------------------------------------------
    
    for(int i = 1; i < Nvec; i++) { 
      
      // ------------------------------------------------------------
      // Parital Pivot: Find the row with the largest element in the
      // ith-column and make that the i-th row. This swaps rows.
      // so I don't need to reorder the unknowns, but I do need 
      // to reorder the b_local
      // ------------------------------------------------------------
      Double maxnorm = norm2(M_local(i,i));
      int maxrow = i;
      
      // Compare norms with other elements in column j for row i+1.N
      for(int row=i+1; row < Nvec; row++) {
	if ( toBool( norm2(M_local(row,i)) > maxnorm ) ) {
	  // Norm of M_local(j,i) is bigger, store it as the maximum
	  // and store its index
	  maxnorm = norm2(M_local(row,i));
	  maxrow = row;
	}
      }
      
      // If the element with maximum norm is not in row i, swap
      // its row with row i
      if( maxrow != i ) {
	
	DComplex tmp;
	
	// Swap rows i and maxindex
	for(int j=0; j < Nvec; j++ ) {
	  tmp = M_local(i, j);
	  M_local(i,j) = M_local(maxrow, j);
	  M_local(maxrow, j) = tmp;
	}
	
	// Swap elems of b
	tmp = b_local[i];
	b_local[i] = b_local[maxrow];
	b_local[maxrow] = tmp;
      }
      
      // --------------------------------------------------------
      // End of pivoting code
      // --------------------------------------------------------
      
      
      // --------------------------------------------------------
      // Work out elements of L & U in place in M_local for row i
      // --------------------------------------------------------
      for(int j=0; j < i; j++) { 
	
	DComplex sum_LU = DComplex(0);;
	for(int k = 0; k < j; k++) {
	  sum_LU += M_local(i,k)*M_local(k,j);
	}
	
	M_local(i,j) -= sum_LU;
	M_local(i,j) /= M_local(j,j);
      }
      
      for(int j=i; j < Nvec; j++) { 
	DComplex sum_LU = DComplex(0);
	for(int k = 0; k < i; k++) { 
	  sum_LU += M_local(i,k)*M_local(k,j);
	}
	M_local(i,j) -= sum_LU;
      }

    }
    
    // ----------------------------------------------------
    // LU Decomp finished. M_local now holds the 
    //   U matrix in its diagonal and superdiagonal elements
    //   and the subdiagonal elements of the L matrix in its
    //   subdiagonal. Recall that the Diagonal elements of L 
    //   are chosen to be 1
    // -----------------------------------------------------
    
    // Solve L y = b by forward substitution
    multi1d<DComplex> y(Nvec);
    
    y[0] = b[0];
    for(int i=1; i < Nvec; i++) { 
      y[i] = b_local[i];
      for(int j=0; j < i; j++) { 
	y[i] -= M_local(i,j)*y[j];
      }
    }
    
    // Solve U a = y by back substitution
    a[Nvec-1] = y[Nvec-1] / M_local(Nvec-1, Nvec-1);
    
    for(int i = Nvec-2; i >= 0; i--) { 
      DComplex tmpcmpx = y[i];
      for(int j=i+1; j < Nvec; j++) { 
	tmpcmpx -= M_local(i,j)*a[j];
      }
      a[i] = tmpcmpx/M_local(i,i);
    }
  
    END_CODE();
  }


} // End Namespace
