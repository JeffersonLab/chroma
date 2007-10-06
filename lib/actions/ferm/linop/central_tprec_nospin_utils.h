// -*- C++ -*-
// $Id: central_tprec_nospin_utils.h,v 1.2 2007-10-06 15:33:09 edwards Exp $
/*! \file
 *  \brief Support for time preconditioning
 */

#ifndef CENTRAL_TPREC_NOSPIN_UTILS_H
#define CENTRAL_TPREC_NOSPIN_UTILS_H

#include "chromabase.h"

namespace Chroma 
{

  //! Support for time preconditioning
  /*!
   * \ingroup linop
   */
  namespace CentralTPrecNoSpinUtils 
  {
    typedef PScalar<PColorMatrix<RComplex<REAL>, Nc> > CMat;              // Useful type: ColorMat with no Outer<>
    typedef PSpinVector<PColorVector<RComplex<REAL>, Nc>, (Ns>>1) > HVec_site;   // Useful type: Half Vec with no Outer<>

    /*! \ingroup linop */
    inline 
    void TOp(LatticeHalfFermion& chi, 
	     const LatticeHalfFermion& psi,
	     const multi1d<LatticeColorMatrix>& u,
	     const multi2d<int>& tsite,
	     const Real& fact,
	     enum PlusMinus isign) {
      const int t_index = 3; // Fixed
      int Nt = tsite.size1();
      
      for(int site=0; site < tsite.size2(); site++) { 
	
	// tsite[t] now holds the index for the elem() so that we can cycle
	// through the time
	
	// Now we need to apply different matrices 
	switch(isign) { 
	case PLUS:
	  {
	    // Now we need to apply the matrix
	    //
	    // [ chi( t_0 ) ]  [ fact    -U(x,t_0)  0 ...                   ] [ psi( t_0 )    ]
	    // [ chi( t_1 ) ]  [  0         fact   -U(x,t_1) 0 ...          ] [ psi( t_1 )    ]
	    // [  ....      ]= [  ...        0      fact     -U(x, t_Nt-2)  ] [ psi( ... )    ]
	    // [ chi(t_Nt-1)]  [ -U(x, t_Nt-1)    0  ....         fact      ] [ psi( t_Nt-1 ) ]
	    //
	    // where t_i = tsite( i )
	    //
	    //  This need a loop over spatial sites
	    
	    // Last row
	    chi.elem(tsite(site,Nt-1)) = fact.elem()*psi.elem( tsite(site,Nt-1) );
	    chi.elem(tsite(site,Nt-1)) -= u[t_index].elem(tsite(site,Nt-1)) * psi.elem( tsite(site,0) );
	    
	    // Rest of the rows
	    for(int t=Nt-2; t >= 0; t--) { 
	      chi.elem( tsite(site,t) )  = fact.elem()*psi.elem(tsite(site,t));
	      chi.elem( tsite(site,t) ) -= u[t_index].elem(tsite(site,t)) * psi.elem(tsite(site,t+1) );
	    }
	  }
	  break;
	case MINUS:
	  {
	    // Now we need to apply the matrix
	    //
	    // [ chi( t_0 ) ]  [ fact               0             0 ...  -U^{+}(x,t_Nt-1) ] [ psi( t_0 )    ]
	    // [ chi( t_1 ) ]  [ -U^{+}(x,t_0)     fact           0  ...          0       ] [ psi( t_1 )    ]
	    // [  ....      ]= [   0           -U^{+}(x,t_1)     fact             0       ] [ psi( ... )    ]
	    // [ chi(t_Nt-1)]  [   0                        -U^{+}(x,t_Nt-2)    fact      ] [ psi( t_Nt-1 ) ]
	    //
	    // where t_i = tsite(site, i )
	    //
	    //  This need a loop over spatial sites
	    
	    // Last row
	    
	    chi.elem(tsite(site,0)) = fact.elem()*psi.elem(tsite(site,0));
	    chi.elem(tsite(site,0)) -= adj( u[t_index].elem(tsite(site,Nt-1)) ) * psi.elem(tsite(site,Nt-1));
	    
	    for(int t=1; t < Nt; t++) { 
	      
	      chi.elem(tsite(site,t)) = fact.elem() *psi.elem(tsite(site,t));
	      chi.elem(tsite(site,t))  -= adj( u[t_index].elem( tsite(site,t-1) ) ) * psi.elem(tsite(site,t-1) );
	    }
	  }
	  break;
	default:
	  QDPIO::cout << "Unknown Sign " << endl;
	  QDP_abort(1);
	}

      }     
    }



    /*! \ingroup linop */
    inline 
    void invTOp(LatticeHalfFermion& chi, 
		const LatticeHalfFermion& psi,
		const multi1d<LatticeColorMatrix>& u,
		const multi2d<int>& tsite,
		const multi2d<CMat>& P_mat,
		const multi2d<CMat>& P_mat_dag,
		const multi1d<CMat>& Q_mat_inv,
		const multi1d<CMat>& Q_mat_dag_inv,
		const Real& invfact,
		enum PlusMinus isign)
    {
      const int t_index=3;
      int Nt=tsite.size1();

      for(int site=0; site < tsite.size2(); site++) { 
	
	// I can index the psi and u directly, but I can keep the temporary in just a site's worth.
	// Better for cacheing maybe.
	//	  multi1d<HVec_site> z_strip(Nt);
	
	// First apply (T_0)^{-1}  -- easy to do with backward / forward (for dagger) subsitution
	//
	// z = T_{0}^{-1} psi 
	switch(isign) { 
	case PLUS:
	  {
	    // Now we need to solve
	    //
	    // [ fact    -U(x,t_0)  0 ...                   ] [ chi( t_0 )    ]   [ psi( t_0 )    ]
	    // [  0         fact   -U(x,t_1) 0 ...          ] [ chi( t_1 )    ]   [ psi( t_1 )    ]
	    // [  ...        0      fact     -U(x, t_Nt-2)  ] [ chi( ... )    ] = [ psi( ... )    ]
	    // [  0              0  ....         fact       ] [ chi( t_Nt-1 ) ]   [ psi( t_Nt-1 ) ]
	    //
	    // where t_i = tsite(z,y,x)( i )
	    //
	    // by backsubstitution
	    
	    // Last row
	    chi.elem( tsite(site,Nt-1) )= invfact.elem()*psi.elem( tsite(site,Nt-1) );
	    
	    
	    // Rest of the rows
	    for(int t=Nt-2; t >= 0; t--) { 
	      chi.elem( tsite(site,t) )  = psi.elem( tsite(site,t) );
	      chi.elem( tsite(site,t) ) += u[t_index].elem( tsite(site,t) ) * chi.elem( tsite(site,t+1) );
	      chi.elem( tsite(site,t) ) *= invfact.elem();
	      
	    }
	    
	    // Now  get  ( 1 - P Q W^\dag ) z  ( SMW Formula )
	    // P and Q precomputed in constructor as P_mat and Q_mat_inv
	    
	    // Compute z - P Q W^\dag z
	    // W^\dag = (1, 0, 0, ...) 
	    // W^\dag z = z_0
	    // Let x = Q z_0 = Q W z
	    HVec_site x_site = Q_mat_inv(site) * chi.elem( tsite(site,0) );
	    
	    // Now ( 1 - P Q W z ) = z - P x 
	    for(int t=0; t < Nt; t++) { 
	      chi.elem( tsite(site,t) ) -= P_mat(site,t) * x_site;
	    }
	    
	    // Done! Not as painful as I thought
	  }
	  break;
	case MINUS:
	  {
	    // Now we need to apply the matrix
	    //
	    // [ fact               0             0 ...           0       ] [ chi( t_0 )    ]     [ psi(t_0)    ]
	    // [ -U^{+}(x,t_0)     fact           0  ...          0       ] [ chi( t_1 )    ]     [ psi(t_1)    ]
	    // [   0           -U^{+}(x,t_1)     fact             0       ] [ chi( ... )    ]  =  [ psi(...)    ]
	    // [   0                        -U^{+}(x,t_Nt-2)    fact      ] [ chi( t_Nt-1 ) ]     [ psi(t_Nt-1) ]
	    //
	    // where t_i = tsite(z,y,x)( i )
	    //
	    //  This need a loop over spatial sites
	    
	    // Last row
	    chi.elem( tsite(site,0) )= invfact.elem()*psi.elem( tsite(site,0) );
	    
	    for(int t=1; t < Nt; t++) { 
	      chi.elem( tsite(site,t) ) = psi.elem( tsite(site,t) );
	      chi.elem( tsite(site,t) ) += adj( u[t_index].elem(tsite(site,t-1) )  ) * chi.elem( tsite(site,t-1) );
	      chi.elem( tsite(site,t) ) *= invfact.elem();
	    }
	    // Compute z - P^dag Q^dag W^dag z
	    // W^\dag = ( 0, ...., 1) => \W^\dag z = z_Nt-1
	    
	    HVec_site x_site = Q_mat_dag_inv(site) * chi.elem( tsite(site,Nt-1) );
	    
	    // Now ( 1 - P Q W z) = z - P x 
	    for(int t=0; t < Nt; t++) { 
	      chi.elem( tsite(site,t) ) -= P_mat_dag(site,t)*x_site;
	    }
	    
	    // Done! Not as painful as I thought
	  }
	  break;
	default:
	  QDPIO::cout << "Unknown Sign " << endl;
	  QDP_abort(1);
	}       
      }
    }
		

    /*! \ingroup linop */
    inline 
    void invert3by3( CMat& M_inv, const CMat& M )
    { 
      START_CODE();
    
      int Nvec = 3;
      
      // Copy b and M so that we can change elements of the copies
      // during the pivoting and the LU decomposition. We can save
      // memory by just destroying M, but we never really expect
      // it to be really big: At most 20x20 I should imagine.
      
      PScalar<PColorMatrix<RComplex<REAL>, Nc> > M_tmp;
      for(int i=0; i < 3; i++) { 
	for(int j=0; j < 3; j++) { 
	  M_tmp.elem().elem(i,j).real() = 0;
	  M_tmp.elem().elem(i,j).imag() = 0;
	}
	M_tmp.elem().elem(i,i).real() = 1;
      }
      
      PScalar<PColorMatrix<RComplex<REAL>, Nc> > M_local = M;
      
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
	REAL maxnorm = M_local.elem().elem(i,i).real()*M_local.elem().elem(i,i).real()
	  + M_local.elem().elem(i,i).imag()*M_local.elem().elem(i,i).imag();
	
	int maxrow = i;
	
	// Compare norms with other elements in column j for row i+1.N
	for(int row=i+1; row < Nvec; row++) {
	  REAL normcheck =  M_local.elem().elem(row,i).real()*M_local.elem().elem(row,i).real()
	    + M_local.elem().elem(row,i).imag()*M_local.elem().elem(row,i).imag();
	  
	  if ( toBool( normcheck > maxnorm ) ) {
	    // Norm of M_local(j,i) is bigger, store it as the maximum
	    // and store its index
	    maxnorm = normcheck;
	    maxrow = row;
	  }
	}
	
	// If the element with maximum norm is not in row i, swap
	// its row with row i
	if( maxrow != i ) {
	  Complex tmp;
	  
	  // Swap rows i and maxindex
	  for(int j=0; j < Nvec; j++ ) {
	    
	    tmp.elem().elem().elem() = M_local.elem().elem(i, j);
	    M_local.elem().elem(i,j) = M_local.elem().elem(maxrow, j);
	    M_local.elem().elem(maxrow, j) = tmp.elem().elem().elem();
	    
	    // Swap elems of Minx
	    tmp.elem().elem().elem() = M_tmp.elem().elem(i,j);
	    M_tmp.elem().elem(i,j) = M_tmp.elem().elem(maxrow,j);
	    M_tmp.elem().elem(maxrow,j) = tmp.elem().elem().elem();
	    
	  }
	  
	}
	
	// --------------------------------------------------------
	// End of pivoting code
	// --------------------------------------------------------
	
	
	// --------------------------------------------------------
	// Work out elements of L & U in place in M_local for row i
	// --------------------------------------------------------
	for(int j=0; j < i; j++) { 
	  
	  Complex sum_LU(0);
	  
	  for(int k = 0; k < j; k++) {
	    sum_LU.elem().elem().elem() += M_local.elem().elem(i,k)*M_local.elem().elem(k,j);
	  }
	  
	  M_local.elem().elem(i,j) -= sum_LU.elem().elem().elem();
	  M_local.elem().elem(i,j) /= M_local.elem().elem(j,j);
	}
	
	for(int j=i; j < Nvec; j++) { 
	  Complex sum_LU(0);
	  for(int k = 0; k < i; k++) { 
	    sum_LU.elem().elem().elem() += M_local.elem().elem(i,k)*M_local.elem().elem(k,j);
	  }
	  M_local.elem().elem(i,j) -= sum_LU.elem().elem().elem();
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
      multi1d< Complex > y(Nvec);
      
      for(int k=0; k < Nvec; k++) { 
	
	y[0].elem().elem().elem() = M_tmp.elem().elem(0,k);
	for(int i=1; i < Nvec; i++) { 
	  y[i].elem().elem().elem() = M_tmp.elem().elem(i,k);
	  for(int j=0; j < i; j++) { 
	    y[i].elem().elem().elem() -= M_local.elem().elem(i,j)*y[j].elem().elem().elem();
	  }
	}
	
	// Solve U a = y by back substitution
	M_inv.elem().elem(Nvec-1,k) = y[Nvec-1].elem().elem().elem() / M_local.elem().elem(Nvec-1, Nvec-1);
	
	for(int i = Nvec-2; i >= 0; i--) { 
	  Complex tmpcmpx = y[i];
	  for(int j=i+1; j < Nvec; j++) { 
	    tmpcmpx.elem().elem().elem() -= M_local.elem().elem(i,j)*M_inv.elem().elem(j,k);
	  }
	  M_inv.elem().elem(i,k) = tmpcmpx.elem().elem().elem() /M_local.elem().elem(i,i);
	}
	
      }
    }


  } // Namespace 
  
} // Namespace chroma

#endif
