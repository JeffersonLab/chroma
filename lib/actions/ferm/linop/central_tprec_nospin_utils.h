// -*- C++ -*-
// $Id: central_tprec_nospin_utils.h,v 1.6 2008-05-07 01:12:18 bjoo Exp $
/*! \file
 *  \brief Support for time preconditioning
 */

#ifndef CENTRAL_TPREC_NOSPIN_UTILS_H
#define CENTRAL_TPREC_NOSPIN_UTILS_H
#include "qdp_config.h"

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4
#include "chromabase.h"

namespace Chroma 
{

  //! Support for time preconditioning
  /*!
   * \ingroup linop
   */
  namespace CentralTPrecNoSpinUtils 
  {
    typedef PScalar<PColorMatrix<RComplex<REAL>, 3> > CMat;              // Useful type: ColorMat with no Outer<>
    typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, 2 > HVec_site;   // Useful type: Half Vec with no Outer<>

    /*! \ingroup linop */
    inline 
    void TOp(LatticeHalfFermion& chi, 
	     const LatticeHalfFermion& psi,
	     const multi1d<LatticeColorMatrix>& u,
	     const multi2d<int>& tsite,
	     const Real& fact,
	     enum PlusMinus isign, 
	     const bool schroedingerTP=false) {
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
	    // If boundaries are not Schroedinger add the lower left piece
	    if( ! schroedingerTP ) {
	      chi.elem(tsite(site,Nt-1)) -= u[t_index].elem(tsite(site,Nt-1)) * psi.elem( tsite(site,0) );
	    }
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

	    // If not Schroedinger, then subtrac off upper right contribution
	    if( !schroedingerTP ) {
	      chi.elem(tsite(site,0)) -= adj( u[t_index].elem(tsite(site,Nt-1)) ) * psi.elem(tsite(site,Nt-1));
	    }

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
		enum PlusMinus isign,
		const int  t_max,
		const bool schroedingerTP=false)
    {
#ifndef QDP_IS_QDPJIT
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
	
	    // NB: For Schroedinger boundaries we don't need 
	    // To do the Woodbury piece since the matrix is strictly
	    // bidiagonal (no wraparound piece that needs Woodbury correction)
	    if( !schroedingerTP ) { 
	      // Now  get  ( 1 - P Q W^\dag ) z  ( SMW Formula )
	      // P and Q precomputed in constructor as P_mat and Q_mat_inv
	
	      // Compute z - P Q W^\dag z
	      // W^\dag = (1, 0, 0, ...) 
	      // W^\dag z = z_0
	      // Let x = Q z_0 = Q W z
	      HVec_site x_site = Q_mat_inv(site) * chi.elem( tsite(site,0) );
	      
	      // Now ( 1 - P Q W z ) = z - P x

	      // OK Basically we can save flops here.
	      // We only need to add on Px if P not numerically zero
	      int loop_end;
	      if( Nt-1-t_max <= 0 ) {
		loop_end=-1;
	      }
	      else { 
		loop_end=Nt-1-t_max;
	      }
		
	      for(int t=Nt-1; t > loop_end ; t--) { 
		  chi.elem( tsite(site,t) ) -= P_mat(site,t) * x_site;
	      }
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

	    // NB: For Schroedinger boundaries we don't need 
	    // To do the Woodbury piece since the matrix is strictly
	    // bidiagonal (no wraparound piece that needs Woodbury correction)
	    if( !schroedingerTP ) { 
	      // Compute z - P^dag Q^dag W^dag z
	      // W^\dag = ( 0, ...., 1) => \W^\dag z = z_Nt-1
	      
	      HVec_site x_site = Q_mat_dag_inv(site) * chi.elem( tsite(site,Nt-1) );
	      
	      // Now ( 1 - P Q W z) = z - P x 
	      int loop_end;
	      if ( t_max >= Nt-1 ) {
		loop_end=Nt;
	      }
	      else { 
		loop_end=t_max;
	      }

	      for(int t=0; t < loop_end; t++) { 
		  chi.elem( tsite(site,t) ) -= P_mat_dag(site,t)*x_site;
	      }
	    }
	    // Done! Not as painful as I thought
	  }
	  break;
	default:
	  QDPIO::cout << "Unknown Sign " << endl;
	  QDP_abort(1);
	}       
      }
#endif
    }
		

    /*! \ingroup linop */
    inline 
    void invert3by3( CMat& M_inv, const CMat& M )
    { 
#ifndef QDP_IS_QDPJIT
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
#endif
    }

    inline
    Double logDet(const CMat& M) {
#ifndef QDP_IS_QDPJIT
     // Possibly this is the Dumb way but it is only a small matrix
     // This is to be done by the matrix of cofactors:
     //
     //   [  M_00  M_01  M_02  ]                  
     //   [  M_10  M_11  M_12  ] 
     //   [  M_20  M_21  M_22  ] 
     // 
     // I express by the top row so M_00 det( A ) - M_01 det(B) + M_02 det(C)
     //
     //       [ M_11 M_12 ]        [ M_10  M_12 ]        [ M_10  M_11 ]
     //   A = [           ]   B =  [            ]    C = [            ]
     //       { M_21 M_22 ]        [ M_20  M_22 ]        [ M_20  M_21 ]
     //
     //  and the 2 by 2 determinant is (eg for A) det(A)= M_11 M_22 - M_21 M_12
     //
     //  Since our Matrix M has to be of the form ( 1 + P_0 )^\dagger (1 + P_0) 
     //  we KNOW that the determinant has to be real, so just take real part and return.

     //PScalar< PScalar< PScalar< RComplex<REAL> > > > ret_val;
     Complex ret_val;

     ret_val.elem().elem().elem() = M.elem().elem(0,0)*( M.elem().elem(1,1)*M.elem().elem(2,2) - M.elem().elem(2,1)*M.elem().elem(1,2) );
     ret_val.elem().elem().elem() -=  M.elem().elem(0,1)*( M.elem().elem(1,0)*M.elem().elem(2,2) - M.elem().elem(2,0)*M.elem().elem(1,2) );
     ret_val.elem().elem().elem() +=  M.elem().elem(0,2)*( M.elem().elem(1,0)*M.elem().elem(2,1) - M.elem().elem(2,0)*M.elem().elem(1,1) );
     
     return Double(log(real(ret_val)));
#endif
   }
    

    inline
    void derivLogDet(multi1d<LatticeColorMatrix>& F, 
		     const multi1d<LatticeColorMatrix>& U,
		     const multi1d<CMat>& Q,
		     const multi2d<int>& tsites,
		     const int t_dir,
		     const Real NdPlusM,
		     const bool schroedingerTP=false)
    {
#ifndef QDP_IS_QDPJIT


      F.resize(Nd);

      // Initial setup. Zero out all the non-time directions
      for(int mu=0; mu < Nd; mu++) {
	if( mu != t_dir ) {
	  F[mu] = zero;
	}
      }

      if( ! schroedingerTP ) { 
	
	LatticeColorMatrix T1;
	Real invmass = Real(1)/ NdPlusM;
	Real factor;
	
	// Work out factor
	int Nspace=tsites.size2();
	int Nt = tsites.size1();
	
	factor=Real(2);
	for(int t=0; t < Nt; t++) {
	  factor*=invmass;
	}
	
	for(int site=0; site < Nspace; site++) {
	  
	  // Compute force for all timeslices T
	  for(int t=0; t < Nt; t++) { 
	    
	    // Initialize temporary
	    CMat F_tmp;
	    
	    // This creates a unit matrix
	    for(int r=0; r < 3; r++) { 
	      for(int c=0; c < 3; c++) { 
		F_tmp.elem().elem(r,c).real() = 0;
		F_tmp.elem().elem(r,c).imag() = 0;
	      }
	    }
	    
	    for(int r=0; r < 3; r++) { 
	      F_tmp.elem().elem(r,r).real() = 1;
	      F_tmp.elem().elem(r,r).imag() = 0;
	    }
	    
	    CMat F_tmp2;
	    
	    // Do the U-s from t+1 to Nt unless you already are Nt-1
	    for(int j=t+1; j < Nt; j++)  {
	      F_tmp2 = F_tmp;
	      int s = tsites(site,j);
	      F_tmp = F_tmp2 * U[t_dir].elem(s);
	    }
	    
	    // Insert the Q
	    F_tmp2 = F_tmp;
	    F_tmp = F_tmp2 * Q[site];
	    
	    // Do the U-s from zero up to t-1
	    for(int j=0; j < t; j++) { 
	      F_tmp2 = F_tmp;
	      int s = tsites(site,j);
	      F_tmp = F_tmp2*U[t_dir].elem(s);
	    }
	    // Copy force to temporary
	    T1.elem(tsites(site,t)) = F_tmp;
	    
	  }
	}
	
	// Create the force
	F[t_dir] = factor*T1;
      }
      else { 
	// Schroedinger BC-s. Matrix is upper/lower bidiagonal
	// and the determinant is a constant independent of the
	// gauge fields.-In this case we ignore the constant
	// and generate no_force
	F[t_dir] = zero;
      }
#endif
    }


    // For use by the checkerboarded version...
    inline
    void derivLogDet(multi1d<LatticeColorMatrix>& F, 
		     const multi1d<LatticeColorMatrix>& U,
		     const multi2d<CMat>& Q,
		     const multi3d<int>& tsites,
		     const int t_dir,
		     const Real NdPlusM,
		     const bool schroedingerTP=false)
    {
#ifndef QDP_IS_QDPJIT


      F.resize(Nd);

      // Initial setup. Zero out all the non-time directions
      for(int mu=0; mu < Nd; mu++) {
	if( mu != t_dir ) {
	  F[mu] = zero;
	}
      }

      if( ! schroedingerTP ) { 
	
	LatticeColorMatrix T1;
	Real invmass = Real(1)/ NdPlusM;
	Real factor;
	
	// Work out factor
	int Nspace=tsites.size2();
	int Nt = tsites.size1();
	
	factor=Real(2);
	for(int t=0; t < Nt; t++) {
	  factor*=invmass;
	}
	
	for(int cb3=0; cb3 < 2; cb3++) { 
	  for(int site=0; site < Nspace; site++) {
	    
	    // Compute force for all timeslices T
	    for(int t=0; t < Nt; t++) { 
	      
	      // Initialize temporary
	      CMat F_tmp;
	      
	      // This creates a unit matrix
	      for(int r=0; r < 3; r++) { 
		for(int c=0; c < 3; c++) { 
		  F_tmp.elem().elem(r,c).real() = 0;
		  F_tmp.elem().elem(r,c).imag() = 0;
		}
	      }
	      
	      for(int r=0; r < 3; r++) { 
		F_tmp.elem().elem(r,r).real() = 1;
		F_tmp.elem().elem(r,r).imag() = 0;
	      }
	      
	      CMat F_tmp2;
	      
	      // Do the U-s from t+1 to Nt unless you already are Nt-1
	      for(int j=t+1; j < Nt; j++)  {
		F_tmp2 = F_tmp;
		int s = tsites(cb3,site,j);
		F_tmp = F_tmp2 * U[t_dir].elem(s);
	      }
	      
	      // Insert the Q
	      F_tmp2 = F_tmp;
	      F_tmp = F_tmp2 * Q[cb3][site];
	      
	      // Do the U-s from zero up to t-1
	      for(int j=0; j < t; j++) { 
		F_tmp2 = F_tmp;
		int s = tsites(cb3,site,j);
		F_tmp = F_tmp2*U[t_dir].elem(s);
	      }
	      // Copy force to temporary
	      T1.elem(tsites(cb3,site,t)) = F_tmp;
	      
	    }
	  }
	}
	// Create the force
	F[t_dir] = factor*T1;
      }
      else { 
	// Schroedinger BC-s. Matrix is upper/lower bidiagonal
	// and the determinant is a constant independent of the
	// gauge fields.-In this case we ignore the constant
	// and generate no_force
	F[t_dir] = zero;
      }
#endif
    }
  } // Namespace 
  
} // Namespace chroma

#endif
#endif
#endif

#endif
