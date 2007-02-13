// $Id: unprec_s_cprec_t_wilson_linop_w.cc,v 1.1 2007-02-13 22:25:04 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_s_cprec_t_wilson_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 


  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	  gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   * \param aniso   anisotropy struct   	       (Read)
   */
  void UnprecSCprecTWilsonLinOp::create(Handle< FermState<T,P,Q> > fs_,
					const Real& Mass_,
					const AnisoParam_t& anisoParam_)
  {
    START_CODE();

    // Check we are in 4D
    if ( Nd != 4 ) { 
      QDPIO::cout << "This class (UnprecSCprecTWilsonLinOp) only works in 4D" << endl;
      QDP_abort(1);
    }
    
    // Check Aniso Direction -- has to be 3.
    if ( anisoParam_.t_dir != 3 ) { 
      QDPIO::cout << "This class (UnprecSCprecTWilsonLinOp) is hardwired for t_dir=3"<< endl;
      QDP_abort(1);
    }

    // Check that the time direction is local:
    const multi1d<int>& s_size =  QDP::Layout::subgridLattSize();  // Local Lattice
    const multi1d<int>& t_size =  QDP::Layout::lattSize(); // Total Latt Size
    if( t_size[3] != s_size[3] ) { 
      QDPIO::cout << "This class (UnprecSCprecTWilsonLinOp) needs time to be local" << endl;
      QDP_abort(1);
    }

    // Store gauge state etc
    fs = fs_;
    aniso = anisoParam_;
    Mass = Mass_;

    // Work out aniso factors
    Real ff = where(aniso.anisoP, aniso.nu / aniso.xi_0, Real(1));
    fact = 1 + (Nd-1)*ff + Mass;
    invfact = Real(1)/fact;
    u = fs->getLinks();

    // Incorporate into links
    if (aniso.anisoP)
    {
      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != aniso.t_dir)
	  u[mu] *= ff;
      }
    }

    // These will be useful for indexing things
    x_index =0;
    y_index =1;
    z_index =2;
    t_index =3;

    // These will be useful for controlling loops.
    Nx = QDP::Layout::subgridLattSize()[x_index];
    Ny = QDP::Layout::subgridLattSize()[y_index];
    Nz = QDP::Layout::subgridLattSize()[z_index];
    Nt = QDP::Layout::subgridLattSize()[t_index];

    // Resize lookup table: for spatial site z,y,x we will get back an array of length t giving the site indices
    // for the t sites, with that spatial index.
    //
    tsite.resize(Nz,Ny,Nx);
   
    // Loop over all spatial sites
    for(int z=0; z < Nz; z++) { 
      for(int y=0; y < Ny; y++) { 
	for(int x=0; x < Nx; x++) { 
	  tsite(z,y,x).resize(Nt);

	  // Assemble the site indices for all timeslices of this coordinate
	  multi1d<int> tsite_indices(Nt);
	  multi1d<int> coords(Nd);
	  coords[0]=x;
	  coords[1]=y;
	  coords[2]=z;
	  for(int t=0; t < Nt; t++) {
	    coords[3] = t;
	    tsite(z,y,x)[t] = QDP::Layout::linearSiteIndex(coords);
	  }
	}
      }
    }

    // P and P_mat dag are needed for the Woodbury 
    // (P_mat for inverting T, P_mat_dag for inverting T_dag - NB P_mat_dag != (P_mat)^dagger
    P_mat.resize(Nz,Ny,Nx);
    P_mat_dag.resize(Nz,Ny,Nx);
    Q_mat_inv.resize(Nz,Ny,Nx);
    Q_mat_dag_inv.resize(Nz,Ny,Nx);

    for(int z=0; z < Nz; z++) { 
      for(int y=0; y < Ny; y++) { 
	for(int x=0; x < Nx; x++) { 

	  P_mat(z,y,x).resize(Nt);
	  P_mat_dag(z,y,x).resize(Nt);

	  Real minvfact = Real(-1)*invfact;

	  // Compute P by backsubsitution - See Balint's notes eq 38-42
	  P_mat(z,y,x)[Nt-1] = u[t_index].elem( tsite(z,y,x)[Nt-1] );
	  P_mat(z,y,x)[Nt-1] *= minvfact.elem();


	  for(int t=Nt-2; t >=0; t--) { 
	    P_mat(z,y,x)[t] =  u[t_index].elem( tsite(z,y,x)[t]  )  * P_mat(z,y,x)[t+1];
	    P_mat(z,y,x)[t] *= invfact.elem();
	  }

	  // Compute the dagger. Opposite order (forward sub) similara to eq 38-42
	  // in notes
	  P_mat_dag(z,y,x)[0] = adj( u[t_index].elem( tsite(z,y,x)[Nt-1]) );
	  P_mat_dag(z,y,x)[0] *= minvfact.elem();
	  for(int t=1; t < Nt; t++) { 
	    P_mat_dag(z,y,x)[t] = adj( u[t_index].elem(tsite(z,y,x)[t-1]) )*P_mat_dag(z,y,x)[t-1] ;
	    P_mat_dag(z,y,x)[t] *= invfact.elem();
	  }


	  // Compute Q = (1 + W^\dag P)^{-1}, W = [1, 0, 0, 0...] => (1 + P_{0})^{-1} eq: 43
	  // NB: This is not necessarily SU(3) now, so we can't just take the dagger to get the inverse
	  Real one=Real(1);
	  CMat tmp = one.elem() + P_mat(z,y,x)[0];
	  invert3by3( Q_mat_inv(z,y,x), tmp );

	  // Compute Q_dag = 1 + W^\dag P^\dag, W = [ 0, ..., 1 ]  = > Q_dag = P^\dag [Nt-1]
	  // Similar to eq 43
	  // NB: This is not necessarily SU(3) now, so we can't just take the dagger to get the inverse
	  tmp = one.elem() + P_mat_dag(z,y,x)[Nt-1];
	  invert3by3( Q_mat_dag_inv(z,y,x), tmp);
	}
      }
    }

    
    END_CODE();
  }

  //! Apply (C_L)^{-1}
  void UnprecSCprecTWilsonLinOp::invCLeftLinOp(T& chi, 
					       const T& psi, 
					       enum PlusMinus isign) const

  {
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T;

    // This is just a way of doing P+ psi - decompose and reconstruct.
    // It may be more straightforward to just write a projector ... -- later
    tmp_plus  = spinProjectDir3Plus(psi);
    chi = spinReconstructDir3Plus(tmp_plus);

    // Here I project - apply TOp to only the halfector - then reconstruct
    tmp_minus = spinProjectDir3Minus(psi);
    TOp(tmp_T, tmp_minus, isign);
    chi += spinReconstructDir3Minus(tmp_T);
    chi *= Real(0.5);
  } 

  //! Apply (C_R)^{-1}
  void UnprecSCprecTWilsonLinOp::invCRightLinOp(T& chi, 
						const T& psi, 
						enum PlusMinus isign) const 
  {
    // Right op is P- + P+ T^{\dagger} so I need the other isign from what I am given
    enum PlusMinus other_sign = (isign == PLUS ? MINUS : PLUS) ;
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T;

    // This does the P- modulo a factor of 1/2
    // Rather than spooling through 2 half vectors I could just do a 
    // ProjectDir3Minus rather than going through half vectirs
    tmp_minus  = spinProjectDir3Minus(psi);
    chi = spinReconstructDir3Minus(tmp_minus);

    // This does the P+ T^\dagger modulo a factor of 1/2
    tmp_plus = spinProjectDir3Plus(psi);

    // Apply T
    TOp(tmp_T, tmp_plus, other_sign);
    chi += spinReconstructDir3Plus(tmp_T);

    chi *= Real(0.5); //The overall factor of 1/2
  } 


  //! Apply C_L
  void UnprecSCprecTWilsonLinOp::cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const
  {
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T;

    tmp_plus  = spinProjectDir3Plus(psi);
    chi = spinReconstructDir3Plus(tmp_plus);

    tmp_minus = spinProjectDir3Minus(psi);
    invTOp(tmp_T, tmp_minus, isign);
    chi += spinReconstructDir3Minus(tmp_T);
    chi *= Real(0.5);
  }

  //! Apply C_R
  void UnprecSCprecTWilsonLinOp::cRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const
  {
    enum PlusMinus other_sign = isign == PLUS ? MINUS : PLUS ;
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T;

    tmp_minus  = spinProjectDir3Minus(psi);
    chi = spinReconstructDir3Minus(tmp_minus);

    tmp_plus = spinProjectDir3Plus(psi);
    invTOp(tmp_T, tmp_plus, other_sign);
    chi += spinReconstructDir3Plus(tmp_T);
    chi *= Real(0.5);
  }

  //! Apply the the space block onto a source vector
  // Apply -1/2 Dslash_s
  // Scabbed this from dslash operator. Could consider
  // using the level 3 one here.
  void UnprecSCprecTWilsonLinOp::spaceLinOp(T& chi, 
					    const T& psi, 
					    enum PlusMinus isign) const {


    switch(isign) {
    case PLUS:  
      {
	LatticeHalfFermion tmp;
	LatticeHalfFermion tmp2;
	
	tmp  = spinProjectDir0Minus(psi);
	tmp2 = shift(tmp, FORWARD, 0);
	chi = spinReconstructDir0Minus(u[0]*tmp2);
	
	tmp   = spinProjectDir1Minus(psi);
	tmp2  = shift(tmp, FORWARD, 1);
	chi  += spinReconstructDir1Minus(u[1]*tmp2);
	
	tmp   = spinProjectDir2Minus(psi);
	tmp2  = shift(tmp, FORWARD, 2);
	chi  += spinReconstructDir2Minus(u[2]*tmp2);
	
	tmp  = adj(u[0])*spinProjectDir0Plus(psi);
	tmp2 = shift(tmp, BACKWARD, 0);
	chi += spinReconstructDir0Plus(tmp2);
	
	tmp  = adj(u[1])*spinProjectDir1Plus(psi);
	tmp2 = shift(tmp, BACKWARD, 1);
	chi += spinReconstructDir1Plus(tmp2);
	
	tmp  = adj(u[2])*spinProjectDir2Plus(psi);
	tmp2 = shift(tmp, BACKWARD, 2);
	chi += spinReconstructDir2Plus(tmp2);
      }
      break;
      
    case MINUS:
      {
	LatticeHalfFermion tmp;
	LatticeHalfFermion tmp2;
	
	tmp  = spinProjectDir0Plus(psi);
	tmp2 = shift(tmp, FORWARD, 0);
	chi  = spinReconstructDir0Plus(u[0]*tmp2);
	
	tmp  = spinProjectDir1Plus(psi);
	tmp2 = shift(tmp, FORWARD, 1);
	chi += spinReconstructDir1Plus(u[1]*tmp2);

	tmp = spinProjectDir2Plus(psi);
	tmp2 = shift(tmp, FORWARD, 2);
	chi += spinReconstructDir2Plus(u[2]*tmp2);
	
	tmp  = adj(u[0])*spinProjectDir0Minus(psi);
	tmp2 = shift(tmp, BACKWARD, 0);
	chi += spinReconstructDir0Minus(tmp2);
	
	tmp  = adj(u[1])*spinProjectDir1Minus(psi);
	tmp2 = shift(tmp, BACKWARD, 1);
	chi += spinReconstructDir1Minus(tmp2);
	
	tmp  = adj(u[2])*spinProjectDir2Minus(psi);
	tmp2 = shift(tmp, BACKWARD, 2);
	chi += spinReconstructDir2Minus(tmp2);
	
      }

      break;
  
    default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
    }

    chi *= Real(-0.5);
  }


  //! Apply T or T^\dagger
  inline
  void UnprecSCprecTWilsonLinOp::TOp(LatticeHalfFermion& chi, 
				     const LatticeHalfFermion& psi, 
				     enum PlusMinus isign) const 
  {


    
    // Loop over all spatial sites
    for(int z=0; z < Nz; z++) { 
      for(int y=0; y < Ny; y++) { 
	for(int x=0; x < Nx; x++) { 
	  
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
	      chi.elem(tsite(z,y,x)[Nt-1]) = fact.elem()*psi.elem(tsite(z,y,x)[Nt-1]);
	      chi.elem(tsite(z,y,x)[Nt-1]) -= u[t_index].elem(tsite(z,y,x)[Nt-1]) * psi.elem(tsite(z,y,x)[0]);
	      
	      // Rest of the rows
	      for(int t=Nt-2; t >= 0; t--) { 
		chi.elem(tsite(z,y,x)[t])  = fact.elem()*psi.elem(tsite(z,y,x)[t]);
		chi.elem(tsite(z,y,x)[t]) -= u[t_index].elem(tsite(z,y,x)[t]) * psi.elem(tsite(z,y,x)[t+1]);
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
	      // where t_i = tsite(z,y,x)( i )
	      //
	      //  This need a loop over spatial sites
	      
	      // Last row
	      
	      chi.elem(tsite(z,y,x)[0]) = fact.elem()*psi.elem(tsite(z,y,x)[0]);
	      chi.elem(tsite(z,y,x)[0]) -= adj( u[t_index].elem(tsite(z,y,x)[Nt-1]) ) * psi.elem(tsite(z,y,x)[Nt-1]);
	      
	      for(int t=1; t < Nt; t++) { 
		
		chi.elem(tsite(z,y,x)[t]) = fact.elem() *psi.elem(tsite(z,y,x)[t]);
		chi.elem(tsite(z,y,x)[t])  -= adj( u[t_index].elem(tsite(z,y,x)[t-1]) ) * psi.elem(tsite(z,y,x)[t-1]);
	      }
	    }
	    break;
	  default:
	    QDPIO::cout << "Unknown Sign " << endl;
	    QDP_abort(1);
	  }
	  
	} // x loop
      } // y loop
    } // z loop
  }


  // Apply T^{-1} or ( T^\dagger )^{-1} using Sherman Morrison Woodbury
  inline
  void UnprecSCprecTWilsonLinOp::invTOp(LatticeHalfFermion& chi, 
					const LatticeHalfFermion& psi, 
					enum PlusMinus isign) const {

    for(int z=0; z < Nz; z++) { 
      for(int y=0; y < Ny; y++) { 
	for(int x=0; x < Nx; x++) {

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
	      chi.elem( tsite(z,y,x)[Nt-1] )= invfact.elem()*psi.elem( tsite(z,y,x)[Nt-1] );
	      
	      
	      // Rest of the rows
	      for(int t=Nt-2; t >= 0; t--) { 
		chi.elem( tsite(z,y,x)[t] )  = psi.elem( tsite(z,y,x)[t] );
		chi.elem( tsite(z,y,x)[t] ) += u[t_index].elem( tsite(z,y,x)[t] ) * chi.elem( tsite(z,y,x)[t+1] );
		chi.elem( tsite(z,y,x)[t] ) *= invfact.elem();
		
	      }
	  
	      // Now  get  ( 1 - P Q W^\dag ) z  ( SMW Formula )
	      // P and Q precomputed in constructor as P_mat and Q_mat_inv

	      // Compute z - P Q W^\dag z
	      // W^\dag = (1, 0, 0, ...) 
	      // W^\dag z = z_0
	      // Let x = Q z_0 = Q W z
	      HVec_site x_site = Q_mat_inv(z,y,x) * chi.elem( tsite(z,y,x)[0] );

	      // Now ( 1 - P Q W z ) = z - P x 
	      for(int t=0; t < Nt; t++) { 
		chi.elem( tsite(z,y,x)[t] ) -= P_mat(z,y,x)[t] * x_site;
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
	      chi.elem( tsite(z,y,x)[0] )= invfact.elem()*psi.elem( tsite(z,y,x)[0] );
	      
	      for(int t=1; t < Nt; t++) { 
		chi.elem( tsite(z,y,x)[t] ) = psi.elem( tsite(z,y,x)[t] );
		chi.elem( tsite(z,y,x)[t] ) += adj( u[t_index].elem(tsite(z,y,x)[t-1] )  ) * chi.elem( tsite(z,y,x)[t-1] );
		chi.elem( tsite(z,y,x)[t] ) *= invfact.elem();
	      }
	      // Compute z - P^dag Q^dag W^dag z
	      // W^\dag = ( 0, ...., 1) => \W^\dag z = z_Nt-1

	      HVec_site x_site = Q_mat_dag_inv(z,y,x) * chi.elem( tsite(z,y,x)[Nt-1] );

	      // Now ( 1 - P Q W z) = z - P x 
	      for(int t=0; t < Nt; t++) { 
		chi.elem( tsite(z,y,x)[t] ) -= P_mat_dag(z,y,x)[t]*x_site;
	      }

	      // Done! Not as painful as I thought
	    }
	    break;
	  default:
	    QDPIO::cout << "Unknown Sign " << endl;
	    QDP_abort(1);
	  }       

	} // x
      } // y
    } // z
    
  } // end of function

  inline
  void UnprecSCprecTWilsonLinOp::invert3by3( CMat& M_inv, const CMat& M ) const
  {
    START_CODE();
    
    int Nvec = 3;
    
    // Copy b and M so that we can change elements of the copies
    // during the pivoting and the LU decomposition. We can save
    // memory by just destroying M, but we never really expect
    // it to be really big: At most 20x20 I should imagine.

    PScalar<PColorMatrix<RComplex<REAL>, 3> > M_tmp;
    for(int i=0; i < 3; i++) { 
      for(int j=0; j < 3; j++) { 
	M_tmp.elem().elem(i,j).real() = 0;
	M_tmp.elem().elem(i,j).imag() = 0;
      }
      M_tmp.elem().elem(i,i).real() = 1;
    }
    
    PScalar<PColorMatrix<RComplex<REAL>, 3> > M_local = M;
    
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


} // End Namespace Chroma
