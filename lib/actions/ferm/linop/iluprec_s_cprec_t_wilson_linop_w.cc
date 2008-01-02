// $Id: iluprec_s_cprec_t_wilson_linop_w.cc,v 1.4 2008-01-02 17:05:54 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "chromabase.h"
#include "actions/ferm/linop/iluprec_s_cprec_t_wilson_linop_w.h"
#include "actions/ferm/linop/central_tprec_nospin_utils.h"

using namespace QDP::Hints;

namespace Chroma 
{ 


  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	  gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   * \param aniso   anisotropy struct   	       (Read)
   */
  void ILUPrecSCprecTWilsonLinOp::create(Handle< FermState<T,P,Q> > fs_,
					const Real& Mass_,
					const AnisoParam_t& anisoParam_)
  {
    START_CODE();

    // Check we are in 4D
    if ( Nd != 4 ) { 
      QDPIO::cout << "This class (ILUPrecSCprecTWilsonLinOp) only works in 4D" << endl;
      QDP_abort(1);
    }
    
    // Check Aniso Direction -- has to be 3.
    if ( anisoParam_.t_dir != 3 ) { 
      QDPIO::cout << "This class (ILUPrecSCprecTWilsonLinOp) is hardwired for t_dir=3"<< endl;
      QDP_abort(1);
    }

    // Check that the time direction is local:
    const multi1d<int>& s_size =  QDP::Layout::subgridLattSize();  // Local Lattice
    const multi1d<int>& t_size =  QDP::Layout::lattSize(); // Total Latt Size
    if( t_size[3] != s_size[3] ) { 
      QDPIO::cout << "This class (ILUPrecSCprecTWilsonLinOp) needs time to be local" << endl;
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

    // These will be useful for controlling loops.
    int t_index=3;
    int Nx = QDP::Layout::subgridLattSize()[0];
    int Ny = QDP::Layout::subgridLattSize()[1];
    int Nz = QDP::Layout::subgridLattSize()[2];
    int Nt = QDP::Layout::subgridLattSize()[3];

    // Resize lookup table: for spatial site z,y,x we will get back an array of length t giving the site indices
    // for the t sites, with that spatial index.
    //
    int Nspaceby2 = Nx*Ny*Nz/2;  // This always works because Nx has to be even...

    tsite.resize(2, Nspaceby2, Nt);


    int ssite_even=0;
    int ssite_odd=0;

    int color;
    int site;

    // Loop over all spatial sites
    for(int z=0; z < Nz; z++) { 
      for(int y=0; y < Ny; y++) { 
	for(int x=0; x < Nx; x++) { 
	  // Assemble the site indices for all timeslices of this coordinate

	  multi1d<int> coords(Nd);
	  coords[0]=x;
	  coords[1]=y;
	  coords[2]=z;
	  color = (x + y + z) & 1;

	  if ( color == 0 ) {     
	    for(int t=0; t < Nt; t++) {
	      coords[3] = t;
	      tsite(0,ssite_even, t) = QDP::Layout::linearSiteIndex(coords);
	    }
	    ssite_even++;

	  }
	  else { 
	    for(int t=0; t < Nt; t++) {
	      coords[3] = t;
	      tsite(1,ssite_odd, t) = QDP::Layout::linearSiteIndex(coords);
	    }
	    ssite_odd++;
	  }
	}
      }
    }

    // P and P_mat dag are needed for the Woodbury 
    // (P_mat for inverting T, P_mat_dag for inverting T_dag - NB P_mat_dag != (P_mat)^dagger
    P_mat.resize(2, Nspaceby2, Nt);
    P_mat_dag.resize(2, Nspaceby2, Nt);

    Q_mat_inv.resize(2,Nspaceby2);
    Q_mat_dag_inv.resize(2,Nspaceby2);

    logDetTSq = zero;
    for(int cb3=0; cb3 < 2; cb3++) { 
      for(int site=0; site < Nspaceby2; site++) { 
	
	
	Real minvfact = Real(-1)*invfact;
	
	// Compute P by backsubsitution - See Balint's notes eq 38-42
	P_mat(cb3, site, Nt-1) = u[t_index].elem( tsite(cb3,site,Nt-1) );
	P_mat(cb3, site, Nt-1) *= minvfact.elem();
      
      
	for(int t=Nt-2; t >=0; t--) { 
	  P_mat(cb3, site, t) =  u[t_index].elem( tsite(cb3, site,t)  )  * P_mat(cb3, site,t+1);
	  P_mat(cb3, site, t) *= invfact.elem();
	}
	
	// Compute the dagger. Opposite order (forward sub) similara to eq 38-42
	// in notes
	P_mat_dag(cb3, site, 0) = adj( u[t_index].elem( tsite(cb3, site, Nt-1) ) );
	P_mat_dag(cb3, site, 0) *= minvfact.elem();

	for(int t=1; t < Nt; t++) { 
	  P_mat_dag(cb3, site,t) = adj( u[t_index].elem(  tsite(cb3, site, t-1) ) )*P_mat_dag(cb3,site, t-1) ;
	  P_mat_dag(cb3, site,t) *= invfact.elem();
	}
      
      
	// Compute Q = (1 + W^\dag P)^{-1}, W = [1, 0, 0, 0...] => (1 + P_{0})^{-1} eq: 43
	// NB: This is not necessarily SU(3) now, so we can't just take the dagger to get the inverse
	Real one=Real(1);
	CMat one_plus_P0 = one.elem() + P_mat(cb3, site,0);
	CentralTPrecNoSpinUtils::invert3by3( Q_mat_inv(cb3, site), one_plus_P0 );
	
	// Compute Q_dag = 1 + W^\dag P^\dag, W = [ 0, ..., 1 ]  = > Q_dag = P^\dag [Nt-1]
	// Similar to eq 43
	// NB: This is not necessarily SU(3) now, so we can't just take the dagger to get the inverse
	CMat one_plus_P0_dag = one.elem() + P_mat_dag(cb3, site,Nt-1);
	CentralTPrecNoSpinUtils::invert3by3( Q_mat_dag_inv(cb3, site), one_plus_P0_dag);

	CMat prod = one_plus_P0_dag*one_plus_P0;
	logDetTSq += CentralTPrecNoSpinUtils::logDet(prod);
      }
    }

    Dw3D.create(fs_, anisoParam_);    
    END_CODE();
  }

  //! Apply (C_L)^{-1}
  void ILUPrecSCprecTWilsonLinOp::TPlusOp(T& chi, 
					  const T& psi, 
					  enum PlusMinus isign,
					  int cb3d) const

  {
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T;

    // This is just a way of doing P+ psi - decompose and reconstruct.
    // It may be more straightforward to just write a projector ... -- later
    tmp_plus[ rb3[ cb3d ] ]  = spinProjectDir3Plus(psi);
    chi[ rb3[ cb3d ] ] = spinReconstructDir3Plus(tmp_plus);

    // Here I project - apply TOp to only the halfector - then reco nstruct
    tmp_minus[ rb3[ cb3d ] ] = spinProjectDir3Minus(psi);


    // Use shared routine to apply T or T^\dagger.
    // Pass in u, and tsite for the subset
    CentralTPrecNoSpinUtils::TOp(tmp_T, 
				 tmp_minus, 
				 u,
				 tsite[cb3d],
				 fact,
				 isign);

    chi[rb3[ cb3d ]]  += spinReconstructDir3Minus(tmp_T);
    chi[rb3[ cb3d ]]  *= Real(0.5);
  } 

  //! Apply (C_R)^{-1}
  void ILUPrecSCprecTWilsonLinOp::TMinusOp(T& chi, 
					   const T& psi, 
					   enum PlusMinus isign,
					   int cb3d) const 
  {
    // Right op is P- + P+ T^{\dagger} so I need the other isign from what I am given
    enum PlusMinus other_sign = (isign == PLUS ? MINUS : PLUS) ;
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T;

    // This does the P- modulo a factor of 1/2
    // Rather than spooling through 2 half vectors I could just do a 
    // ProjectDir3Minus rather than going through half vectirs
    tmp_minus[rb3[cb3d]]  = spinProjectDir3Minus(psi);
    chi[rb3[cb3d]] = spinReconstructDir3Minus(tmp_minus);

    // This does the P+ T^\dagger modulo a factor of 1/2
    tmp_plus[rb3[cb3d]] = spinProjectDir3Plus(psi);


    // Use shared routine to apply T or T^\dagger.
    // Pass in u, and tsite for the subset
    CentralTPrecNoSpinUtils::TOp(tmp_T, 
				 tmp_plus, 
				 u,
				 tsite[cb3d],
				 fact,
				 other_sign);

    chi[rb3[cb3d]] += spinReconstructDir3Plus(tmp_T);

    chi[rb3[cb3d]] *= Real(0.5); //The overall factor of 1/2
  } 


  //! Apply C_L
  void ILUPrecSCprecTWilsonLinOp::invTPlusOp(T& chi, const T& psi, enum PlusMinus isign,
					     int cb3d) const
  {
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T;

    tmp_plus[rb3[cb3d]]  = spinProjectDir3Plus(psi);
    chi[rb3[cb3d]] = spinReconstructDir3Plus(tmp_plus);

    tmp_minus[rb3[cb3d]] = spinProjectDir3Minus(psi);

    // Use shared routine to apply (T)^{-1} or  (T^\dagger)^{-1}.
    // Pass in u, and tsite, P, P^\dag, Q and Q^\dag for the subset.
    CentralTPrecNoSpinUtils::invTOp( tmp_T,
				     tmp_minus, 
				     u,
				     tsite[cb3d],
				     P_mat[cb3d],
				     P_mat_dag[cb3d],
				     Q_mat_inv[cb3d],
				     Q_mat_dag_inv[cb3d],
				     invfact,
				     isign);


    chi[rb3[cb3d]] += spinReconstructDir3Minus(tmp_T);
    chi[rb3[cb3d]] *= Real(0.5);
  }

  //! Apply C_R
  void ILUPrecSCprecTWilsonLinOp::invTMinusOp(T& chi, const T& psi, enum PlusMinus isign,
					      int cb3d) const
  {
    enum PlusMinus other_sign = isign == PLUS ? MINUS : PLUS ;
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T;

    tmp_minus[rb3[cb3d]]  = spinProjectDir3Minus(psi);
    chi[rb3[cb3d]] = spinReconstructDir3Minus(tmp_minus);

    tmp_plus[rb3[cb3d]] = spinProjectDir3Plus(psi);

    // Use shared routine to apply (T)^{-1} or  (T^\dagger)^{-1}.
    // Pass in u, and tsite, P, P^\dag, Q and Q^\dag for the subset.
    CentralTPrecNoSpinUtils::invTOp( tmp_T,
				     tmp_plus, 
				     u,
				     tsite[cb3d],
				     P_mat[cb3d],
				     P_mat_dag[cb3d],
				     Q_mat_inv[cb3d],
				     Q_mat_dag_inv[cb3d],
				     invfact,
				     other_sign);


    chi[rb3[cb3d]] += spinReconstructDir3Plus(tmp_T);
    chi[rb3[cb3d]] *= Real(0.5);
  }

  void  ILUPrecSCprecTWilsonLinOp::derivCRight(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const 
  {
    T T_1, T_2, T_3, T_4;
    switch(isign) {
    case PLUS: 
      {
	// Simple Case: -T^o^\dagger_2 \dot(Dslash^oe_S) T^e_1
	
	// [ T^e_1 ]   [ (T-)^{-1}e  0 ] [ Y_e ]   [ (T-)^{-1}_e Y_e ]
	// [       ] = [               ] [     ] = [                 ]
	// { T^o_1 ]   [  0          1 ] [ Y_o ]   [   Y_o           ]
	invTMinusOp(T_1, Y,  PLUS, 0);
	T_1[ rb3[1] ] = Y;

	// [ T^e_2 ]   [ 1           0   ] [ X_e ]  [   X_e             ]
	// [       ] = [                 ] [     ] =[                   ]  
	// [ T^o_2 ]   [ 0  (T-)^{-dag}o ] [ X_o ]  [ (T-)^{-dag}_o X_o ]
	invTMinusOp(T_2, X, MINUS, 1);
	T_2[ rb3[0] ] = X;

	derivSpaceOp(ds_u, T_2, T_1, PLUS, 1);
	for(int mu=0; mu < Nd; mu++) { 
	  ds_u[mu] *= Real(-1);
	}
      }
      break;
    case MINUS:
      {
	// [ T^e_1 ]   [ 1           0   ] [ Y_e ]  [   Y_e             ]
	// [       ] = [                 ] [     ] =[                   ]  
	// [ T^o_1 ]   [ 0  (T-)^{-dag}o ] [ Y_o ]  [ (T-)^{-dag}_o Y_o ]
	invTMinusOp(T_1, Y, MINUS, 1);
	T_1[ rb3[0] ] = Y;


	// [ T^e_2 ]   [ (T-)^{-1}e  0 ] [ X_e ]   [ (T-)^{-1}_e X_e ]
	// [       ] = [               ] [     ] = [                 ]
	// { T^o_2 ]   [  0          1 ] [ X_o ]   [   X_o           ]
	invTMinusOp(T_2, X, PLUS, 0);
	T_2[ rb3[1] ] = X;

	// [ T^e_3 ]   [ 1   - D^{oe}^\dag_s ][ T^e_1 ]
	// [       ] = [                     ][       ]
	// [ T^o_3 ]   [ 0          1        ][ T^o_1 ]
	// We only need the even part tho, and put Y_o in the bottom
	T T_tmp;
	spaceLinOp(T_tmp, T_1, MINUS, 0);
	T_3[rb3[0]] = T_1 - T_tmp; // Mhalf from the (-1/2)dslash
	T_3[rb3[1]] = Y;


	// [ T^e_4 ]   [ 1           0  ][ T^e_2 ]
	// [       ] = [                ][       ]
	// [ T^o_4 ]   [ -D^{eo}_s   1  ][ T^o_2 ]
	// We only need the odd part tho, and put X_o in the top
	spaceLinOp(T_tmp, T_2, PLUS, 1);
	T_4[rb3[1]] = T_2 - T_tmp;       // Mhalf from the (-1/2)dslash
	T_4[rb3[0]] = X;


	derivInvTMinusOp(ds_u, T_4, T_3, MINUS);

	P ds_tmp;

	derivSpaceOp(ds_tmp, T_2, T_1, MINUS, 0);
	for(int mu=0; mu < Nd; mu++) {
	  ds_u[mu] -= ds_tmp[mu];           
	}

      }
      break;
    default:
      QDPIO::cerr << "Bad Case. Should never get here" << endl;
      QDP_abort(1);
    }
  }

  void  ILUPrecSCprecTWilsonLinOp::derivCLeft(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const
  {
    T T_1, T_2, T_3, T_4;

    switch(isign) {
    case PLUS: 
      {
	// [ T^e_1 ]   [ 1      0        ] [ Y_e ]  [   Y_e             ]
	// [       ] = [                 ] [     ] =[                   ]  
	// [ T^o_1 ]   [ 0   (T+)^{-1}_o ] [ Y_o ]  [ (T+)^{-1}_o Y_o   ]
	invTPlusOp(T_1, Y, PLUS, 1);
	T_1[ rb3[0] ] = Y;


	// [ T^e_2 ]   [ (T+)^{-dag}e  0 ] [ X_e ]   [ (T+)^{-dag}_e X_e ]
	// [       ] = [                 ] [     ] = [                   ]
	// { T^o_2 ]   [  0            1 ] [ X_o ]   [   X_o             ]
	invTPlusOp(T_2, X, MINUS, 0);
	T_2[ rb3[1] ] = X;

	// [ T^e_3 ]   [ 1   - D^{eo}_s ][ T^e_1 ]
	// [       ] = [                ][       ]
	// [ T^o_3 ]   [ 0          1   ][ T^o_1 ]
	// We only need the even part tho, and put Y_o in the bottom
	T T_tmp;
	spaceLinOp(T_tmp, T_1, PLUS, 0);
	T_3[rb3[0]] = T_1 - T_tmp;  // mhalf is from (-1/2)Dslash
	T_3[rb3[1]] = Y;


	// [ T^e_4 ]   [ 1                0  ][ T^e_2 ]
	// [       ] = [                     ][       ]
	// [ T^o_4 ]   [ -D^{oe \dag}_s   1  ][ T^o_2 ]
	// We only need the odd part tho, and put X_o in the top
	spaceLinOp(T_tmp, T_2, MINUS, 1);
	T_4[rb3[1]] = T_2 - T_tmp; // mhalf is from (-1/2) Dslash
	T_4[rb3[0]] = X;

	P ds_tmp;
	derivSpaceOp(ds_tmp, T_2, T_1, PLUS, 0);
	derivInvTPlusOp(ds_u, T_4, T_3, PLUS);
	for(int mu=0; mu < Nd; mu++) { 
	  ds_u[mu] -= ds_tmp[mu];
	}
      }
      break;

    case MINUS:
      {
	// Simple Case: -T^o^\dagger_2 \dot(Dslash^oe_S) T^e_1
	
	// [ T^e_1 ]   [ (T+)^{-\dag}e  0 ] [ Y_e ]   [ (T+)^{-\dag}_e Y_e ]
	// [       ] = [                  ] [     ] = [                 ]
	// { T^o_1 ]   [  0             1 ] [ Y_o ]   [   Y_o           ]
	invTPlusOp(T_1, Y,  MINUS, 0);
	T_1[ rb3[1] ] = Y;

	// [ T^e_2 ]   [ 1           0   ] [ X_e ]  [   X_e             ]
	// [       ] = [                 ] [     ] =[                   ]  
	// [ T^o_2 ]   [ 0  (T+)^{-1}o   ] [ X_o ]  [ (T+)^{-1}_o X_o ]
	invTPlusOp(T_2, X, PLUS, 1);
	T_2[ rb3[0] ] = X;

	derivSpaceOp(ds_u, T_2, T_1, MINUS, 1);
	for(int mu=0; mu < Nd; mu++) { 
	  ds_u[mu] *= Real(-1);  // (Combinme factor of -1 and -(1/2)
	}
      }
      break;
    default:
      QDPIO::cerr << "Bad Case. Should never get here" << endl;
      QDP_abort(1);
    }

  }
  
  void ILUPrecSCprecTWilsonLinOp::derivInvTPlusOp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const 
  
  {
  
    // Same code as C_L in unprec case
    ds_u.resize(Nd);
    for(int mu=0; mu < 3; mu++) { 
      ds_u[mu] = zero;
    }
    
    if (isign == PLUS) { 
      
      LatticeHalfFermion tmp1, tmp2;
      T T1, T2; 
      
      tmp1=spinProjectDir3Minus(Y);
      for(int cb3d=0; cb3d < 2; cb3d++) { 
	CentralTPrecNoSpinUtils::invTOp( tmp2,
					 tmp1, 
					 u,
					 tsite[cb3d],
					 P_mat[cb3d],
					 P_mat_dag[cb3d],
					 Q_mat_inv[cb3d],
					 Q_mat_dag_inv[cb3d],
					 invfact,
					 PLUS);
      }

      T1  = spinReconstructDir3Minus(tmp2);
           
      tmp1=spinProjectDir3Minus(X);

      for(int cb3d=0; cb3d < 2; cb3d++) { 
      CentralTPrecNoSpinUtils::invTOp( tmp2,
				       tmp1, 
				       u,
				       tsite[cb3d],
				       P_mat[cb3d],
				       P_mat_dag[cb3d],
				       Q_mat_inv[cb3d],
				       Q_mat_dag_inv[cb3d],
				       invfact,
				       MINUS);
      }

      // Two factors of 0.5 from the projectors.
      // Most cost efficient to apply them together to the half vector...
      tmp2 *= Real(0.25);

      T2  = spinReconstructDir3Minus(tmp2);
     
      
      LatticeFermion T3 = shift(T1, FORWARD, 3);
      
      // A minus from the fact that its dT^{-1}
      // A minus from the fact that the shift has a -
      // Makes  a + 
      ds_u[3] = traceSpin(outerProduct(T3, T2));
      
    }
    else {
      ds_u[3] = zero;
    }

  } 

  void ILUPrecSCprecTWilsonLinOp::derivInvTMinusOp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const 
  {
    // Same code as derivC_R 
    ds_u.resize(Nd);
    for(int mu=0; mu < 3; mu++) { 
      ds_u[mu] = zero;
    }
    
    if (isign == MINUS) { 
      
      LatticeHalfFermion tmp1, tmp2;
      T T1, T2; 
      
      tmp1=spinProjectDir3Plus(Y);
      for(int cb3d=0; cb3d < 2; cb3d++) { 
      CentralTPrecNoSpinUtils::invTOp( tmp2,
				       tmp1, 
				       u,
				       tsite[cb3d],
				       P_mat[cb3d],
				       P_mat_dag[cb3d],
				       Q_mat_inv[cb3d],
				       Q_mat_dag_inv[cb3d],
				       invfact,
				       PLUS);
      }      
      T1  = spinReconstructDir3Plus(tmp2);
      
      tmp1=spinProjectDir3Plus(X);
      for(int cb3d=0; cb3d < 2; cb3d++) { 

      CentralTPrecNoSpinUtils::invTOp( tmp2,
				       tmp1, 
				       u,
				       tsite[cb3d],
				       P_mat[cb3d],
				       P_mat_dag[cb3d],
				       Q_mat_inv[cb3d],
				       Q_mat_dag_inv[cb3d],
				       invfact,
				       MINUS);
      }
      // Two factors of 0.5 from the projectors.
      // Most cost efficient to apply them together to the half vector...
      tmp2 *= Real(0.25);
      
      T2  = spinReconstructDir3Plus(tmp2);
      
      LatticeFermion T3 = shift(T1, FORWARD, 3);
      
      // A minus from the fact that its dT^{-1}
      // A minus from the fact that the shift has a -
      // Makes  a + 
      ds_u[3] = traceSpin(outerProduct(T3, T2));
      
    }
    else { 
      ds_u[3]= zero;
    }


  }

} // End Namespace Chroma

#endif
#endif
#endif
