// $Id: unprec_s_cprec_t_wilson_linop_w.cc,v 1.7 2008-05-07 01:12:19 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */


#include "qdp_config.h"

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4

#include "chromabase.h"
#include "actions/ferm/linop/unprec_s_cprec_t_wilson_linop_w.h"
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
  void UnprecSCprecTWilsonLinOp::create(Handle< FermState<T,P,Q> > fs_,
					const Real& Mass_,
					const AnisoParam_t& anisoParam_)
  {
#ifndef QDP_IS_QDPJIT
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

    u = fs->getLinks();
    if (aniso.anisoP)
    {
      // Work out aniso factors
      Real ff = where(aniso.anisoP, aniso.nu / aniso.xi_0, Real(1));
      fact = 1 + (Nd-1)*ff + Mass;

      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != aniso.t_dir)
	  u[mu] *= ff;
      }
    }
    else { 
      fact = Nd + Mass;
    }

    invfact = Real(1)/fact;
    // These will be useful for controlling loops.
    int t_index=3;
    int Nx = QDP::Layout::subgridLattSize()[0];
    int Ny = QDP::Layout::subgridLattSize()[1];
    int Nz = QDP::Layout::subgridLattSize()[2];
    int Nt = QDP::Layout::subgridLattSize()[3];

    // Resize lookup table: for spatial site z,y,x we will get back an array of length t giving the site indices
    // for the t sites, with that spatial index.
    //
    int Nspace = Nx*Ny*Nz;
    tsite.resize(Nspace, Nt);

    int ssite=0;
    // Loop over all spatial sites
    for(int z=0; z < Nz; z++) { 
      for(int y=0; y < Ny; y++) { 
	for(int x=0; x < Nx; x++) { 
	  // Assemble the site indices for all timeslices of this coordinate
	  multi1d<int> coords(Nd);
	  coords[0]=x;
	  coords[1]=y;
	  coords[2]=z;
	  for(int t=0; t < Nt; t++) {
	    coords[3] = t;
	    tsite(ssite,t) = QDP::Layout::linearSiteIndex(coords);
	  }
	  ssite++;
	}
      }
    }

    // Figure out whether we are Schroedinger in time:
    const FermBC<T,P,Q>& fbc = getFermBC();
    if( fbc.nontrivialP() ) {
      try {
	const SchrFermBC& schr_ref = dynamic_cast<const SchrFermBC&>(fbc);
	schrTP = ( schr_ref.getDir() == tDir() );
      }
      catch(std::bad_cast) { 
	schrTP = false;
      }

    }
    else { 
      // fbc are not nontrivial
      schrTP = false;
    }

    logDetTSq = zero;
    
    if( !schrTP ) {
      // P and P_mat dag are needed for the Woodbury 
      // (P_mat for inverting T, P_mat_dag for inverting T_dag - NB P_mat_dag != (P_mat)^dagger
      P_mat.resize(Nspace, Nt);
      P_mat_dag.resize(Nspace, Nt);
      Q_mat_inv.resize(Nspace);
      Q_mat_dag_inv.resize(Nspace);
      
      
      for(int site=0; site < Nspace; site++) { 
	Real minvfact = Real(-1)*invfact;
	
	// Compute P by backsubsitution - See Balint's notes eq 38-42
	P_mat(site,Nt-1) = u[t_index].elem( tsite(site,Nt-1) );
	P_mat(site,Nt-1) *= minvfact.elem();
	
	
	for(int t=Nt-2; t >=0; t--) { 
	  P_mat(site,t) =  u[t_index].elem( tsite(site,t)  )  * P_mat(site,t+1);
	  P_mat(site,t) *= invfact.elem();
	}
	
	// Compute the dagger. Opposite order (forward sub) similara to eq 38-42
	// in notes
	P_mat_dag(site,0) = adj( u[t_index].elem( tsite(site,Nt-1) ) );
	P_mat_dag(site,0) *= minvfact.elem();
	for(int t=1; t < Nt; t++) { 
	  P_mat_dag(site,t) = adj( u[t_index].elem(  tsite(site,t-1) ) )*P_mat_dag(site, t-1) ;
	  P_mat_dag(site,t) *= invfact.elem();
	}
      
      
	// Compute Q = (1 + W^\dag P)^{-1}, W = [1, 0, 0, 0...] => (1 + P_{0})^{-1} eq: 43
	// NB: This is not necessarily SU(3) now, so we can't just take the dagger to get the inverse
	Real one=Real(1);
	CMat one_plus_P0 = one.elem() + P_mat(site,0);
	CentralTPrecNoSpinUtils::invert3by3( Q_mat_inv(site), one_plus_P0 );
	
	// Compute Q_dag = 1 + W^\dag P^\dag, W = [ 0, ..., 1 ]  = > Q_dag = P^\dag [Nt-1]
	// Similar to eq 43
	// NB: This is not necessarily SU(3) now, so we can't just take the dagger to get the inverse
	CMat one_plus_P0_dag = one.elem() + P_mat_dag(site,Nt-1);
	CentralTPrecNoSpinUtils::invert3by3( Q_mat_dag_inv(site), one_plus_P0_dag);
	
	
	CMat prod = one_plus_P0_dag*one_plus_P0;
	logDetTSq += CentralTPrecNoSpinUtils::logDet(prod);
      }
    }
    else { 
      // Schroedinger Time Case. Matrix(dagger) is strictly upper(lower)
      // bidiagonal - Determinant is therefore a product of the 
      // diagonal elements which is just a factor independent
      // of the gauge fields: (Nd + M)^{Nt}. Since this is a 
      // constant it doesn't need to be simulated and we can
      // happily return logDetTSq=0
      logDetTSq = 0;
    }
    Dw3D.create( fs_, anisoParam_);

    END_CODE();
#endif
  }

  //! Apply (C_L)^{-1}
  // 
  //   C_L^{-1} = P_{+} + P_{-} T
  //
  void UnprecSCprecTWilsonLinOp::invCLeftLinOp(T& chi, 
					       const T& psi, 
					       enum PlusMinus isign) const

  {

    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T=zero;


    // This is just a way of doing P+ psi - decompose and reconstruct.
    // It may be more straightforward to just write a projector ... -- later
    tmp_plus  = spinProjectDir3Plus(psi);
    chi = spinReconstructDir3Plus(tmp_plus);

    // Here I project - apply TOp to only the halfector - then reco nstruct
    tmp_minus = spinProjectDir3Minus(psi);


    // Use shared routine to apply T or T^\dagger.
    // Pass in u, and tsite for the subset
    CentralTPrecNoSpinUtils::TOp(tmp_T, 
				 tmp_minus, 
				 u,
				 tsite,
				 fact,
				 isign,
				 schroedingerTP());

    chi += spinReconstructDir3Minus(tmp_T);
    chi *= Real(0.5);

    getFermBC().modifyF(chi);

  } 

  //! Apply (C_R)^{-1}
  //
  //     C_R^{-1} = P_{-} + P_{+} T^{\dagger} 
  //
  void UnprecSCprecTWilsonLinOp::invCRightLinOp(T& chi, 
						const T& psi, 
						enum PlusMinus isign) const 
  {

    // Right op is P- + P+ T^{\dagger} so I need the other isign from what I am given
    enum PlusMinus other_sign = (isign == PLUS ? MINUS : PLUS) ;
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T=zero;


    // This does the P- modulo a factor of 1/2
    // Rather than spooling through 2 half vectors I could just do a 
    // ProjectDir3Minus rather than going through half vectirs
    tmp_minus  = spinProjectDir3Minus(psi);
    chi = spinReconstructDir3Minus(tmp_minus);

    // This does the P+ T^\dagger modulo a factor of 1/2
    tmp_plus = spinProjectDir3Plus(psi);


    // Use shared routine to apply T or T^\dagger.
    // Pass in u, and tsite for the subset
    CentralTPrecNoSpinUtils::TOp(tmp_T, 
				 tmp_plus, 
				 u,
				 tsite,
				 fact,
				 other_sign,
				 schroedingerTP());

    chi += spinReconstructDir3Plus(tmp_T);

    chi *= Real(0.5); //The overall factor of 1/2
    getFermBC().modifyF(chi);
  } 


  //! Apply C_L
  //
  //    C_L = P_{+} + P_{-} T^{-1} 
  //
  void UnprecSCprecTWilsonLinOp::cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const
  {

    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T=zero;

    //  2 P_{+} = ( 1 + gamma_3 )
    tmp_plus  = spinProjectDir3Plus(psi);
    chi = spinReconstructDir3Plus(tmp_plus);

    // 2 P_{-} = (1 - gamma_3 )
    tmp_minus = spinProjectDir3Minus(psi);

    // Use shared routine to apply (T)^{-1} or  (T^\dagger)^{-1}.
    // Pass in u, and tsite, P, P^\dag, Q and Q^\dag for the subset.
    CentralTPrecNoSpinUtils::invTOp( tmp_T,
				     tmp_minus, 
				     u,
				     tsite,
				     P_mat,
				     P_mat_dag,
				     Q_mat_inv,
				     Q_mat_dag_inv,
				     invfact,
				     isign,
				     getTMax(),
				     schroedingerTP());

    // Reconstruct
    chi += spinReconstructDir3Minus(tmp_T);

    // Overall factor of 2 to turn (1 +/- gamma_3) into projector P_{+/-}
    //  No sense to fold it into the half vector because it 
    //  would also be needed for the half vector in the P_{+} piece
    //  so overall cost is still 1 full vector of multiply
    chi *= Real(0.5);
    getFermBC().modifyF(chi);

  }

  //! Apply C_R
  //
  //    C_L = P_{-} + P_{+} [ T^{\dagger} ]^{-1}
  //
  void UnprecSCprecTWilsonLinOp::cRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const
  {

    enum PlusMinus other_sign = isign == PLUS ? MINUS : PLUS ;
    LatticeHalfFermion tmp_minus;
    LatticeHalfFermion tmp_plus;
    LatticeHalfFermion tmp_T = zero;


    // 2*P_{-} 
    tmp_minus  = spinProjectDir3Minus(psi);
    chi = spinReconstructDir3Minus(tmp_minus);

    // 2*P_{+} 
    tmp_plus = spinProjectDir3Plus(psi);

    // Use shared routine to apply (T)^{-1} or  (T^\dagger)^{-1}.
    // Pass in u, and tsite, P, P^\dag, Q and Q^\dag for the subset.
    CentralTPrecNoSpinUtils::invTOp( tmp_T,
				     tmp_plus, 
				     u,
				     tsite,
				     P_mat,
				     P_mat_dag,
				     Q_mat_inv,
				     Q_mat_dag_inv,
				     invfact,
				     other_sign,
				     getTMax(),
				     schroedingerTP());

    // Reconstruct to full vector
    chi += spinReconstructDir3Plus(tmp_T);

    // Overall factor of 2 to turn (1 +/- gamma_3) into projector P_{+/-}
    //  No sense to fold it into the half vector because it 
    //  would also be needed for the half vector in the P_{+} piece
    //  so overall cost is still 1 full vector of multiply
    chi *= Real(0.5);
    getFermBC().modifyF(chi);

  }

  //! Apply the the space block onto a source vector
  //  Call to 3D Dslash.  - overall factor of -0.5from -(1/2) Dslash
  void 
  UnprecSCprecTWilsonLinOp::spaceLinOp(T& chi, const T& psi, enum PlusMinus isign) const
  {
    Real mhalf=Real(-0.5);
    Dw3D.apply(chi, psi, isign,0);
    Dw3D.apply(chi, psi, isign,1);
    chi *= mhalf;
    getFermBC().modifyF(chi);
  }

  //! Derivative of the derivCLeft
  // 
  //    X^\dagger dC_l Y 
  //    = X^\dagger d [  P_{-} T^{-1} ] Y 
  //    = - X^\dagger  P{-} T^{-1} dT T^{-1} Y 
  //    = - X^\dagger  P{-} P_{-} T^{-1} dT T^{-1} Y  since P_{-} P_{-} = P_{-}
  //    = - X^\dagger  T^{-1} P_{-} dT P_{-} T^{-1} Y since Ps and Ts commute
  //    = - L^\dagger dT R
  // 
  //   with R = P_{-} T^{-1} Y and L = P_{-} T^{-dagger} X
  //
  //   T contains only U terms (not U-dagger) and dT/dU R = -shift(R, FORWARD, t)
  //   T^\dagger contains no U terms so           dT^\dagger / dU = zero
  //
  //   similarly to how we do dslash-es. dT/dU^\dagger terms contribute only to the 
  //   hermitian conjugate terms which are automagically taken care of when we do the
  //   TAProj() later on
  void   
  UnprecSCprecTWilsonLinOp::derivCLeft(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const
  {
    ds_u.resize(Nd);
    for(int mu=0; mu < 3; mu++) { 
      ds_u[mu] = zero;
    }
    
    if (isign == PLUS) { 
      
      LatticeHalfFermion tmp1, tmp2;
      T T1, T2; 
      
      tmp1=spinProjectDir3Minus(Y);
      CentralTPrecNoSpinUtils::invTOp( tmp2,
				       tmp1, 
				       u,
				       tsite,
				       P_mat,
				       P_mat_dag,
				       Q_mat_inv,
				       Q_mat_dag_inv,
				       invfact,
				       PLUS,
				       getTMax(),
				       schroedingerTP());
      
      T1  = spinReconstructDir3Minus(tmp2);
           
      tmp1=spinProjectDir3Minus(X);
      CentralTPrecNoSpinUtils::invTOp( tmp2,
				       tmp1, 
				       u,
				       tsite,
				       P_mat,
				       P_mat_dag,
				       Q_mat_inv,
				       Q_mat_dag_inv,
				       invfact,
				       MINUS,
				       getTMax(),
				       schroedingerTP());

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
    getFermBC().zero(ds_u);
    
  }


  //! Derivative of the derivCRight
  // 
  //    X^\dagger dC_R Y 
  //    = X^\dagger d [  P_{+} T^{-\dagger} ] Y 
  //    = - X^\dagger  P{+} T^{-\dagger} dT^\dagger T^{-\dagger} Y 
  //    = - X^\dagger  P{+} P_{+} T^{-\dagger } dT^\dagger T^{-\dagger} Y  since P_{+} P_{+} = P_{+}
  //    = - X^\dagger  T^{-\dagger} P_{+} dT^\dagger P_{+} T^{-\dagger} Y since Ps and Ts commute
  //    = - L^\dagger dT R
  // 
  //   with R = P_{+} T^{-\dagger} Y and L = P_{+} T^{-1} X
  //
  //   T^\dagger contains only U-dagger terms and dT^\dagger/dU  = 0
  //   T contains only U terms so  (dT^\dagger)^\dagger / dU R = dT/dU R = -shift(R, FORWARD, 1)
  //
  //   similarly to how we do dslash-es.
  //   hermitian conjugate terms which are automagically taken care of when we do the
  //   TAProj() later on  
  // 
  //   Actually this is just derivCLeft with PLUS --> MINUS in the Ifs
  //  and P_{-} <-> P_{+}
  // 
  void   
  UnprecSCprecTWilsonLinOp::derivCRight(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const
  {
    ds_u.resize(Nd);
    for(int mu=0; mu < 3; mu++) { 
      ds_u[mu] = zero;
    }
    
    if (isign == MINUS) { 
      
      LatticeHalfFermion tmp1, tmp2;
      T T1, T2; 
      
      tmp1=spinProjectDir3Plus(Y);
      CentralTPrecNoSpinUtils::invTOp( tmp2,
				       tmp1, 
				       u,
				       tsite,
				       P_mat,
				       P_mat_dag,
				       Q_mat_inv,
				       Q_mat_dag_inv,
				       invfact,
				       PLUS,
				       getTMax(),
				       schroedingerTP());
      
      T1  = spinReconstructDir3Plus(tmp2);
      
      tmp1=spinProjectDir3Plus(X);
      CentralTPrecNoSpinUtils::invTOp( tmp2,
				       tmp1, 
				       u,
				       tsite,
				       P_mat,
				       P_mat_dag,
				       Q_mat_inv,
				       Q_mat_dag_inv,
				       invfact,
				       MINUS,
				       getTMax(),
				       schroedingerTP());

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
    getFermBC().zero(ds_u);
  
  }

  //! derivSpace
  // 
  //  Easy peasy... call the force in the dslash.
  void   
  UnprecSCprecTWilsonLinOp::derivSpace(P& ds_u, const T& X, const T& Y, 
				       enum PlusMinus isign) const 
  {
    Real mhalf=Real(-0.5);
    ds_u.resize(Nd);
    
    Dw3D.deriv(ds_u, X, Y, isign);
    for(int mu = 0; mu < 3 ; mu++) {
      ds_u[mu] *= mhalf;
    }
    
    ds_u[3] = zero;
    getFermBC().zero(ds_u);
  }
  
  //! Apply the d/dt of the preconditioned linop
  //
  //  derivative of  C_L D_s C_R is evaluated by the chain rule:
  //
  //              = dC_L D_s C_R
  //              + C_L dD_s C_R
  //              + C_L D_s dC_R
  // 
  //   but we know that dC_R/ dU = 0 since C_R contains only U^\dagger terms
  //   so we can drop that term. Analogously in the daggered case we can drop
  //   the term containing C_L^\dagger...
  void   
  UnprecSCprecTWilsonLinOp::deriv(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const
  {
    T T1,T2,T3;
    
    P ds_tmp;
    ds_tmp.resize(Nd);
    
    
    if ( isign == PLUS ) { 
      
      // Derivative 
      cRightLinOp(T3, Y, PLUS);  // T3 = C_r Y
      spaceLinOp(T1, T3, PLUS);  // T1 = D_s C_r Y
      derivCLeft(ds_u, X, T1, PLUS);  // X^\dag dC_l T1

      // X^\dag dC_l D_s C_r Y
      
      cLeftLinOp(T2, X, MINUS);   // T2 = C_L^\dag X => T2^dag = X^\dag C_L
      derivSpace(ds_tmp, T2, T3, PLUS);  // X^\dag C_l dD_s C_r Y
      ds_u += ds_tmp;
    }
    else {
      
      
      // Derivative of the daggered
      cLeftLinOp(T3, Y, MINUS);   // T3 = C_l^\dag Y
      spaceLinOp(T1, T3, MINUS);   // T1 = D_s^\dag C_l^\dag Y
      derivCRight(ds_u, X, T1, MINUS); // X^\dag dC_r^\dag D_s^\dag C_l^\dag Y
      
      cRightLinOp(T2, X, PLUS);    // T2 = C_r X
      derivSpace(ds_tmp, T2, T3, MINUS); // X^\dag C_r^\dag dD_s^\dag C_l^\dag Y
      ds_u += ds_tmp;
      
    }
    getFermBC().zero(ds_u);
  }


  // Derivative of  trace log (T^\dagger T) -- Forward to central tprec utils
  void 
  UnprecSCprecTWilsonLinOp::derivLogDetTDagT(P& ds_u, 
					     enum PlusMinus isign) const
  {

    // Derivative of a Hermitian quantity so ignore isign?
    // Initial development -- set to 0
    ds_u.resize(Nd);
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu] = zero;
    }

    CentralTPrecNoSpinUtils::derivLogDet(ds_u, 
					 u,
					 Q_mat_inv,
					 tsite,
					 tDir(),
					 fact,
					 schroedingerTP());
    getFermBC().zero(ds_u);
		    
  }
  
} // End Namespace Chroma

#endif
#endif
#endif
