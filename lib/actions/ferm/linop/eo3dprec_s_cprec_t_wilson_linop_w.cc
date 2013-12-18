// $Id: eo3dprec_s_cprec_t_wilson_linop_w.cc,v 1.4 2008-05-07 01:12:18 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "chromabase.h"
#include "actions/ferm/linop/eo3dprec_s_cprec_t_wilson_linop_w.h"
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
  void EO3DPrecSCprecTWilsonLinOp::create(Handle< FermState<T,P,Q> > fs_,
					const Real& Mass_,
					const AnisoParam_t& anisoParam_)
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    // Check we are in 4D
    if ( Nd != 4 ) { 
      QDPIO::cout << "This class (EO3DPrecSCprecTWilsonLinOp) only works in 4D" << endl;
      QDP_abort(1);
    }
    
    // Check Aniso Direction -- has to be 3.
    if ( anisoParam_.t_dir != 3 ) { 
      QDPIO::cout << "This class (EO3DPrecSCprecTWilsonLinOp) is hardwired for t_dir=3"<< endl;
      QDP_abort(1);
    }

    // Check that the time direction is local:
    const multi1d<int>& s_size =  QDP::Layout::subgridLattSize();  // Local Lattice
    const multi1d<int>& t_size =  QDP::Layout::lattSize(); // Total Latt Size
    if( t_size[3] != s_size[3] ) { 
      QDPIO::cout << "This class (EO3DPrecSCprecTWilsonLinOp) needs time to be local" << endl;
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
	CMat tmp = one.elem() + P_mat(cb3, site,0);
	CentralTPrecNoSpinUtils::invert3by3( Q_mat_inv(cb3, site), tmp );
	
	// Compute Q_dag = 1 + W^\dag P^\dag, W = [ 0, ..., 1 ]  = > Q_dag = P^\dag [Nt-1]
	// Similar to eq 43
	// NB: This is not necessarily SU(3) now, so we can't just take the dagger to get the inverse
	tmp = one.elem() + P_mat_dag(cb3, site,Nt-1);
	CentralTPrecNoSpinUtils::invert3by3( Q_mat_dag_inv(cb3, site), tmp);
      }
    }


    Dw3D.create(fs_, anisoParam_);

    END_CODE();
#endif
  }

  //! Apply (C_L)^{-1}
  void EO3DPrecSCprecTWilsonLinOp::invCLeftLinOp(T& chi, 
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
    getFermBC().modifyF(chi, rb3[cb3d]);
  } 

  //! Apply (C_R)^{-1}
  void EO3DPrecSCprecTWilsonLinOp::invCRightLinOp(T& chi, 
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
    getFermBC().modifyF(chi, rb3[cb3d]);
  } 


  //! Apply C_L
  void EO3DPrecSCprecTWilsonLinOp::cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign,
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
				     isign,
				     getTMax());


    chi[rb3[cb3d]] += spinReconstructDir3Minus(tmp_T);
    chi[rb3[cb3d]] *= Real(0.5);
    getFermBC().modifyF(chi, rb3[cb3d]);
  }

  //! Apply C_R
  void EO3DPrecSCprecTWilsonLinOp::cRightLinOp(T& chi, const T& psi, enum PlusMinus isign,
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
				     other_sign,
				     getTMax());


    chi[rb3[cb3d]] += spinReconstructDir3Plus(tmp_T);
    chi[rb3[cb3d]] *= Real(0.5);
    getFermBC().modifyF(chi, rb3[cb3d]);
  }

} // End Namespace Chroma

#endif
#endif
#endif
