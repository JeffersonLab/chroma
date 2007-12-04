// $Id: unprec_s_cprec_t_wilson_linop_w.cc,v 1.3 2007-12-04 16:04:42 bjoo Exp $
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
      CMat tmp = one.elem() + P_mat(site,0);
      CentralTPrecNoSpinUtils::invert3by3( Q_mat_inv(site), tmp );
      
      // Compute Q_dag = 1 + W^\dag P^\dag, W = [ 0, ..., 1 ]  = > Q_dag = P^\dag [Nt-1]
      // Similar to eq 43
      // NB: This is not necessarily SU(3) now, so we can't just take the dagger to get the inverse
      tmp = one.elem() + P_mat_dag(site,Nt-1);
      CentralTPrecNoSpinUtils::invert3by3( Q_mat_dag_inv(site), tmp);
    }

    Dw3D.create( fs_, anisoParam_);

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

    // Here I project - apply TOp to only the halfector - then reco nstruct
    tmp_minus = spinProjectDir3Minus(psi);


    // Use shared routine to apply T or T^\dagger.
    // Pass in u, and tsite for the subset
    CentralTPrecNoSpinUtils::TOp(tmp_T, 
				 tmp_minus, 
				 u,
				 tsite,
				 fact,
				 isign);

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


    // Use shared routine to apply T or T^\dagger.
    // Pass in u, and tsite for the subset
    CentralTPrecNoSpinUtils::TOp(tmp_T, 
				 tmp_plus, 
				 u,
				 tsite,
				 fact,
				 other_sign);

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
				     isign);


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
				     other_sign);


    chi += spinReconstructDir3Plus(tmp_T);
    chi *= Real(0.5);
  }


} // End Namespace Chroma

#endif
#endif
#endif
