// $Id: ilu2prec_s_cprec_t_clover_linop_w.cc,v 3.1 2008-10-08 19:40:17 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Clover linear operator
 */


#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "chromabase.h"
#include "actions/ferm/linop/ilu2prec_s_cprec_t_clover_linop_w.h"
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
  void ILU2PrecSCprecTCloverLinOp::create(Handle< FermState<T,P,Q> > fs_,
					const CloverFermActParams& param_)
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    // Copy Params
    param = param_;

    // Check we are in 4D
    if ( Nd != 4 ) { 
      QDPIO::cout << "This class (ILU2PrecSCprecTCloverLinOp) only works in 4D" << endl;
      QDP_abort(1);
    }
    
    // Check Aniso Direction -- has to be 3.
    if ( param.anisoParam.anisoP== true && param.anisoParam.t_dir != 3 ) { 
      QDPIO::cout << "This class (ILU2PrecSCprecTCloverLinOp) is hardwired for t_dir=3"<< endl;
      QDP_abort(1);
    }

    // Check that the time direction is local:
    const multi1d<int>& s_size =  QDP::Layout::subgridLattSize();  // Local Lattice
    const multi1d<int>& t_size =  QDP::Layout::lattSize(); // Total Latt Size
    if( t_size[3] != s_size[3] ) { 
      QDPIO::cout << "This class (ILU2PrecSCprecTCloverLinOp) needs time to be local" << endl;
      QDP_abort(1);
    }

    // Store gauge state etc
    fs = fs_;
    
    // Work out aniso factors
    Real ff = where(param.anisoParam.anisoP, param.anisoParam.nu / param.anisoParam.xi_0, Real(1));
    fact = 1 + (Nd-1)*ff + param.Mass;
    QDPIO::cout << "Factor=" << fact << "  log10(Factor)="<< log10(fact) << endl;

    invfact = Real(1)/fact;
    u = fs->getLinks();

    // Incorporate into links
    if (param.anisoParam.anisoP)
    {
      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != param.anisoParam.t_dir)
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
    QDPIO::cout << "Got here " << endl << flush;
    if( !schrTP ) { 
      // If we are using the max_norm tric. Compute the t_needed

      if( param.max_norm_usedP) {
	Real tmp1=param.max_norm/sqrt(Real(3));
	Real t = -( log(tmp1)/log(fact));
	QDPIO::cout << "Cutoff t is " << t<< endl;
	t_max=(int)toDouble(ceil(t));
	QDPIO::cout << "T_max="<<t_max << endl;
      }
      else {
	t_max=Nt;
	QDPIO::cout << "Not using cutoff trick. Setting T_max="<<t_max<<endl;
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
	    if( t > Nt-1-t_max) {
	      P_mat(cb3, site, t) =  u[t_index].elem( tsite(cb3, site,t)  )  * P_mat(cb3, site,t+1);
	      P_mat(cb3, site, t) *= invfact.elem();
	    }
	    else { 
	      zero_rep(P_mat(cb3, site, t));
	    }
	  }
	  // Compute the dagger. Opposite order (forward sub) similara to eq 38-42
	   // in notes
	  P_mat_dag(cb3, site, 0) = adj( u[t_index].elem( tsite(cb3, site, Nt-1) ) );
	  P_mat_dag(cb3, site, 0) *= minvfact.elem();
	 
	  for(int t=1; t < Nt; t++) { 
	    if( t < t_max ) { 
	      P_mat_dag(cb3, site,t) = adj( u[t_index].elem(  tsite(cb3, site, t-1) ) )*P_mat_dag(cb3,site, t-1) ;
	      P_mat_dag(cb3, site,t) *= invfact.elem();
	    }
	    else {
	      zero_rep(P_mat_dag(cb3, site,t));
	    }
	    
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
    }
    else {
      t_max = Nt;
    }
      
    // Create Clover term = A + factor (we don't want the diag mass bit)
    APlusFact.create(fs, param);
    Dw3D.create(fs_, param_.anisoParam);

#if 1
    QDPIO::cout << "schrTP is " << schrTP << endl;
    if(schrTP==false) {
     // Numerical experiments
      for(int i=0; i < Nt; i++) { 
      Double fnorm=Double(0);
      Double fnorm2=Double(0);
      
      OScalar<CMat> f; f.elem() = P_mat(0,0,i);
      fnorm=sqrt(norm2(f));
      OScalar<CMat> f2; f2.elem() = P_mat_dag(0,0,i);
      fnorm2=sqrt(norm2(f2));
          
      QDPIO::cout << "cb=0 x=0 t="<<i<<" normP="<<fnorm<<" normPdag="<<fnorm2 << endl;
    }
    }
#endif		 

#endif
    END_CODE();
  }


} // End Namespace Chroma

#endif
#endif
#endif
