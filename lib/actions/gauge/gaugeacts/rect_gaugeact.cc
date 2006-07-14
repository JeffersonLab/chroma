// $Id: rect_gaugeact.cc,v 3.1 2006-07-14 20:33:17 bjoo Exp $
/*! \file
 *  \brief Rectangle gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace RectGaugeActEnv 
  {
    //! Callback
    GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
											    const std::string& path) 
    {
      return new RectGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			      RectGaugeActParams(xml, path));
    }

    const std::string name = "RECT_GAUGEACT";
    const bool registered = TheGaugeActFactory::Instance().registerObject(name, 
									  createGaugeAct);
  };


  // Param constructor/reader
  RectGaugeActParams::RectGaugeActParams(XMLReader& xml_in, const std::string& path) {
    XMLReader paramtop(xml_in, path);

    try {
      
      // Read AnisoParam if it exists 
      if( paramtop.count("AnisoParam") != 0 ) {
	read(paramtop, "AnisoParam", aniso);
      }
      else { 
	// Set aniso to not be anisotropic if it doesn't
	aniso.anisoP = false;
	aniso.t_dir = Nd -1;
	aniso.xi_0 = 1;
	aniso.nu = 1;
      }


      if ( aniso.anisoP == false ) { 
	
	// If we are not anisotropic look for coeff 
	// and set both coeff_s = coeff_t = coeff, do temporal = true

	read(paramtop, "coeff", coeff_s);
	coeff_t1 = coeff_s;
	coeff_t2 = coeff_s;
	
      }
      else { 

	// If we are ansotropic, look for coeff_s and coeff_t
	read(paramtop, "coeff_s", coeff_s);
	read(paramtop, "coeff_t1", coeff_t1);
	read(paramtop, "coeff_t2", coeff_t2);

      }

      // In either anisotropic or isotropic case we can optionally
      // specify whether we want to omit temporal rectangles
      if( paramtop.count("no_temporal_2link") == 0 ) {
	// Default: always do temporal
	no_temporal_2link = false;
      }
      else { 
	// If it is specified in the file, then read it.
	read(paramtop, "no_temporal_2link", no_temporal_2link);
      }

    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }

  // Read params
  void read(XMLReader& xml, const string& path, RectGaugeActParams& p) 
  {
    RectGaugeActParams tmp(xml, path);
    p=tmp;
  }


  //! Compute staple
  /*!
   * \param u_staple   result      ( Write )
   * \param state      gauge field ( Read )
   * \param mu         direction for staple ( Read )
   * \param cb         subset on which to compute ( Read )
   */
  void
  RectGaugeAct::staple(LatticeColorMatrix& u_staple,
		       const Handle< GaugeState<P,Q> >& state,
		       int mu, int cb) const
  {
    QDPIO::cerr << "RectGaugeAct::staple() - not converted from szin" << endl;
    QDP_abort(1);

    END_CODE();
  }



  //! Computes the derivative of the fermionic action respect to the link field
  /*!
   *         |  dS
   * ds_u -- | ----   ( Write )
   *         |  dU  
   *
   * \param ds_u       result      ( Write )
   * \param state      gauge field ( Read )
   */
  void
  RectGaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
		      const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    ds_u.resize(Nd);

    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_3;
    LatticeColorMatrix tmp_4;
    LatticeColorMatrix tmp_5;
    LatticeColorMatrix tmp_6;
    LatticeColorMatrix tmp_tot;
    LatticeColorMatrix tmp_munu;

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    ds_u = zero;

    // It is 1/(4Nc) to account for normalisation relevant to fermions
    // in the taproj, which is a factor of 2 different from the 
    // one used here.

    Real coeff_tmp = Real(-1)/Real(2*Nc);

    for(int mu = 0; mu < Nd; ++mu) {

      tmp_tot = zero;	

      Real coeff_mu;	

      // The way this routine works:
      // 
      // Compute the contribs of the 2x1 rectangles to ds_u[mu]
      // Compute the contribs of the 1x2 rectangles to ds_u[nu]
      
      // So if mu == t_dir ds_u[mu] only gets contribs from loops 
      // that have an extent 2 in the t_dir ==> coeff_mu = coeff_t1
      if ( mu == params.aniso.t_dir ) {
	coeff_mu = params.coeff_t1;
      }
      else {
	coeff_mu = params.coeff_s;
      }
      
      
      bool skip;

      // This is a flag to guard us from doing extra work. 
      // When mu == t_dir and no_temporal_2link == true
      // this should be set to false

      if( mu == params.aniso.t_dir && params.no_temporal_2link == true ) {
	skip = true;
      }
      else {
	skip = false;
      }
      
      
      /*  2*1 rectangles */
      /* calculate double_mu links */
      tmp_3 = u[mu] * shift(u[mu], FORWARD, mu);
      
      for(int cb=0; cb < 2; cb++) {
	
	
	tmp_6[rb[1-cb]] = tmp_3 * coeff_tmp;

	for(int j=0, nu=0; nu < Nd; nu++) {
	  // The way this routine works:
	  // 
	  // Compute the contribs of the 2x1 rectangles to ds_u[mu]
	  // Compute the contribs of the 1x2 rectangles to ds_u[nu]
	  
	  // So if mu == t_dir ds_u[mu] only gets contribs from loops 
	  // that have an extent 2 in the t_dir ==> coeff_mu = coeff_t1
	  
	  // Likewise if nu == t_dir d_u[nu] only gets contribs from 
	  // loops that have an extent 1 in the t_dir .==> coeff_nu = coeff_t2
	  // In all other cases, mu and nu are spatial so we use the 
	  // coeff_s;
	  Real coeff_nu;
	  
	  if( nu == params.aniso.t_dir ) { 
	    coeff_nu = params.coeff_t2;
	  }
	  else { 
	    coeff_nu = params.coeff_s;
	  }
	  
	  
	  
	  if(nu == mu)
	    continue;
	  
	  /* forward plaquette */
	  // Doing this for ds_u[nu]
	  //
	  //     <--- <-----
	  //     |           ^
	  //     |           . U_{nu} (x)  (coeff_nu) 
	  //     V           . 
	  //     ---> -----> . 
	  tmp_1[rb[cb]] = u[nu] * shift(tmp_6, FORWARD, nu);
	  tmp_5[rb[1-cb]] = shift(adj(tmp_1)*tmp_3, BACKWARD, mu);
	  
	  ds_u[nu][rb[cb]] += coeff_nu*u[nu]*shift(tmp_5, BACKWARD, mu);
	  
	  
	  
	  tmp_5[rb[1-cb]] = shift(u[nu], FORWARD, mu);
	  
	  /* at this point we add the nu contribution directly to ds_u */
	  /* we could make tmp_tot carry a direction index in order to avoid that */
	  if(j++ == 0) {
	    
	    
	    //        ^---> --->
	    //        |        |
	    // tmp2 = |        |
	    //                 V 	  
	    tmp_2[rb[cb]] = tmp_1 * adj(shift(tmp_5, FORWARD, mu));
	    
	    if( !skip ) {
	      // Start off the sum of 2x1 rectangle staples if 
	      // we are doing those
	      tmp_4[rb[cb]] = tmp_2;
	    }
	    //          -----> --->
	    //          ^         |
	    // U(x,nu)  .         |     (coeff[nu][mu]
	    //          .         |
	    //          <--- <--- V
	    ds_u[nu][rb[cb]] += coeff_nu*tmp_2*adj(tmp_3);
	  }
	  else {
	    
	    //         ^---> --->
	    //         |        |
	    // tmp4 += |        |
	    //                 V
	    
	    tmp_2[rb[cb]] = tmp_1 * adj(shift(tmp_5, FORWARD, mu));
	    
	    if( !skip ) { 
	      // Add to the sum of 2x1 rectangle staples if we are doing
	      // those
	      tmp_4[rb[cb]] += tmp_2; /* sum to the staple */
	    }
	    
	    ds_u[nu][rb[cb]] += coeff_nu*tmp_2 * adj(tmp_3);
	  }
	  
	  /* backward plaquette */
	  tmp_1[rb[1-cb]] = adj(u[nu]) * tmp_6;
	  tmp_2[rb[cb]] = shift(u[nu], FORWARD, mu);
	  tmp_5[rb[1-cb]] = shift(tmp_2, FORWARD, mu);
	  
	  /* add   holds: |        ^ to tmp 4
	     |        |  
	     V --> --->          */
	  
	  if( !skip ) {
	    // Add downwards staple to the sum of 2x1 staples if 
	    // we are doing those
	    tmp_4[rb[cb]] += shift(tmp_1*tmp_5, BACKWARD, nu);
	  }
	  
	  
	} // End of nu loop
	
	// At this point tmp 4 holds sum_{nu} 2x1 flat staples
	// and we need to partially 'close them off' with U{mu} if we 
	// are doing that direction
	if( !skip ) {
	  tmp_tot[rb[cb]] += shift(u[mu], FORWARD, mu) * adj(tmp_4);
	  tmp_tot[rb[1-cb]] += shift(adj(tmp_4)*u[mu], BACKWARD, mu);
	}
      } // end of cb loop
      
	/* ds_u =  u(x,mu) * tmp_tot */
      if( !skip ) { 
	// Here we close off the things in tmp_tot_completely.
	ds_u[mu] += coeff_mu*u[mu] * tmp_tot;
      }
      
    }
  
  
    // Zero the force on any fixed boundaries
    getGaugeBC().zero(ds_u);
    
    END_CODE();
  }


  // Get the gauge action
  //
  // S = -(coeff/(Nc) Sum Re Tr Rect
  //
  // w_rect is defined in MesPlq as
  //
  // w_rect =( 2/(V*Nd*(Nd-1)*Nc)) * Sum Re Tr Rect
  //
  // so 
  // S = -coeff * (V*Nd*(Nd-1)/2) w_rect 
  //   = -coeff * (V*Nd*(Nd-1)/2)*(2/(V*Nd*(Nd-1)*Nc))* Sum Re Tr Rect
  //   = -coeff * (1/(Nc)) * Sum Re Tr Rect

Double
RectGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
{
  START_CODE();
  
  const multi1d<LatticeColorMatrix>& u = state->getLinks();
  
  LatticeColorMatrix tmp_0;
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;
  LatticeColorMatrix tmp_tot;
  LatticeReal lgimp = zero;
  
  // 2x1 rectangle piece
  
  
  for(int mu = 1; mu < Nd; ++mu)  {

    for(int nu = 0; nu < mu; ++nu) { 
	
      Real coeff = params.coeff_s; // The weight of the loop in the action 
      // Initialize to spatial
      
      for(int cb = 0; cb < 2; ++cb)  {
	
	// Going t do 1x2 loops (mu, nu, nu)
	// 
	//      <---^
	//      |   |    nu
	//      |   |     ^
	//      V   ^     |
	//      |   |     |
	//      |   |     |----> mu
	//      V--->
	
	// Mu is the short direction. If it is temporal
	// we use the coefficient for temporal 1x2 loops (c2t)
	//
	// if on the other hand nu (the long direction) is temporal
	// we use the coefficient for temporal 2x1 loops (c1t)
	// This overrides the setting above. Luckily 
	// the way the sum is done, we can never have mu == nu
	// so the cases don't conflict
	//
	// Additionally if nu is temporal, and no_temporal2_link
	// is turned on we skip this loop entirely
	bool skip = ( nu == params.aniso.t_dir && params.no_temporal_2link == true );
	if ( !skip ) {
	  
	  if( mu == params.aniso.t_dir ) {
	    coeff = params.coeff_t2;    // Mu is length 1 in time
	  }
	  else if( nu == params.aniso.t_dir ) {
	    coeff = params.coeff_t1;   
	  }
	  else {
	    coeff = params.coeff_s;
	  }
	  
	  /* tmp_1(x) = u(x-nu,nu) */
	  tmp_1[rb[cb]] = shift(u[nu], BACKWARD, nu);
	  
	  /* tmp_0 = tmp_1 * u(x,nu) = u(x,nu)*u(x+nu,nu) */
	  tmp_0[rb[cb]] = tmp_1 * u[nu];
	  
	  /* tmp_2(x) = tmp_0(x+mu) */
	  tmp_2[rb[1-cb]] = shift(tmp_0, FORWARD, mu);
	  
	  /* tmp_1(x) = u(x+nu,mu) */
	  tmp_1[rb[1-cb]] = shift(u[mu], FORWARD, nu);
	  
	  /* tmp_0 = tmp_2 * tmp_1_dag = u(x+mu-nu,nu)*u(x+mu,nu)*u_dag(x+nu,mu) */
	  tmp_0[rb[1-cb]] = tmp_2 * adj(tmp_1);
	  
	  /* tmp_1 = tmp_0 * u_dag(x,nu) */
	  /*       = u(x+mu-nu,nu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	  tmp_1[rb[1-cb]] = tmp_0 * adj(u[nu]);
	  
	  /* tmp_2(x) = tmp_1(x+nu) */
	  tmp_2[rb[cb]] = shift(tmp_1, FORWARD, nu);
	  
	  /* tmp_tot = tmp_2 * u_dag(x,nu) */
	  /*         = u(x+mu,nu)*u(x+mu+nu,nu)*u_dag(x+2*nu,mu) */
	  /*          *u_dag(x+nu,nu)*u_dag(x,nu) */
	  
	  
	  tmp_tot[rb[cb]] = coeff * tmp_2 * adj(u[nu]);
	}
	else { 
	  // If the case above is skipped we should zero tmp_tot[rb[cb]]
	  tmp_tot[rb[cb]] = zero; 
	}
	
	
	
	// Going t do 2x1 loops (mu, mu, nu)
	// 
	//                   ^ nu
	//  <---- <----^     |
	//  |          |     |
	//  |          |     ---> mu
	//  V ---> ---->
	
	
	// Now Mu is the long direction. If it is temporal
	// we use the coefficient for temporal 2x1 loops (c1t)
	// Nu is the short direction. If it is temporal
	// we use the coefficient for temporal 1x2 loops (c2t)
	// This overrides the setting above. Luckily 
	// the way the sum is done, we can never have mu == nu
	// so the cases don't conflict
	// Additionally if mu is temporal, and no_temporal2_link
	// is turned on we skip this loop entirely
	
	skip = ( mu == params.aniso.t_dir && params.no_temporal_2link == true);
	if (! skip ) {
	  
	  if( mu == params.aniso.t_dir ) {
	    coeff = params.coeff_t1;   
	  }
	  else if( nu == params.aniso.t_dir ) {
	    coeff = params.coeff_t2;   
	  }
	  else {
	    coeff = params.coeff_s;
	  }
	  
	  /* tmp_1(x) = u(x+mu,nu) */
	  tmp_1[rb[1-cb]] = shift(u[nu], FORWARD, mu);
	  
	  /* tmp_0 = u(x,nu) * tmp_1 = u(x,mu)*u(x+mu,nu) */
	  tmp_0[rb[1-cb]] = u[mu] * tmp_1;
	  
	  /* tmp_2(x) = tmp_0(x-nu) */
	  tmp_2[rb[cb]] = shift(tmp_0, BACKWARD, nu);
	  
	  /* tmp_0 = tmp_2 * u_dag(x,mu) = u(x-nu,mu)*u(x+mu-nu,nu)*u_dag(x,mu) */
	  tmp_0[rb[cb]] = tmp_2 * adj(u[mu]);
	  
	  /* tmp_1(x) = tmp_0(x+mu) */
	  tmp_1[rb[1-cb]] = shift(tmp_0, FORWARD, mu);
	  
	  /* tmp_0 = tmp_1 * u_dag(x,mu) */
	  /*       = u(x+mu-nu,mu)*u(x+2*mu-nu,nu)*u_dag(x+mu,mu)*u_dag(x,mu) */
	  tmp_0[rb[1-cb]] = tmp_1 * adj(u[mu]);
	  
	  /* tmp_2(x) = tmp_0(x+nu) */
	  tmp_2[rb[cb]] = shift(tmp_0, FORWARD, nu);
	  
	  /* tmp_tot += tmp_2 * u_dag(x,nu) */
	  /*         += u(x+mu,mu)*u(x+2*mu,nu)*u_dag(x+mu+nu,mu)*u_dag(x+nu,mu) */
	  /*           *u_dag(x,nu) */
	  tmp_tot[rb[cb]] += coeff * tmp_2 * adj(u[nu]);
	} // End of skip
	    
	/* lgimp += c_1(mu,nu)*trace(u(x,mu) * tmp_tot) */
	lgimp[rb[cb]] += real(trace(u[mu] * tmp_tot));
      }
    } // end nu loop
  } // end mu loop
  
  Double S_rect = sum(lgimp);
  S_rect *= -Real(1) / Real(Nc);   // note sign here
  
  END_CODE();
  
  return S_rect;
} 


// Isotropic case
RectGaugeAct::RectGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_,
			     const Real& coeff_s_) : cgs(cgs_)
{
  // Setup the params structure
  params.coeff_s = params.coeff_t1 = params.coeff_t2 = coeff_s_;
  params.no_temporal_2link = false;
  params.aniso.anisoP = false;
  params.aniso.t_dir = Nd -1;
  params.aniso.xi_0 = 1;
  params.aniso.nu = 1;    
}

// Anisotropic case
RectGaugeAct::RectGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_,
			   const Real& coeff_s_,
			   const Real& coeff_t1_,
			   const Real& coeff_t2_,
			   const bool no_temporal_2link_,
			   const AnisoParam_t& aniso_) : cgs(cgs_) 
{
  // Setup the params structure
  params.coeff_s = coeff_s_ / aniso_.xi_0 ;
  params.coeff_t1 = coeff_t1_ * aniso_.xi_0 ;
  params.coeff_t2 = coeff_t2_ * aniso_.xi_0;
  params.no_temporal_2link = no_temporal_2link_ ;
  params.aniso = aniso_ ;          // Rely on struct copy constructor here

}

//! Read rectangle coefficient from a param struct
RectGaugeAct::RectGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
			   const RectGaugeActParams& p) : cgs(cgs_), params(p) {}

}

