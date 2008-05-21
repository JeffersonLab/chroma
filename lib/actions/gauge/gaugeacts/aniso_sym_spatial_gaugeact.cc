// $Id: aniso_sym_spatial_gaugeact.cc,v 3.3 2008-05-21 17:07:50 bjoo Exp $
/*! \file
 *  \brief Anisotropic gaugeact useful for spectrum from hep-lat/9911003
 *
 *  Tree-level LW with tapole improvement, missing 1x2 in time, also including
 *  2-plaq term. Taken from Morningstar-Peardon, hep-lat/9911003
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/aniso_sym_spatial_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "actions/gauge/gaugeacts/aniso_sym_shared_functions.h"

#include <cstdio>
using namespace std;

namespace Chroma
{
 
  namespace AnisoSymSpatialGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new AnisoSymSpatialGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
				       AnisoSymSpatialGaugeActParams(xml, path));
    }

    const std::string name = "ANISO_SYM_SPATIAL_GAUGEACT";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeActFactory::Instance().registerObject(name, createGaugeAct);
	registered = true;
      }
      return success;
    }

  }

  Double AnisoSymSpatialGaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
  {
    START_CODE();

    LatticeReal lgimp=zero;
    const multi1d<LatticeColorMatrix>& u_bc = state->getLinks();

    
    for(int mu = 0; mu < Nd; mu++) { 
      for(int nu = 0 ; nu < Nd; nu++) { 
	if ( ( mu != nu ) && (mu != param.aniso.t_dir) 
	     && ( nu != param.aniso.t_dir ) ) {
      
	  AnisoSym::S_part(mu, 
			   nu, 
			   param.aniso.t_dir,
			   plaq_c_s, 
			   rect_c_s, 
			   true,
			   lgimp,
			   u_bc);

	}
      }
    }
    
    // If user gave a zero point energy, subtract it off
    if( param.use_subtraction  ) { 
      LatticeReal ff;
      ff = param.sub_zero;
      lgimp -= ff;
    }

    // Sum action and normalize out.
    Double ret_val = sum(lgimp);
    ret_val *= -Double(1)/Double(Nc);

    END_CODE();

    return ret_val;
  }


  //! Compute dS/dU
  void AnisoSymSpatialGaugeAct::deriv(multi1d<LatticeColorMatrix>& result,
	     const Handle< GaugeState<P,Q> >& state) const 
  {
    START_CODE();

    const multi1d<LatticeColorMatrix>& u_bc = state->getLinks();
    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    result.resize(Nd);
    for(int mu =0; mu < Nd; mu++) { 
      result[mu] = zero;
      ds_tmp[mu] = zero;
    }
    



    for(int mu = 0; mu < Nd; mu++) { 
      for(int nu = 0 ; nu < Nd; nu++) { 

	// mu and nu both have to be spatial and not equal to each other
	// It is OK to accumulate into ds_tmp
	if( ( mu != nu ) && (mu != param.aniso.t_dir) 
	    && (nu != param.aniso.t_dir ) ) {


	  AnisoSym::deriv_part(mu, 
			       nu, 
			       param.aniso.t_dir,
			       plaq_c_s, 
			       rect_c_s, 
			       true,
			       ds_tmp,
			       u_bc);


	  
	}
      }
    }

    // Close loops
    for(int mu=0 ; mu < Nd; mu++) { 
      result[mu] = u_bc[mu]*ds_tmp[mu];
    }

    // Apply boundaries
    getGaugeBC().zero(result);
    END_CODE();

  }

  // Private initializer
  void
  AnisoSymSpatialGaugeAct::init(void)
  {
    START_CODE();

    // Do the plaquette first. Spatial and temporal coeffs
    // anisotropy multiplied in in the terms constructor

    // Various tadpole things
    // spatial powers
    Real u_s_2 = param.u_s * param.u_s;
    Real u_s_4 = u_s_2 * u_s_2;
    Real u_s_6 = u_s_4 * u_s_2;

    // Compute coefficient of spatial plaquettes and rectangles
    plaq_c_s = param.beta * Real(5)/( Real(3)* u_s_4 );
    rect_c_s =  - param.beta / ( Real(12)*u_s_6 );

    // Now take care of anisotropy
    if ( param.aniso.anisoP == true ) { 
      plaq_c_s /= param.aniso.xi_0;
      rect_c_s /= param.aniso.xi_0;
    }
    
    END_CODE();
  } 

}

