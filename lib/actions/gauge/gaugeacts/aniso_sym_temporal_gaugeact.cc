// $Id: aniso_sym_temporal_gaugeact.cc,v 3.1 2007-05-22 14:19:42 bjoo Exp $
/*! \file
 *  \brief  Temporal part of Anisotropic Tree leve LW gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/aniso_sym_temporal_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace AnisoSymTemporalGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new AnisoSymTemporalGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
				       AnisoSymGaugeActParams(xml, path));
    }

    const std::string name = "ANISO_SYM_TEMPORAL_GAUGEACT";

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

  // Private initializer
  void
  AnisoSymTemporalGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
  {
    START_CODE();

    // Do the plaquette first. Spatial and temporal coeffs
    // anisotropy multiplied in in the terms constructor

    // Various tadpole things
    // spatial powers
    Real u_s_2 = param.u_s * param.u_s;
    Real u_s_4 = u_s_2 * u_s_2;
    Real u_s_6 = u_s_4 * u_s_2;
    Real u_s_8 = u_s_4 * u_s_4;

    // temporal powers
    Real u_t_2 = param.u_t * param.u_t;
    Real u_t_4 = u_t_2 * u_t_2;

    // Coefficients for the plaquette term (eq 4 in hep-lat/9911003)
    // Space Space terms do not contribute.
    // Real plaq_c_s = param.beta * Real(5) / ( Real(3) * u_s_4 );
    Real plaq_c_s = Real(0);

    // Parameter needed for gauge act, but should never be used
    Real plaq_c_t = param.beta * Real(4) / ( Real(3) * u_s_2 * u_t_2 );
    
    plaq = new PlaqGaugeAct(cgs, plaq_c_s, plaq_c_t, param.aniso);

    
    // Coefficients for the rectangle 
    // Space-Space terms do not contribute
    //Real rect_c_s = - param.beta / ( Real(12)*u_s_6 );
    Real rect_c_s = Real(0);

    // Loops that are short in the time direction
    // Param needed for rect_gaugeact driver, but never used
    Real rect_c_t_2 = - param.beta / ( Real(12)*u_s_4*u_t_2);
    

    // Loops that are long int the time direction ought to be ommitted
    bool no_temporal_2link = true;
    Real rect_c_t_1 = 0; // Specify a zero coefficient (skipped anyway)

    rect = new RectGaugeAct(cgs, rect_c_s, rect_c_t_1, rect_c_t_2, no_temporal_2link, param.aniso);
    
    END_CODE();
  } 

}

