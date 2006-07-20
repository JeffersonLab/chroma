// $Id: mp_gaugeact.cc,v 1.2 2006-07-20 18:24:15 bjoo Exp $
/*! \file
 *  \brief Tree-level tadpole-improved Luscher-Weisz gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/mp_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace MPGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new MPGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
				MPGaugeActParams(xml, path));
    }

    const std::string name = "MORNINGSTAR_PEARDON_GAUGEACT";
    const bool registered = TheGaugeActFactory::Instance().registerObject(name, 
									  createGaugeAct);
  };


  MPGaugeActParams::MPGaugeActParams(XMLReader& xml_in, const std::string& path) {
    XMLReader paramtop(xml_in, path);

    try {
      read(paramtop, "./beta", beta);
      read(paramtop, "./u0s", u0s);
      read(paramtop, "./u0t", u0t);
      read(paramtop, "./omega", omega);

      //  Read optional anisoParam.
      if (paramtop.count("AnisoParam") != 0) {
	read(paramtop, "AnisoParam", aniso);
      }
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, const string& path, MPGaugeActParams& p) {
    MPGaugeActParams tmp(xml, path);
    p=tmp;
  }


  // Private initializer
  void
  MPGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
  {
    // Do the plaquette first. Spatial and temporal coeffs
    // anisotropy multiplied in in the terms constructor

    // Various tadpole things
    // spatial powers
    Real u_s_2 = param.u0s * param.u0s;
    Real u_s_4 = u_s_2 * u_s_2;
    Real u_s_6 = u_s_4 * u_s_2;
    Real u_s_8 = u_s_4 * u_s_4;

    // temporal powers
    Real u_t_2 = param.u0t * param.u0t;
    Real u_t_4 = u_t_2 * u_t_2;

    // Coefficients for the plaquette term (eq 4 in hep-lat/9911003)
    Real plaq_c_s = param.beta * Real(5) * ( Real(1) + param.omega ) / ( Real(3) * u_s_4 );
    Real plaq_c_t = param.beta * Real(4) / ( Real(3) * u_s_2 * u_t_2 );
    plaq = new PlaqGaugeAct(cgs, plaq_c_s, plaq_c_t, param.aniso);

    // Coefficients for the rectangle 
    Real rect_c_s = - param.beta / ( Real(12)*u_s_6 );
    
    // Loops that are short in the time direction
    Real rect_c_t_2 = - param.beta / ( Real(12)*u_s_4*u_t_2);

    // Loops that are long int the time direction ought to be ommitted
    bool no_temporal_2link = true;
    Real rect_c_t_1 = 0; // Specify a zero coefficient (skipped anyway)


    rect = new RectGaugeAct(cgs, rect_c_s, rect_c_t_1, rect_c_t_2, no_temporal_2link, param.aniso);
    

    // Coefficient of 2 plaquette spatial adjoint like thingie
    Real coeff_2plaq = Real(-5)*param.beta*param.omega/(Real(3)*u_s_8);
    two_plaq = new SpatialTwoPlaqGaugeAct(cgs, coeff_2plaq, param.aniso);

  } 

}

