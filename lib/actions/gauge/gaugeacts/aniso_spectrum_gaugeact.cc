// $Id: aniso_spectrum_gaugeact.cc,v 1.7 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Anisotropic gaugeact useful for spectrum from hep-lat/9911003
 *
 *  Tree-level LW with tapole improvement, missing 1x2 in time, also including
 *  2-plaq term. Taken from Morningstar-Peardon, hep-lat/9911003
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/aniso_spectrum_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace AnisoSpectrumGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, 
		 multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
							       const std::string& path) 
    {
      return new AnisoSpectrumGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
				       AnisoSpectrumGaugeActParams(xml, path));
    }

    const std::string name = "ANISO_SPECTRUM_GAUGEACT";

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


  AnisoSpectrumGaugeActParams::AnisoSpectrumGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try 
    {
      read(paramtop, "beta", beta);
      read(paramtop, "u_s", u_s);
      read(paramtop, "u_t", u_t);
      read(paramtop, "omega", omega);
      read(paramtop, "AnisoParam", aniso);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }


  void read(XMLReader& xml, const string& path, AnisoSpectrumGaugeActParams& p) 
  {
    AnisoSpectrumGaugeActParams tmp(xml, path);
    p=tmp;
  }

  void write(XMLWriter& xml, const string& path, const AnisoSpectrumGaugeActParams& param) 
  {
    push(xml, path);

    write(xml, "beta", param.beta);
    write(xml, "u_s", param.u_s);
    write(xml, "u_t", param.u_t);
    write(xml, "omega", param.omega);
    write(xml, "AnisoParam", param.aniso);

    pop(xml);
  }


  // Private initializer
  void
  AnisoSpectrumGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
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
    Real plaq_c_s = param.beta * Real(5) * ( Real(1) + param.omega ) / ( Real(3) * u_s_4 );
    Real plaq_c_t = param.beta * Real(4) / ( Real(3) * u_s_2 * u_t_2 );

    //    plaq = new PlaqGaugeAct(cgs, plaq_c_s, plaq_c_t, param.aniso);

    // Coefficient of 2 plaquette spatial adjoint like thingie
    Real coeff_2plaq = Real(-5)*param.beta*param.omega/(Real(3)*u_s_8);
    plaq_plus_two_plaq = new PlaqPlusSpatialTwoPlaqGaugeAct(cgs, 
							    plaq_c_s,
							    plaq_c_t,
							    coeff_2plaq, 
							    param.aniso);
    // Coefficients for the rectangle 
    Real rect_c_s = - param.beta / ( Real(12)*u_s_6 );
    
    // Loops that are short in the time direction
    Real rect_c_t_2 = - param.beta / ( Real(12)*u_s_4*u_t_2);

    // Loops that are long int the time direction ought to be ommitted
    bool no_temporal_2link = true;
    Real rect_c_t_1 = 0; // Specify a zero coefficient (skipped anyway)

    rect = new RectGaugeAct(cgs, rect_c_s, rect_c_t_1, rect_c_t_2, no_temporal_2link, param.aniso);
    
    END_CODE();
  } 

}

