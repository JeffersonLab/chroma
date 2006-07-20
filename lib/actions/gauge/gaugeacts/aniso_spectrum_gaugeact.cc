// $Id: aniso_spectrum_gaugeact.cc,v 1.1 2006-07-20 18:40:35 edwards Exp $
/*! \file
 *  \brief Anisotropic gaugeact useful for spectrum from hep-lat/9911003
 *
 *  Tree-level LW with tapole improvement, missing 1x2 in time, also including
 *  2-plaq term. Taken from Morningstar-Peardon, hep-lat/9911003
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/aniso_spectrum_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_aggregate.h"

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
    const bool registered = TheGaugeActFactory::Instance().registerObject(name, 
									  createGaugeAct);
  };


  AnisoSpectrumGaugeActParams::AnisoSpectrumGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try {
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


  // Private initializer
  void
  AnisoSpectrumGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
  {
    // Fold in normalizations and create action
    // NOTE: the 5/3 is folded into beta, hence divided out of c1 and c2
    Real c0 = param.beta;
    plaq = new PlaqGaugeAct(cgs,c0,param.aniso);

    Real c1 = -c0 * (1/(20*param.u_s*param.u_s));
    rect = new RectGaugeAct(cgs,c1);

    Real c2 = -c0;
    plaq_sq = new SpatialTwoPlaqGaugeAct(cgs,param.omega,param.aniso);
  } 

}

