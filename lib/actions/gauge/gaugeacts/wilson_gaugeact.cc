// $Id: wilson_gaugeact.cc,v 3.1 2006-04-19 02:29:45 edwards Exp $
/*! \file
 *  \brief Wilson gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/wilson_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace WilsonGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
											    const std::string& path) 
    {
      return new WilsonGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
				WilsonGaugeActParams(xml, path));
    }

    const std::string name = "WILSON_GAUGEACT";
    const bool registered = TheGaugeActFactory::Instance().registerObject(name, 
									  createGaugeAct);
  };


  WilsonGaugeActParams::WilsonGaugeActParams(XMLReader& xml_in, const std::string& path) {
    XMLReader paramtop(xml_in, path);

    try {
      read(paramtop, "./beta", beta);

      //  Read optional anisoParam.
      if (paramtop.count("AnisoParam") != 0) 
	read(paramtop, "AnisoParam", aniso);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, const string& path, WilsonGaugeActParams& p) {
    WilsonGaugeActParams tmp(xml, path);
    p=tmp;
  }


  // Private initializer
  void
  WilsonGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
  {
    // Fold in normalizations and create action
    Real coeff = param.beta;
    plaq = new PlaqGaugeAct(cgs,coeff,param.aniso);
  } 

}

