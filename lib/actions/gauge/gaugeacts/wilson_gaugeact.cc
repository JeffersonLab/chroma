// $Id: wilson_gaugeact.cc,v 2.0 2005-09-25 21:04:31 edwards Exp $
/*! \file
 *  \brief Wilson gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/wilson_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"

namespace Chroma
{
 
  namespace WilsonGaugeActEnv 
  { 
    GaugeAction* createGaugeAct(XMLReader& xml, const std::string& path) 
    {
      return new WilsonGaugeAct(GaugeTypeGaugeBCEnv::reader(xml, path), 
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
  WilsonGaugeAct::init(Handle< GaugeBC > gbc)
  {
    // Fold in normalizations and create action
    Real coeff = param.beta;
    plaq = new PlaqGaugeAct(gbc,coeff,param.aniso);
  } 

}

