// $Id: lw_tree_gaugeact.cc,v 1.1 2005-01-13 04:30:51 edwards Exp $
/*! \file
 *  \brief Tree-level tadpole-improved Luscher-Weisz gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/lw_tree_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"

namespace Chroma
{
 
  namespace LWTreeGaugeActEnv 
  { 
    GaugeAction* createGaugeAct(XMLReader& xml, const std::string& path) 
    {
      return new LWTreeGaugeAct(GaugeTypeGaugeBCEnv::reader(xml, path), 
				LWTreeGaugeActParams(xml, path));
    }

    const std::string name = "LW_TREE_GAUGEACT";
    const bool registered = TheGaugeActFactory::Instance().registerObject(name, 
									  createGaugeAct);
  };


  LWTreeGaugeActParams::LWTreeGaugeActParams(XMLReader& xml_in, const std::string& path) {
    XMLReader paramtop(xml_in, path);

    try {
      read(paramtop, "./beta", beta);
      read(paramtop, "./u0", u0);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, const string& path, LWTreeGaugeActParams& p) {
    LWTreeGaugeActParams tmp(xml, path);
    p=tmp;
  }


  // Private initializer
  void
  LWTreeGaugeAct::init(Handle< GaugeBC > gbc)
  {
    // Fold in normalizations and create action
    // NOTE: the 5/3 is folded into beta, hence divided out of c1
    Real c0 = beta;
    plaq = new PlaqGaugeAct(gbc,c0);

    Real c1 = -(1/(20*u0*u0))*beta;
    rect = new RectGaugeAct(gbc,c1);
  } 

}

