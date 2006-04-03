// $Id: lw_tree_gaugeact.cc,v 3.0 2006-04-03 04:58:54 edwards Exp $
/*! \file
 *  \brief Tree-level tadpole-improved Luscher-Weisz gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/lw_tree_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace LWTreeGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
											    const std::string& path) 
    {
      return new LWTreeGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
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
  LWTreeGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
  {
    // Fold in normalizations and create action
    // NOTE: the 5/3 is folded into beta, hence divided out of c1
    AnisoParam_t aniso;  // empty aniso
    Real c0 = beta;
    plaq = new PlaqGaugeAct(cgs,c0,aniso);

    Real c1 = -c0 * (1/(20*u0*u0));
    rect = new RectGaugeAct(cgs,c1);
  } 

}

