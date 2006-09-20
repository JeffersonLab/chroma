// $Id: lw_1loop_gaugeact.cc,v 3.3 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief 1-loop tadpole-improved Luscher-Weisz gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/lw_1loop_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace LW1LoopGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
											    const std::string& path) 
    {
      return new LW1LoopGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
				 LW1LoopGaugeActParams(xml, path));
    }

    const std::string name = "LW_1LOOP_GAUGEACT";

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


  LW1LoopGaugeActParams::LW1LoopGaugeActParams(XMLReader& xml_in, const std::string& path) {
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

  void read(XMLReader& xml, const string& path, LW1LoopGaugeActParams& p) {
    LW1LoopGaugeActParams tmp(xml, path);
    p=tmp;
  }


  // Private initializer
  void
  LW1LoopGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
  {
    START_CODE();

    // Fold in normalizations and create action
    // NOTE: the 5/3 is folded into beta, hence divided out of c1 and c2
    AnisoParam_t aniso;  // empty aniso
    Real alpha_s = -4.0 * log(u0) / 3.06839;

    Real c0 = beta;
    plaq = new PlaqGaugeAct(cgs,c0,aniso);

    Real c1 = -c0 * (1 + 0.4805*alpha_s) / (20*u0*u0);
    rect = new RectGaugeAct(cgs,c1);

    Real c2 = -c0 * 0.03325 * alpha_s / (u0*u0);
    pg = new PgGaugeAct(cgs,c2);
    
    END_CODE();
  } 

}

