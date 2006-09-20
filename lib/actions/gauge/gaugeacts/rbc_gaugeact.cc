// $Id: rbc_gaugeact.cc,v 3.3 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief RG style plaquette + rectangle gauge action following RBC conventions
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/rbc_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace RBCGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
											    const std::string& path) 
    {
      return new RBCGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			     RBCGaugeActParams(xml, path));
    }

    const std::string name = "RBC_GAUGEACT";

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


  RBCGaugeActParams::RBCGaugeActParams(XMLReader& xml_in, const std::string& path) {
    XMLReader paramtop(xml_in, path);

    try {
      read(paramtop, "beta", beta);
      read(paramtop, "c1", c1);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }

  void read(XMLReader& xml, const string& path, RBCGaugeActParams& p) {
    RBCGaugeActParams tmp(xml, path);
    p=tmp;
  }


  // Private initializer
  void
  RBCGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
  {
    START_CODE();

    // Fold in normalizations and create action
    // Here, the rectangle weight is relative to the plaquette
    AnisoParam_t aniso;  // empty aniso
    Real coeff0 = beta*(1 - 8*c1);
    plaq = new PlaqGaugeAct(cgs,coeff0,aniso);

    Real coeff1 = beta*c1;
    rect = new RectGaugeAct(cgs,coeff1);
    
    END_CODE();
  } 

}

