// -*- C++ -*-
// $Id: gauge_monomial.cc,v 3.1 2006-04-19 02:28:38 edwards Exp $
/*! \file
 *  \brief Generic gauge action monomial wrapper
 */

#include "chromabase.h"

#include "update/molecdyn/monomial/gauge_monomial.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "actions/gauge/gaugeacts/gaugeacts_aggregate.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"

namespace Chroma 
{ 
  namespace GaugeMonomialEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >*
    createMonomial(XMLReader& xml, const string& path) 
    {
      QDPIO::cout << "Create monomial: " << name << endl;

      return new GaugeMonomial(GaugeMonomialParams(xml, path));
    }
    
    const std::string name("GAUGE_MONOMIAL");

    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;

      foo &= GaugeActsEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(name, createMonomial);

      return foo;
    }

    //! Register the gaugeact
    const bool registered = registerAll();

  } //end namespace GaugeMonomialEnv



  // Read the parameters
  GaugeMonomialParams::GaugeMonomialParams(XMLReader& xml_in, const string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try {
      // Read the inverter Parameters
      XMLReader xml_tmp(paramtop, "./GaugeAction");
      std::ostringstream os;
      xml_tmp.print(os);
      gauge_act = os.str();
   
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    QDPIO::cout << "GaugeMonomialParams: read \n" << gauge_act << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    GaugeMonomialParams& params) 
  {
    GaugeMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const GaugeMonomialParams& params) 
  {
    // Not implemented
    QDPIO::cerr << GaugeMonomialEnv::name << ": write not implemented" << endl;
    QDP_abort(1);
  }


  // Constructor
  GaugeMonomial::GaugeMonomial(const GaugeMonomialParams& param_) 
  {
    std::istringstream is(param_.gauge_act);
    XMLReader gaugeact_reader(is);

    // Get the name of the gauge act
    std::string gaugeact_string;
    try { 
      read(gaugeact_reader, "/GaugeAction/Name", gaugeact_string);
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Error grepping the gaugeact name: " << e<<  endl;
      QDP_abort(1);
    }

    // Throw an exception if not found
    gaugeact = TheGaugeActFactory::Instance().createObject(gaugeact_string, 
							   gaugeact_reader, 
							   "/GaugeAction");
  }

} //end namespace Chroma


