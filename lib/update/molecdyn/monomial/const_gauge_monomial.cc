// -*- C++ -*-
// $Id: gauge_monomial.cc,v 3.2 2006-09-20 20:28:05 edwards Exp $
/*! \file
 *  \brief Generic gauge action monomial wrapper
 */

#include "chromabase.h"

#include "update/molecdyn/monomial/const_gauge_monomial.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "actions/gauge/gaugeacts/gaugeacts_aggregate.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"

namespace Chroma 
{ 
  namespace ConstGaugeMonomialEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >*
      createMonomial(XMLReader& xml, const string& path) 
      {
	QDPIO::cout << "Create monomial: " << name << endl;

	return new ConstGaugeMonomial(GaugeMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name("CONST_GAUGE_MONOMIAL");
    
    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= GaugeActsEnv::registerAll();
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  } //end namespace ConstGaugeMonomialEnv

} //end namespace Chroma


