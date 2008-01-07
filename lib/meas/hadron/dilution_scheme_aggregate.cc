// $Id: dilution_scheme_aggregate.cc,v 1.1 2008-01-07 15:21:26 jbulava Exp $
/*! \file
 *  \brief All dilution scheme factories
 */

#include "meas/hadron/dilution_scheme_aggregate.h"

#include "meas/hadron/dilution_quark_source_const_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace DilutionSchemeEnv
  {
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Hadron
	success &= QuarkSourceConstDilutionEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
