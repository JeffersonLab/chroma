// $Id: dilution_scheme_aggregate.cc,v 1.5 2008-04-21 03:19:35 edwards Exp $
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
	success &= DilutionQuarkSourceConstEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
