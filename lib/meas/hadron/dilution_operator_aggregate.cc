// $Id: dilution_operator_aggregate.cc,v 1.1 2007-12-14 06:53:42 edwards Exp $
/*! \file
 *  \brief All dilution operator factories
 */

#include "meas/hadron/dilution_operator_aggregate.h"

#include "meas/hadron/dilution_quark_source_const_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace DilutionOperatorEnv
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
