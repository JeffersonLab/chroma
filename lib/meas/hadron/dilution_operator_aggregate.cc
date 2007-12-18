// $Id: dilution_operator_aggregate.cc,v 1.2 2007-12-18 13:40:25 edwards Exp $
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
	success &= DilutionQuarkSourceConstEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
