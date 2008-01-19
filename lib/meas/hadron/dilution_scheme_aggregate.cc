// $Id: dilution_scheme_aggregate.cc,v 1.4 2008-01-19 02:10:21 jbulava Exp $
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
