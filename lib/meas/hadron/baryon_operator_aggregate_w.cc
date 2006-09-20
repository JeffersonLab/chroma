// $Id: baryon_operator_aggregate_w.cc,v 1.2 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief All baryon operators
 */

#include "meas/hadron/baryon_operator_aggregate_w.h"

#include "meas/hadron/simple_baryon_operator_w.h"
#include "meas/hadron/group_baryon_operator_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace BaryonOperatorEnv
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
	success &= SimpleBaryonOperatorEnv::registerAll();
	success &= GroupBaryonOperatorEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
