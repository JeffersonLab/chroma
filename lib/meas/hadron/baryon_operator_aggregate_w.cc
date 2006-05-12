// $Id: baryon_operator_aggregate_w.cc,v 1.1 2006-05-12 03:38:01 edwards Exp $
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
    bool registerAll() 
    {
      bool success = true;

      // Hadron
      success &= SimpleBaryonOperatorEnv::registered;
      success &= GroupBaryonOperatorEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
