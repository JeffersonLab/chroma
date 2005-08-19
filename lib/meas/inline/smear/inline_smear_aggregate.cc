// $Id: inline_smear_aggregate.cc,v 1.2 2005-08-19 05:32:12 edwards Exp $
/*! \file
 *  \brief Inline smear measurement aggregator
 */

#include "meas/inline/smear/inline_smear_aggregate.h"
#include "meas/inline/smear/inline_ape_smear.h"
#include "meas/inline/smear/inline_hyp_smear.h"
#include "meas/inline/smear/inline_hyp_smear4d.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineSmearAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 
      success &= InlineAPESmearEnv::registered;
      success &= InlineHypSmearEnv::registered;
      success &= InlineHypSmear4dEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
