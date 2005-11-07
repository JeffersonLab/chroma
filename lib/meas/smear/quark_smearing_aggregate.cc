// $Id: quark_smearing_aggregate.cc,v 2.1 2005-11-07 06:40:55 edwards Exp $
/*! \file
 *  \brief All quark smearing
 */

#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/gaus_quark_smearing.h"

namespace Chroma
{

  // Registration aggregator
  namespace PropSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Propagator smearing
      success &= GausPropSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

  // Registration aggregator
  namespace FermSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Wilson-like Fermion smearing
      success &= GausFermSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

  // Registration aggregator
  namespace ColorVecSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Color-vector
      success &= GausColorVecSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
