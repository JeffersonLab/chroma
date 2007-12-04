// $Id: hadron_contract_aggregate.cc,v 3.4 2007-12-04 04:08:53 kostas Exp $
/*! \file
 *  \brief All hadron contraction constructors
 */

#include "meas/hadron/hadron_contract_aggregate.h"

#include "meas/hadron/simple_meson_2pt_w.h"
#include "meas/hadron/delta_2pt_w.h"
//#include "meas/hadron/simple_baryon_2pt_w.h"
//#include "meas/hadron/deriv_meson_2pt_w.h"

#include "meas/hadron/stoch_cond_cont_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace HadronContractEnv
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
	success &= SimpleMeson2PtEnv::registerAll();
	success &= Delta2PtEnv::registerAll();
//	success &= SimpleBaryon2PtEnv::registerAll();
//	success &= DerivMeson2PtEnv::registerAll();
	success &= StochCondContEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
