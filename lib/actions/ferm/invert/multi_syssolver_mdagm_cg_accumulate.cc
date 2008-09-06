// $Id: multi_syssolver_mdagm_cg_accumulate.cc,v 3.2 2008-09-06 18:35:35 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate_factory.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_mdagm_cg_accumulate.h"

namespace Chroma
{

  //! CG2 system solver namespace
            
  namespace MdagMMultiSysSolverAccumulateCGEnv
  {
    //! Callback function
    MdagMMultiSystemSolverAccumulate<LatticeFermion>* createFerm(XMLReader& xml_in,
						       const std::string& path,
						       Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMMultiSysSolverCGAccumulate<LatticeFermion>(A, SysSolverCGParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("CG_INVERTER");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMdagMFermMultiSystemSolverAccumulateFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
}
