// $Id: syssolver_mdagm_mr.cc,v 1.1 2007-04-11 03:42:07 edwards Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by MR
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_mr.h"

namespace Chroma
{

  //! MR system solver namespace
  namespace MdagMSysSolverMREnv
  {
    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverMR<LatticeFermion>(A, SysSolverMRParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("MR_INVERTER");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
}
