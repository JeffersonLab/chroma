// $Id: multi_syssolver_linop_mr.cc,v 1.1 2007-04-11 03:41:36 edwards Exp $
/*! \file
 *  \brief Solve a (M+shift)*psi=chi linear system by MR
 */

#include "actions/ferm/invert/multi_syssolver_linop_factory.h"
#include "actions/ferm/invert/multi_syssolver_linop_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_linop_mr.h"

namespace Chroma
{

  //! MR system solver namespace
  namespace LinOpMultiSysSolverMREnv
  {
    //! Callback function
    LinOpMultiSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						       const std::string& path,
						       Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpMultiSysSolverMR<LatticeFermion>(A, MultiSysSolverMRParams(xml_in, path));
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
	success &= Chroma::TheLinOpFermMultiSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
}
