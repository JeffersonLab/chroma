// $Id: syssolver_linop_mr.cc,v 1.2 2009-04-17 02:05:31 bjoo Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by MR
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_mr.h"

namespace Chroma
{

  //! MR system solver namespace
  namespace LinOpSysSolverMREnv
  {
    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverMR<LatticeFermion>(A, SysSolverMRParams(xml_in, path));
    }

    //! Callback function
    LinOpSystemSolver<LatticeStaggeredFermion>* createStagFerm(
      XMLReader& xml_in,
      const std::string& path,
      Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new LinOpSysSolverMR<LatticeStaggeredFermion>(A, SysSolverMRParams(xml_in, path));
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
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	success &= Chroma::TheLinOpStagFermSystemSolverFactory::Instance().registerObject(name, createStagFerm);
	registered = true;
      }
      return success;
    }
  }
}
