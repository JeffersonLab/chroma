// $Id: syssolver_linop_cg_timing.cc,v 3.2 2009-04-17 02:05:31 bjoo Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_cg_timing.h"

namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverCGTimingEnv
  {
    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverCGTiming<LatticeFermion>(A, SysSolverCGParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("CG_INVERTER_TIMINGS");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
}
