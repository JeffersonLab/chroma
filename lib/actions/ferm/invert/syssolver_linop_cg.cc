// $Id: syssolver_linop_cg.cc,v 3.4 2009-01-26 22:47:05 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_cg.h"

namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverCGEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("CG_INVERTER");

      //! Local registration flag
      bool registered = false;
    }


    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverCG<LatticeFermion>(A, SysSolverCGParams(xml_in, path));
    }

    //! Callback function
    LinOpSystemSolver<LatticeStaggeredFermion>* createStagFerm(XMLReader& xml_in,
							       const std::string& path,
							       Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new LinOpSysSolverCG<LatticeStaggeredFermion>(A, SysSolverCGParams(xml_in, path));
    }

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
