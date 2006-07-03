// $Id: syssolver_linop_cg.cc,v 3.1 2006-07-03 15:26:08 edwards Exp $
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
    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverCG<LatticeFermion>(A, SysSolverCGParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("CG_INVERTER");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }
}
