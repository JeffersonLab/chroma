// $Id: syssolver_polyprec_cg.cc,v 3.1 2006-07-03 15:26:09 edwards Exp $
/*! \file
 *  \brief Solve a PolyPrec*psi=chi linear system by CG1
 */

#include "actions/ferm/invert/syssolver_polyprec_factory.h"
#include "actions/ferm/invert/syssolver_polyprec_aggregate.h"

#include "actions/ferm/invert/syssolver_polyprec_cg.h"

namespace Chroma
{

  //! CG1 system solver namespace
  namespace PolyPrecSysSolverCGEnv
  {
    //! Callback function
    PolyPrecSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						     const std::string& path,
						     Handle< LinearOperator<LatticeFermion> > A)
    {
      return new PolyPrecSysSolverCG<LatticeFermion>(A, SysSolverCGParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("CG_INVERTER");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::ThePolyPrecFermSystemSolverFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }
}
