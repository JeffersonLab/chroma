// $Id: syssolver_polyprec_cg.cc,v 3.2 2006-09-20 20:28:00 edwards Exp $
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

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::ThePolyPrecFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
}
