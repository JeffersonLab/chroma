// $Id: syssolver_linop_cg_array.cc,v 3.1 2006-07-03 15:26:08 edwards Exp $
/*! \file
 *  \brief Solve a CG1 system
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_cg_array.h"

namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverCGArrayEnv
  {
    //! Callback function
    LinOpSystemSolverArray<LatticeFermion>* createFerm(XMLReader& xml_in,
						       const std::string& path,
						       Handle< LinearOperatorArray<LatticeFermion> > A)
    {
      return new LinOpSysSolverCGArray<LatticeFermion>(A, SysSolverCGParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("CG_INVERTER");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheLinOpFermSystemSolverArrayFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }
}
