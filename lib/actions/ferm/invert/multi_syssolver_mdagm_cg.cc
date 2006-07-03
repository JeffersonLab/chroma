// $Id: multi_syssolver_mdagm_cg.cc,v 3.1 2006-07-03 15:26:08 edwards Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_factory.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_mdagm_cg.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMMultiSysSolverCGEnv
  {
    //! Callback function
    MdagMMultiSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						       const std::string& path,
						       Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMMultiSysSolverCG<LatticeFermion>(A, MultiSysSolverCGParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("CG_INVERTER");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheMdagMFermMultiSystemSolverFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }
}
