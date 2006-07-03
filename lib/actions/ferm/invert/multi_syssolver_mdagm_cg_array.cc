// $Id: multi_syssolver_mdagm_cg_array.cc,v 3.1 2006-07-03 15:26:08 edwards Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by multi-shift CG
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_factory.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/multi_syssolver_mdagm_cg_array.h"

namespace Chroma
{

  //! CG system solver namespace
  namespace MdagMMultiSysSolverCGArrayEnv
  {
    //! Callback function
    MdagMMultiSystemSolverArray<LatticeFermion>* createFerm(XMLReader& xml_in,
							    const std::string& path,
							    Handle< LinearOperatorArray<LatticeFermion> > A)
    {
      return new MdagMMultiSysSolverCGArray<LatticeFermion>(A, MultiSysSolverCGParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("CG_INVERTER");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheMdagMFermMultiSystemSolverArrayFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }
}
