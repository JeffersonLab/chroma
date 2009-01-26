// $Id: syssolver_mdagm_cg.cc,v 3.3 2009-01-26 22:47:06 edwards Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_cg.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverCGEnv
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
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverCG<LatticeFermion>(A, SysSolverCGParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
}
