// $Id: syssolver_mdagm_bicgstab.cc,v 3.1 2007-05-01 14:39:13 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_bicgstab.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverBiCGStabEnv
  {
    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverBiCGStab<LatticeFermion>(A, SysSolverBiCGStabParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("BICGSTAB_INVERTER");

    //! Local registration flag
    static bool registered = false;

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
