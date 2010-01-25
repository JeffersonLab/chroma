// $Id: multi_syssolver_mdagm_cg_array.cc,v 3.2 2006-09-20 20:28:00 edwards Exp $
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
							    Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >,
							    Handle< LinearOperatorArray<LatticeFermion> > A)
    {
      return new MdagMMultiSysSolverCGArray<LatticeFermion>(A, MultiSysSolverCGParams(xml_in, path));
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
	success &= Chroma::TheMdagMFermMultiSystemSolverArrayFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
}
