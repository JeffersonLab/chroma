// $Id: syssolver_mdagm_richardson_multiprec_clover.cc,v 3.1 2009-04-17 02:05:32 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/syssolver_richardson_clover_params.h"
#include "actions/ferm/invert/syssolver_mdagm_richardson_multiprec_clover.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverRichardsonCloverEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("RICHARDSON_MP_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }

    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverRichardsonClover(A, state,SysSolverRichardsonCloverParams(xml_in, path));
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
