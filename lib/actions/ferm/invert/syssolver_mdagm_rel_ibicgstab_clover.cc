// $Id: syssolver_mdagm_rel_ibicgstab_clover.cc,v 3.1 2009-07-06 19:02:34 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/syssolver_rel_bicgstab_clover_params.h"
#include "actions/ferm/invert/syssolver_mdagm_rel_ibicgstab_clover.h"

namespace Chroma
{
  namespace MdagMSysSolverReliableIBiCGStabCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("RELIABLE_IBICGSTAB_MP_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverReliableIBiCGStabClover(A, state,SysSolverReliableBiCGStabCloverParams(xml_in, path));
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
