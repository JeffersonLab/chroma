// $Id: syssolver_mdagm_cg_lf_clover.cc,v 3.1 2009-05-20 18:22:34 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/syssolver_richardson_clover_params.h"
#include "actions/ferm/invert/syssolver_mdagm_cg_lf_clover.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverCGLFCloverEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("SINGLE_PREC_CG_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }

    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverCGLFClover(A, state,SysSolverCGCloverParams(xml_in, path));
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
