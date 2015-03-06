// $Id: syssolver_mdagm_quda_clover.cc,v 1.6 2009-10-09 13:59:46 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/qphix/syssolver_qphix_clover_params.h"
#include "actions/ferm/invert/qphix/syssolver_mdagm_clover_qphix_w.h"
#include "io/aniso_io.h"

#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"

namespace Chroma
{
  namespace MdagMSysSolverQPhiXCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QPHIX_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }





    // Double precision
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverQPhiXClover<LatticeFermion, LatticeColorMatrix>(A, state,SysSolverQPhiXCloverParams(xml_in, path));
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

