// $Id: syssolver_linop_quda_clover.cc,v 1.6 2009-10-09 13:59:46 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/intel_solvers/syssolver_intel_clover_params.h"
#include "actions/ferm/invert/intel_solvers/syssolver_linop_clover_intel_w.h"
#include "io/aniso_io.h"

#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"

namespace Chroma
{
  namespace LinOpSysSolverIntelCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("INTEL_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    LinOpSystemSolver<LatticeFermionF>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > state, 
						  
						  Handle< LinearOperator<LatticeFermionF> > A)
    {
      return new LinOpSysSolverIntelClover<LatticeFermionF, LatticeColorMatrixF>(A, state,SysSolverIntelCloverParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

}

