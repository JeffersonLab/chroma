/*! \file
 *  \brief Solve a CG1 system
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_cg_array.h"

namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverCGArrayEnv
  {
    //! Callback function
    LinOpSystemSolverArray<LatticeFermion>* createFerm(XMLReader& xml_in,
						       const std::string& path,
						       Handle< FermState<
					                                 LatticeFermion, 
						                         multi1d<LatticeColorMatrix>,
						                         multi1d<LatticeColorMatrix> 
						             > 
							  > state, 
						       Handle< LinearOperatorArray<LatticeFermion> > A)
    {
      return new LinOpSysSolverCGArray<LatticeFermion>(A, SysSolverCGParams(xml_in, path));
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
	success &= Chroma::TheLinOpFermSystemSolverArrayFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
}
