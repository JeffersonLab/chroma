/*! \file
 *  \brief Solve a M*psi=chi linear system by Fake
 */
#include "state.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_fake.h"

namespace Chroma
{

  //! Fake system solver namespace
  namespace LinOpSysSolverFakeEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("FAKE_INVERTER");

      //! Local registration flag
      bool registered = false;
    }


    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState<
						                     LatticeFermion, 
						                     multi1d<LatticeColorMatrix>,
						                     multi1d<LatticeColorMatrix> 
					 	  > 
							  > state, 

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverFake<LatticeFermion>(A);
    }

    //! Callback function
    LinOpSystemSolver<LatticeFermionF>* createFermF(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState<
						                     LatticeFermionF, 
						                     multi1d<LatticeColorMatrixF>,
						                     multi1d<LatticeColorMatrixF> 
						  > 
							  > state, 

						  Handle< LinearOperator<LatticeFermionF> > A)
    {
      return new LinOpSysSolverFake<LatticeFermionF>(A);
    }

    //! Callback function
    LinOpSystemSolver<LatticeStaggeredFermion>* createStagFerm(XMLReader& xml_in,
							       const std::string& path,
							       Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new LinOpSysSolverFake<LatticeStaggeredFermion>(A);
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	success &= Chroma::TheLinOpFFermSystemSolverFactory::Instance().registerObject(name, createFermF);
	success &= Chroma::TheLinOpStagFermSystemSolverFactory::Instance().registerObject(name, createStagFerm);
	registered = true;
      }
      return success;
    }
  }
}
