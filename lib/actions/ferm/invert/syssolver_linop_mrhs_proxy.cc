/*! \file
 *  \brief Solve a M*psi=chi linear system
 */
#include "chromabase.h"
#include "state.h"
#include "handle.h"
#include "actions/ferm/invert/syssolver_linop_mrhs_factory.h"
#include "actions/ferm/invert/syssolver_mrhs_proxy_params.h"
#include "actions/ferm/invert/syssolver_mrhs_proxy.h"

using namespace QDP;
using namespace Chroma;
namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverMRHSProxyEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("MULTI_RHS_PROXY_INVERTER");

      //! Local registration flag
      bool registered = false;
    }


    //! Callback function
    LinOpMRHSSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  const FermAct4D< LatticeFermion,
						  	  	  	  	  	  	 multi1d<LatticeColorMatrix>,
												 multi1d<LatticeColorMatrix> >& S_ferm,
						  Handle< FermState<
						                     LatticeFermion, 
						                     multi1d<LatticeColorMatrix>,
						                     multi1d<LatticeColorMatrix> 
					 	  	  	  	  	   > > state)
    {
      return new LinOpMRHSSysSolverProxy<LatticeFermion,multi1d<LatticeColorMatrix>,
              multi1d<LatticeColorMatrix>>(SysSolverMRHSProxyParams(xml_in,path),S_ferm,state);

    }

#if 0
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
      return new LinOpSysSolverCG<LatticeFermionF>(A, SysSolverCGParams(xml_in, path));
    }

    //! Callback function
    LinOpSystemSolver<LatticeStaggeredFermion>* createStagFerm(XMLReader& xml_in,
							       const std::string& path,
							       Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new LinOpSysSolverCG<LatticeStaggeredFermion>(A, SysSolverCGParams(xml_in, path));
    }
#endif

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFermMRHSSystemSolverFactory::Instance().registerObject(name, createFerm);
#if 0
	success &= Chroma::TheLinOpFFermSystemSolverFactory::Instance().registerObject(name, createFermF);
	success &= Chroma::TheLinOpStagFermSystemSolverFactory::Instance().registerObject(name, createStagFerm);
#endif
	registered = true;
      }
      return success;
    }
  }
}
