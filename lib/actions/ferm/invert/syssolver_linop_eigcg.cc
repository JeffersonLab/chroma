// $Id: syssolver_linop_eigcg.cc,v 3.1 2009-01-26 22:47:05 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by EigCG
 */

#include "actions/ferm/invert/syssolver_linop_eigcg.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

namespace Chroma
{

  //! Eigenvector accelerated CG system solver namespace
  namespace LinOpSysSolverEigCGEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("EIG_CG_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      Handle< MdagMSystemSolver<LatticeFermion> > mdagmSysSolver(
	TheMdagMFermSystemSolverFactory::Instance().createObject(name,
								 xml_in,
								 path,
								 A));

      return new LinOpSysSolverEigCG<LatticeFermion>(A, mdagmSysSolver);
    }

#if 0
    //! Callback function
    LinOpSystemSolver<LatticeStaggeredFermion>* createFerm(XMLReader& xml_in,
							   const std::string& path,
							   Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      Handle< MdagMSystemSolver<LatticeStaggeredFermion> > mdagmSysSolver(
	TheMdagMStagFermSystemSolverFactory::Instance().createObject(name,
								     xml_in,
								     path,
								     A));

      return new LinOpSysSolverEigCG<LatticeStaggeredFermion>(A, mdagmSysSolver);
    }
#endif

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
//	success &= Chroma::TheLinOpStagFermSystemSolverFactory::Instance().registerObject(name, createStagFerm);
	registered = true;
      }
      return success;
    }
  }


}
