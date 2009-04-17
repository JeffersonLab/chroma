// $Id: syssolver_linop_bicgstab.cc,v 3.2 2009-04-17 02:05:31 bjoo Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by MR
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_bicgstab.h"

namespace Chroma
{

  //! MR system solver namespace
  namespace LinOpSysSolverBiCGStabEnv
  {
    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverBiCGStab<LatticeFermion>(A, SysSolverBiCGStabParams(xml_in, path));
    }

    //! Callback function
    LinOpSystemSolver<LatticeFermionF>* createFermF(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > state,
						  Handle< LinearOperator<LatticeFermionF> > A)
    {
      return new LinOpSysSolverBiCGStab<LatticeFermionF>(A, SysSolverBiCGStabParams(xml_in, path));
    }

#if 1
    //! Callback function
    LinOpSystemSolver<LatticeStaggeredFermion>* createStagFerm(
      XMLReader& xml_in,
      const std::string& path,
      Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new LinOpSysSolverBiCGStab<LatticeStaggeredFermion>(A, SysSolverBiCGStabParams(xml_in, path));
    }
#endif

    //! Name to be used
    const std::string name("BICGSTAB_INVERTER");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	success &= Chroma::TheLinOpFFermSystemSolverFactory::Instance().registerObject(name, createFermF);
	//	success &= Chroma::TheLinOpStagFermSystemSolverFactory::Instance().registerObject(name, createStagFerm);
	registered = true;
      }
      return success;
    }
  }
}
