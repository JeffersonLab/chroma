// -*- C++ -*-
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_bicgstab.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverBiCGStabEnv
  {
    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverBiCGStab<LatticeFermion>(A, SysSolverBiCGStabParams(xml_in, path));
    }

    MdagMSystemSolver<LatticeFermionF>* createFermF(XMLReader& xml_in,
						    const std::string& path,
						  Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > state, 
						    Handle< LinearOperator<LatticeFermionF> > A)
    {
      return new MdagMSysSolverBiCGStab<LatticeFermionF>(A, SysSolverBiCGStabParams(xml_in, path));
    }

    MdagMSystemSolver<LatticeFermionD>* createFermD(XMLReader& xml_in,
						    const std::string& path,
						  Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > > state, 

						    Handle< LinearOperator<LatticeFermionD> > A)
    {
      return new MdagMSysSolverBiCGStab<LatticeFermionD>(A, SysSolverBiCGStabParams(xml_in, path));
    }

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
	success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	success &= Chroma::TheMdagMFermFSystemSolverFactory::Instance().registerObject(name, createFermF);
	success &= Chroma::TheMdagMFermDSystemSolverFactory::Instance().registerObject(name, createFermD);

	registered = true;
      }
      return success;
    }
  }
}
