// $Id: syssolver_mdagm_ibicgstab.cc,v 3.1 2009-07-02 18:24:52 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_ibicgstab.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverIBiCGStabEnv
  {
    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverIBiCGStab<LatticeFermion>(A, SysSolverBiCGStabParams(xml_in, path));
    }

    MdagMSystemSolver<LatticeFermionF>* createFermF(XMLReader& xml_in,
						    const std::string& path,
						  Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > state, 
						    Handle< LinearOperator<LatticeFermionF> > A)
    {
      return new MdagMSysSolverIBiCGStab<LatticeFermionF>(A, SysSolverBiCGStabParams(xml_in, path));
    }

    MdagMSystemSolver<LatticeFermionD>* createFermD(XMLReader& xml_in,
						    const std::string& path,
						  Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > > state, 

						    Handle< LinearOperator<LatticeFermionD> > A)
    {
      return new MdagMSysSolverIBiCGStab<LatticeFermionD>(A, SysSolverBiCGStabParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("IBICGSTAB_INVERTER");

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
