// $Id: syssolver_mdagm_mr.cc,v 1.2 2009-04-17 02:05:32 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by MR
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_mr.h"

namespace Chroma
{

  //! MR system solver namespace
  namespace MdagMSysSolverMREnv
  {
    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverMR<LatticeFermion>(A, SysSolverMRParams(xml_in, path));
    }

    //! Callback function
    MdagMSystemSolver<LatticeFermionF>* createFermF(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > state, 

						  Handle< LinearOperator<LatticeFermionF> > A)
    {
      return new MdagMSysSolverMR<LatticeFermionF>(A, SysSolverMRParams(xml_in, path));
    }

    //! Callback function
    MdagMSystemSolver<LatticeFermionD>* createFermD(XMLReader& xml_in,
						  const std::string& path
						  Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > > state, 

						  Handle< LinearOperator<LatticeFermionD> > A)
    {
      return new MdagMSysSolverMR<LatticeFermionD>(A, SysSolverMRParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("MR_INVERTER");

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
