// $Id: syssolver_mdagm_cg.cc,v 3.4 2009-04-17 02:05:32 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_cg.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMSysSolverCGEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("CG_INVERTER");

      //! Local registration flag
      bool registered = false;
    }


    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverCG<LatticeFermion>(A, SysSolverCGParams(xml_in, path));
    }

    //! Callback function
    MdagMSystemSolver<LatticeFermionF>* createFermF(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > state, 

						  Handle< LinearOperator<LatticeFermionF> > A)
    {
      return new MdagMSysSolverCG<LatticeFermionF>(A, SysSolverCGParams(xml_in, path));
    }

    //! Callback function
    MdagMSystemSolver<LatticeFermionD>* createFermD(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > > state, 

						  Handle< LinearOperator<LatticeFermionD> > A)
    {
      return new MdagMSysSolverCG<LatticeFermionD>(A, SysSolverCGParams(xml_in, path));
    }

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
