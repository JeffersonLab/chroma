// $Id: syssolver_mdagm_quda_clover.cc,v 1.6 2009-10-09 13:59:46 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/qphix/syssolver_qphix_clover_params.h"
#include "actions/ferm/invert/qphix/syssolver_mdagm_clover_qphix_w.h"
#include "io/aniso_io.h"

#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"

namespace Chroma
{
  namespace MdagMSysSolverQPhiXCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QPHIX_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



#ifndef CHROMA_QPHIX_ARCH_QPX
    // Strictly single precision
    MdagMSystemSolver<LatticeFermionF>* createFermF(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > state, 
						  
						  Handle< LinearOperator<LatticeFermionF> > A)
    {
      return new MdagMSysSolverQPhiXClover<LatticeFermionF, LatticeColorMatrixF>(A, state,SysSolverQPhiXCloverParams(xml_in, path));
    }
#endif

    // Double precision
    MdagMSystemSolver<LatticeFermionD>* createFermD(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > > state, 
						  
						  Handle< LinearOperator<LatticeFermionD> > A)
    {
      return new MdagMSysSolverQPhiXClover<LatticeFermionD, LatticeColorMatrixD>(A, state,SysSolverQPhiXCloverParams(xml_in, path));
    }


    // Double precision
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverQPhiXClover<LatticeFermion, LatticeColorMatrix>(A, state,SysSolverQPhiXCloverParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
#ifndef CHROMA_QPHIX_ARCH_QPX
	success &= Chroma::TheMdagMFermFSystemSolverFactory::Instance().registerObject(name, createFermF);
#endif

	success &= Chroma::TheMdagMFermDSystemSolverFactory::Instance().registerObject(name, createFermD);
	success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

}

