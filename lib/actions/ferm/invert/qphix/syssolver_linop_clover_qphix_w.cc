/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/qphix/syssolver_qphix_clover_params.h"
#include "actions/ferm/invert/qphix/syssolver_linop_clover_qphix_w.h"
#include "io/aniso_io.h"

#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"

namespace Chroma
{
  namespace LinOpSysSolverQPhiXCloverEnv
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
    LinOpSystemSolver<LatticeFermionF>* createFermF(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermionF, multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> > > state, 
						  
						  Handle< LinearOperator<LatticeFermionF> > A)
    {
      return new LinOpSysSolverQPhiXClover<LatticeFermionF, LatticeColorMatrixF>(A, state,SysSolverQPhiXCloverParams(xml_in, path));
    }
#endif

    // Double precision
    LinOpSystemSolver<LatticeFermionD>* createFermD(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermionD, multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> > > state, 
						  
						  Handle< LinearOperator<LatticeFermionD> > A)
    {
      return new LinOpSysSolverQPhiXClover<LatticeFermionD, LatticeColorMatrixD>(A, state,SysSolverQPhiXCloverParams(xml_in, path));
    }


    // Double precision
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverQPhiXClover<LatticeFermion, LatticeColorMatrix>(A, state,SysSolverQPhiXCloverParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
#ifndef CHROMA_QPHIX_ARCH_QPX
	success &= Chroma::TheLinOpFFermSystemSolverFactory::Instance().registerObject(name, createFermF);
#endif

	success &= Chroma::TheLinOpDFermSystemSolverFactory::Instance().registerObject(name, createFermD);
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

}

