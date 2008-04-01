// $Id: syssolver_linop_OPTeigcg.cc,v 1.2 2008-04-01 04:02:28 kostas Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#include <qdp-lapack.h>

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_OPTeigcg.h"
#include "actions/ferm/invert/inv_eigcg2.h"
#include "actions/ferm/invert/norm_gram_schm.h"


namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverOptEigCGEnv
  {
    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverOptEigCG<LatticeFermion>(A, SysSolverOptEigCGParams(xml_in, path));
    }

#if 0
    //! Callback function
    LinOpSystemSolver<LatticeStaggeredFermion>* createStagFerm(
      XMLReader& xml_in,
      const std::string& path,
      Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new LinOpSysSolverOptEigCG<LatticeStaggeredFermion>(A, SysSolverOptEigCGParams(xml_in, path));
    }
#endif

    //! Name to be used
    const std::string name("EIG_CG_INVERTER");

    //! Local registration flag
    static bool registered = false;

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


  //! Anonymous namespace holding method
  namespace
  {
    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    template<typename T>
    SystemSolverResults_t sysSolver(T& psi, const T& chi, 
				    const LinearOperator<T>& A, 
				    const LinearOperator<T>& MdagM, 
				    const SysSolverOptEigCGParams& invParam)
    {
      START_CODE();

      LinAlg::OptEigInfo& EigInfo = TheNamedObjMap::Instance().getData< LinAlg::OptEigInfo >(invParam.eigen_id);

      Subset s = A.subset() ;

      Complex_C *work=NULL  ;
      Complex_C *V=NULL     ;
      Complex_C *ework=NULL ;
      Complex_C *X ;
      Complex_C *B ;
      
      if(s.hasOrderedRep()){
	X = (Complex_C *) &psi.elem(s.start()).elem().elem().real();
	B = (Complex_C *) &chi.elem(s.start()).elem().elem().real();
      }
      else{//need to copy
	//X = allocate space for them
	//B =  allocate space for them...
	QDPIO::cout<<"OPPS! I have no implemented OPT_EigCG for Linops with non contigius subset\n";
      }
      Complex_C *evecs = (Complex_C *) &EigInfo.evecs[0] ;
      Complex_C *evals = (Complex_C *) &EigInfo.evals[0] ;
      Complex_C *H  = (Complex_C *) &EigInfo.H[0] ;
      Complex_C *HU = (Complex_C *) &EigInfo.HU[0] ;
      
      IncrEigpcg(EigInfo.N, EigInfo.lde, 1, X, B, &EigInfo.ncurEvals, EigInfo.evals.size(), 
		 evecs, evals, H, HU, 
		 MatrixMatvec, NULL, NULL, work, V, ework, esize, 
		 invParam.RsdCG, &invParam.restartTol, invParam.normAestimate, 
		 invParam.updateRestartTol, invParam.MaxCG, invParam.PrintLevel, 
		 invParam.Neig, invParam.Nmax, stdout);
      
      if(!s.hasOrderedRep()){
	QDPIO::cout<<"OPPS! I have no implemented OPT_EigCG for Linops with non contigius subset\n";
      }
      END_CODE();

      return res;
    }

  } // anonymous namespace


  //
  // Wrappers
  //

  // LatticeFermionF
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigCG<LatticeFermionF>::operator()(LatticeFermionF& psi, const LatticeFermionF& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }

#if 0
  //OPT EigCG does not  work with double prec Lattice Fermions
  // LatticeFermionD
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigCG<LatticeFermionD>::operator()(LatticeFermionD& psi, const LatticeFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }


  // Not quite ready yet for these - almost there
  // LatticeStaggeredFermionF
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigCG<LatticeStaggeredFermionF>::operator()(LatticeStaggeredFermionF& psi, const LatticeStaggeredFermionF& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }

  // LatticeStaggeredFermionD
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigCG<LatticeStaggeredFermionD>::operator()(LatticeStaggeredFermionD& psi, const LatticeStaggeredFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }
#endif

}
