/*! \file
 *  \QUDA MULTIGRID Clover solver.
 */
// comment
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_clover_quda_multigrid_w.h"
#include "io/aniso_io.h"


#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"
// QUDA Headers
#include <quda.h>
// #include <util_quda.h>
#include "actions/ferm/invert/quda_solvers/quda_mg_utils.h"

namespace Chroma
{
  namespace LinOpSysSolverQUDAMULTIGRIDCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QUDA_MULTIGRID_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverQUDAMULTIGRIDClover(A, state,SysSolverQUDAMULTIGRIDCloverParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

  SystemSolverResults_t 
  LinOpSysSolverQUDAMULTIGRIDClover::qudaInvert(const CloverTermT<T, U>& clover,
				                const CloverTermT<T, U>& invclov,
				       		const T& chi_s,
				       		T& psi_s) const{

    SystemSolverResults_t ret;

    T mod_chi;

    // Copy source into mod_chi, and zero the off-parity
    mod_chi[rb[0]] = zero;
  

    if( invParam.asymmetricP ) { 
      //
      // symmetric
      // Solve with M_symm = 1 - A^{-1}_oo D A^{-1}ee D 
      //
      // Chroma M =  A_oo ( M_symm )
      //
      //  So  M x = b => A_oo (M_symm) x = b 
      //              =>       M_symm x = A^{-1}_oo b = chi_mod
      invclov.apply(mod_chi, chi_s, PLUS, 1);
    }
    else {
      mod_chi[rb[1]] = chi_s;
    }

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
    void* spinorIn =(void *)&(mod_chi.elem(rb[1].start()).elem(0).elem(0).real());
    void* spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
    // void* spinorIn = GetMemoryPtr( mod_chi.getId() );
    // void* spinorOut = GetMemoryPtr( psi_s.getId() );
    void* spinorIn;
    void* spinorOut;
    GetMemoryPtr2(spinorIn,spinorOut,mod_chi.getId(),psi_s.getId())

#endif

    // Do the solve here 
    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();
    invertQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param);
    swatch1.stop();


    QDPIO::cout << solver_string<< "time="<< quda_inv_param.secs <<" s" ;
    QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" ; 
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<std::endl;

    ret.n_count =quda_inv_param.iter;
    ret.resid = quda_inv_param.true_res;
    return ret;

  }
  

}

