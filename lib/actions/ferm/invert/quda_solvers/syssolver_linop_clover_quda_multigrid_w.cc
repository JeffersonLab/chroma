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
  LinOpSysSolverQUDAMULTIGRIDClover::qudaInvert(const CloverTermT<T, U>::Type_t& clover,
				                const CloverTermT<T, U>::Type_t& invclov,
				       		const T& chi_s,
				       		T& psi_s) const{

    SystemSolverResults_t ret;

    void *spinorIn;

    T mod_chi;
    if ( quda_inv_param.matpc_type == QUDA_MATPC_ODD_ODD_ASYMMETRIC ) {
      // asymmetric 
      //
      // Solve A_oo - D A^{-1}_ee D -- chroma conventions.
      // No need to transform source
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
      spinorIn =(void *)&(chi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
      spinorIn = QDPCache::Instance().getDevicePtr( chi_s.getId() );
      QDPIO::cout << "MDAGM spinor in = " << spinorIn << "\n";
#endif
    }
    else if( quda_inv_param.matpc_type == QUDA_MATPC_ODD_ODD) { 
      //
      // symmetric
      // Solve with M_symm = 1 - A^{-1}_oo D A^{-1}ee D 
      //
      // Chroma M =  A_oo ( M_symm )
      //
      //  So  M x = b => A_oo (M_symm) x = b 
      //              =>       M_symm x = A^{-1}_oo b = chi_mod
      invclov.apply(mod_chi, chi_s, PLUS, 1);
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
      spinorIn =(void *)&(mod_chi.elem(rb[1].start()).elem(0).elem(0).real());
#else
      spinorIn = QDPCache::Instance().getDevicePtr( mod_chi.getId() );
      QDPIO::cout << "MDAGM spinor in = " << spinorIn << "\n";
#endif
    }
    else { 
      QDPIO::cout << "MATPC Type not allowed." << std::endl;
      QDPIO::cout << " Allowed are: QUDA_MATPC_ODD_ODD_ASYMMETRIC or QUDA_MATPC_ODD_ODD" << std::endl;
      QDP_abort(1);
    }

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
    void* spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
    void* spinorOut = QDPCache::Instance().getDevicePtr( psi_s.getId() );
    QDPIO::cout << "MDAGM spinor out = " << spinorOut << "\n";
#endif

    // Do the solve here 
    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();
    invertQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param);
    swatch1.stop();


    QDPIO::cout << "Cuda Space Required" << std::endl;
    QDPIO::cout << "\t Spinor:" << quda_inv_param.spinorGiB << " GiB" << std::endl;
    QDPIO::cout << "\t Gauge :" << q_gauge_param.gaugeGiB << " GiB" << std::endl;
    QDPIO::cout << "\t InvClover :" << quda_inv_param.cloverGiB << " GiB" << std::endl;
    QDPIO::cout << "QUDA_MULTIGRID_"<<solver_string<<"_CLOVER_SOLVER: time="<< quda_inv_param.secs <<" s" ;
    QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" ; 
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<std::endl;

    ret.n_count =quda_inv_param.iter;

    return ret;

  }
  

}

