/*! \file
 *  \QUDA MULTIGRID Wilson solver.
 */
// comment
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_wilson_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_wilson_quda_multigrid_w.h"
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
  namespace LinOpSysSolverQUDAMULTIGRIDWilsonEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QUDA_MULTIGRID_WILSON_INVERTER");

      //! Local registration flag
      bool registered = false;
    }



    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverQUDAMULTIGRIDWilson(A, state,SysSolverQUDAMULTIGRIDWilsonParams(xml_in, path));
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
  LinOpSysSolverQUDAMULTIGRIDWilson::qudaInvert(const T& chi_s,
				       		T& psi_s) const{

    SystemSolverResults_t ret;
   

    T mod_chi;

    // Off parity part may contain junk.
    // Not a problem with BiCGStab etc, where we use an operator defined only on
    // one parity, but in MG we will use an operator defined on both parities.
    // So I will explicitly zero the off parity, into a mod_chi to use as a source.

    mod_chi[rb[0]] = zero; // Zero odd parity
    mod_chi[rb[1]] = chi_s;

    if( invParam.asymmetricP  ) {

      // In this case Chroma and QUDA operators disagree by A_oo = (Mass Term)
      // (or clover term  in the case of clover). SO I need to rescale mod_chi
      
      Real diag_mass;
      {
	// auto is C++11 so I don't have to remember all the silly typenames
	auto wlparams = invParam.WilsonParams;
	
	auto aniso = wlparams.anisoParam;
	
	Real ff = where(aniso.anisoP, Real(1) / aniso.xi_0, Real(1));
	diag_mass = 1 + (Nd-1)*ff + wlparams.Mass;
      }

      mod_chi[rb[1]] /= diag_mass;
    }

    // Set pointers for QUDA
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
    void *spinorIn =(void *)&(mod_chi.elem(rb[1].start()).elem(0).elem(0).real());
    void* spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
    // void*  spinorIn = GetMemoryPtr( mod_chi.getId() );
    // QDPIO::cout << "QUDA_MULTIGRID_QDPJIT spinor in = " << spinorIn << "\n";
    // void* spinorOut = GetMemoryPtr( psi_s.getId() );
    // QDPIO::cout << "QUDA_MULTIGRID_QDPJIT spinor out  = " << spinorOut << "\n";
    void* spinorIn;
    void* spinorOut;
    GetMemoryPtr2(spinorIn,spinorOut,mod_chi.getId(),psi_s.getId());
#endif

    // Do the solve here 
    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();
    invertQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param);
    swatch1.stop();

   
    QDPIO::cout << "QUDA_MULTIGRID_"<<solver_string<<"_WILSON_SOLVER: time="<< quda_inv_param.secs <<" s" ;
    QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" ; 
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<std::endl;

    ret.n_count =quda_inv_param.iter;

    return ret;

  }
  

}

