/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/multi_syssolver_mdagm_factory.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/quda_solvers/multi_syssolver_mdagm_cg_clover_quda_w.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "quda.h"

#include <cstdlib>

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMMultiSysSolverCGQudaCloverEnv
  {
    //! Callback function
    MdagMMultiSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						       const std::string& path,
						       Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						       Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMMultiSysSolverCGQudaClover(A, state,MultiSysSolverQUDACloverParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("MULTI_CG_QUDA_CLOVER_INVERTER");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheMdagMFermMultiSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

  SystemSolverResults_t 
  MdagMMultiSysSolverCGQudaClover::qudaInvertMulti(const T& chi_s,
				       multi1d<T>& psi_s,
				       const multi1d<Real> shifts) const{

    SystemSolverResults_t ret;

    void *spinorIn;

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
      spinorIn =(void *)&(chi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
      // have to do this later
#endif

    void** spinorOut = new void*[ shifts.size() ];
    if (spinorOut == nullptr ) { 
      QDPIO::cerr << "Couldn't allocate spinorOut" << std::endl;
      QDP_abort(1);
    }

    if ( shifts.size() > QUDA_MAX_MULTI_SHIFT ) {
       QDPIO::cerr << "You want more shifts than QUDA_MAX_MULTI_SHIFT" << std::endl;
       QDPIO::cerr << "Requested : " << shifts.size() << " QUDA_MAX_MULTI_SHIFT=" << QUDA_MAX_MULTI_SHIFT << std::endl;
       QDP_abort(1);
    }
 
    psi_s.resize( shifts.size());

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
    for(int s=0; s < shifts.size(); s++) {
      //psi_s[s][ rb[1] ] = zero;
      psi_s[s] = zero; // Sanity check 
      spinorOut[s] = (void *)&(psi_s[s].elem(rb[1].start()).elem(0).elem(0).real());
      quda_inv_param.offset[s] = toDouble(shifts[s]);
   } 
#else
    std::vector<int> ids = {chi_s.getId()};
    for(int s=0; s < shifts.size(); s++) {
      // psi_s[s][ rb[1] ] = zero;
      psi_s[s] = zero; // Sanity Check
      ids.push_back( psi_s[s].getId() );
      quda_inv_param.offset[s] = toDouble(shifts[s]);
    }
    auto dev_ptr = GetMemoryPtr( ids );
    spinorIn = dev_ptr[0];
    for(int s=0; s < shifts.size(); s++) {
      spinorOut[s] = dev_ptr[s+1];
    }
#endif

   quda_inv_param.num_offset = shifts.size();

   if( invParam.RsdTarget.size() == 1 ) { 
     for (int i=0; i< quda_inv_param.num_offset; i++) quda_inv_param.tol_offset[i] = toDouble(invParam.RsdTarget[0]);
   }
   else { 
     for (int i=0; i< quda_inv_param.num_offset; i++) quda_inv_param.tol_offset[i] = toDouble(invParam.RsdTarget[i]);
   }

   // Do the solve here 
    StopWatch swatch1; 
    swatch1.reset();
    swatch1.start();
    QDPIO::cout << "CALLING QUDA SOLVER" << std::endl << std::flush ; 
    invertMultiShiftQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param);
     swatch1.stop();

    // Tidy Up
    delete [] spinorOut;


    QDPIO::cout << "QUDA_"<<solver_string<<"_CLOVER_SOLVER: time="<< quda_inv_param.secs <<" s" ;
    QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" ; 
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<std::endl;

    ret.n_count =quda_inv_param.iter;

    return ret;

  }
  

}

