/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_clover_quda_w.h"
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
    namespace LinOpSysSolverQUDACloverEnv
    {

        //! Anonymous namespace
        namespace
        {
            //! Name to be used
            const std::string name("QUDA_CLOVER_INVERTER");

            //! Local registration flag
            bool registered = false;
        }



        LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
                const std::string& path,
                Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 

                Handle< LinearOperator<LatticeFermion> > A)
        {
            return new LinOpSysSolverQUDAClover(A, state,SysSolverQUDACloverParams(xml_in, path));
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

    SystemSolverResults_t LinOpSysSolverQUDAClover::qudaInvert( const T& chi_s, T& psi_s) const {
        SystemSolverResults_t ret;

        void *spinorIn;
        void *spinorOut;

#ifdef BUILD_QUDA_DEVIFACE_SPINOR
        std::vector<QDPCache::ArgKey> ids;
#endif

        // No need to transform source
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
        spinorIn =(void *)&(chi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
        //spinorIn = GetMemoryPtr( chi_s.);
        //QDPIO::cout << "MDAGM spinor in = " << spinorIn << "\n";
        ids.push_back(chi_s.getId());
#endif


#ifndef BUILD_QUDA_DEVIFACE_SPINOR
        spinorOut =(void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());
#else
        ids.push_back(psi_s.getId());
        auto dev_ptr = GetMemoryPtr(ids);
        spinorIn  = dev_ptr[0];
        spinorOut = dev_ptr[1];

#endif
	QDPIO::cout << "b : norm() = ( " << norm2(chi_s, rb[0])<< " , " << norm2( chi_s, rb[1]) << " ) \n";
        // Do the solve here 
        StopWatch swatch1; 
        swatch1.reset();
        swatch1.start();
        invertQuda(spinorOut, spinorIn, (QudaInvertParam*)&quda_inv_param);
        swatch1.stop();

        QDPIO::cout << "QUDA_"<<solver_string<<"_CLOVER_SOLVER: time="<< quda_inv_param.secs <<" s" ;
        QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" ; 
        QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<std::endl;

        ret.n_count =quda_inv_param.iter;

        return ret;

    }


    void LinOpSysSolverQUDAClover::qudaInvertMultiSrc( const std::vector<std::shared_ptr<T>>& psi_s, 
                                                                        const std::vector<std::shared_ptr<const T>>& chi_s,
                                                                        std::vector<SystemSolverResults_t>& res) const {

        std::vector<void *> spinorIn(chi_s.size());
        std::vector<void *> spinorOut(psi_s.size());

        int N_src = chi_s.size();
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
        // Regular non-qdpjit approach. Just collect the pointers 
        for(int soln=0; soln < chi_s.size(); soln++) { 
            spinorIn[soln] = (void *)&(chi_s[soln]->elem(rb[1].start()).elem(0).elem(0).real());
            spinorOut[soln] = (void *)&(psi_s[soln]->elem(rb[1].start()).elem(0).elem(0).real());
        }
#else
       	std::vector<QDPCache::ArgKey> ids(2*N_src);

        for(int soln=0; soln < N_src; soln++) {
          ids[soln]  = chi_s[soln]->getId();
          ids[N_src+soln] = psi_s[soln]->getId();
        }
       
        // Grab all the keys
        auto dev_ptr = QDP_get_global_cache().get_dev_ptrs( multi1d( ids.data(), ids.size()) );

    
        for(int soln=0; soln < N_src; soln++) {
            spinorIn[soln]  = dev_ptr(soln);
            spinorOut[soln] = dev_ptr(N_src+soln);
        }
#endif

        // Local quda_inv_param (?)
        // Relies on quda_inv_param being just a dumb/struct and or copyable
        QudaInvertParam local_quda_inv_param = quda_inv_param ;

        int totalSubgrids=1;
        const multi1d<int>& machine_size=QDP::Layout::logicalSize();

        for (int i = 0; i < Nd; i++) {
            local_quda_inv_param.split_grid[i] = invParam.GridSplitDims[i];
            totalSubgrids *= invParam.GridSplitDims[i];
            if ( machine_size[i] % invParam.GridSplitDims[i] != 0 ) {
                QDPIO::cerr << "The split-grid-subgrid dimensions must divide the number ranks in each dimension exactly\n";
                QDPIO::cerr << "Currently this is not the case: dim=" << i << " machine_size["<<i<<"] = " << machine_size[i]
                            << " and GridSplitDims[" << i << "] = " << invParam.GridSplitDims[i] <<"\n";
                QDPIO::cerr << "Aborting!\n";
                QDP_abort(1);

            }
        }

        if( (chi_s.size() % totalSubgrids) != 0 ) {
            QDPIO::cerr << "The number of split-grid-subgrids must divide the number of sources exactly\n";
            QDPIO::cerr << "Currently it does not: n_src = " << chi_s.size() 
                            << " and split-grid-subgrids = " << totalSubgrids << "\n";
            QDPIO::cerr << "Aborting!\n";
            QDP_abort(1);
                        
        }
        local_quda_inv_param.num_src = chi_s.size();
        local_quda_inv_param.num_src_per_sub_partition = chi_s.size()/totalSubgrids;

        // Do the solve here 
        StopWatch swatch1; 
        swatch1.reset();
        swatch1.start();

        //  gauge is available in the class definition
        //  however we need to do a bit more work to get the clover:
        //
        invertMultiSrcQuda(spinorOut.data(), spinorIn.data(), &local_quda_inv_param);
        swatch1.stop();

        QDPIO::cout << "QUDA_"<<solver_string<<"_CLOVER_SOLVER: time="<< local_quda_inv_param.secs <<" s" ;
        QDPIO::cout << "\tPerformance="<<  local_quda_inv_param.gflops/local_quda_inv_param.secs<<" GFLOPS" ; 
        QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<std::endl;

        for(int soln = 0; soln < chi_s.size(); soln++) res[soln].n_count =local_quda_inv_param.iter;
        return;
    }

}
