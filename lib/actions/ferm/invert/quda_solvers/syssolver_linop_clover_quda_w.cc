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

        std::vector<void *> spinorIn(psi_s.size());
        std::vector<void *> spinorOut(psi_s.size());

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
        // Regular non-qdpjit approach. Just collect the pointers 
        for(int soln=0; soln < chi_s.size(); soln++) { 
            spinorIn[soln] = (void *)&(chi_s[soln]->elem(rb[1].start()).elem(0).elem(0).real());
            spinorOut[soln] = (void *)&(psi_s[soln]->elem(rb[1].start()).elem(0).elem(0).real());
        }
#else
        // QDP-JIT approach: Collect the cache-keys
        std::vector<QDPCache::ArgKey> ids;
        ids.resize(2*chi_s.size());

        for(int soln=0; soln < chi_s.size(); soln++) {
          ids[2*soln]   = chi_s[soln]->getId();
          ids[2*soln+1] = psi_s[soln]->getId();
        }
       
        // Grab all the keys
        auto dev_ptr = GetMemoryPtr(ids);

    
        for(int soln=0; soln < chi_s.size(); soln++) {
            spinorIn[soln]  = dev_ptr[2*soln];
            spinorOut[soln] = dev_ptr[2*soln+1];
        }
#endif

        // Local quda_inv_param (?)
        // Relies on quda_inv_param being just a dumb/struct and or copyable
        QudaInvertParam local_quda_inv_param = quda_inv_param ;

        for (int i = 0; i < 4; i++) local_quda_inv_param.split_grid[i] = 1;
        local_quda_inv_param.num_src = chi_s.size();
        local_quda_inv_param.num_src_per_sub_partition = chi_s.size(); // Since we have only 1 subpartiton

        // Do the solve here 
        StopWatch swatch1; 
        swatch1.reset();
        swatch1.start();

        //  gauge is available in the class definition
        //  however we need to do a bit more work to get the clover:
        //
        void *clover[2];
        void *cloverInv[2];

#ifndef BUILD_QUDA_DEVIFACE_CLOVER
        clover[0] = (void *)&(packed_clov[0]);
        clover[1] = (void *)&(packed_clov[1]);
        cloverInv[0] = (void *)&(packed_invclov[0]);
        cloverInv[1] = (void *)&(packed_invclov[1]);
#else

        // This is a yucky macro and needs the existence of 'clover' and 'cloverInv' to work
        GetMemoryPtrClover(clov->getOffId(),clov->getDiaId(),invclov->getOffId(),invclov->getDiaId());

#endif
        invertMultiSrcCloverQuda(spinorOut.data(), spinorIn.data(), &local_quda_inv_param, (void *)gauge, (QudaGaugeParam*)&q_gauge_param, clover,
                                       cloverInv);
        swatch1.stop();

        QDPIO::cout << "QUDA_"<<solver_string<<"_CLOVER_SOLVER: time="<< local_quda_inv_param.secs <<" s" ;
        QDPIO::cout << "\tPerformance="<<  local_quda_inv_param.gflops/local_quda_inv_param.secs<<" GFLOPS" ; 
        QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<std::endl;

        for(int soln = 0; soln < chi_s.size(); soln++) res[soln].n_count =local_quda_inv_param.iter;
        return;
    }

}
