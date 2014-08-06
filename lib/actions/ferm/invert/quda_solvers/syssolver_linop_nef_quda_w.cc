// $Id: syssolver_linop_quda_clover.cc,v 1.6 2009-10-09 13:59:46 bjoo Exp $
/*! \file
*  \brief Solve a MdagM*psi=chi linear system by CG2
*/

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_nef_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_nef_quda_w.h"
#include "io/aniso_io.h"


#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/eoprec_nef_linop_array_w.h"
#include "meas/glue/mesplq.h"
// QUDA Headers
#include <quda.h>
// #include <util_quda.h>

namespace Chroma
{
	namespace LinOpSysSolverQUDANEFEnv
	{

		//! Anonymous namespace
		namespace
		{
			//! Name to be used
			const std::string name("QUDA_NEF_INVERTER");

			//! Local registration flag
			bool registered = false;
		}

		LinOpSystemSolverArray<LatticeFermion>* createFerm(XMLReader& xml_in,	
		const std::string& path,
		Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
		Handle< LinearOperatorArray<LatticeFermion> > A)
		{
			return new LinOpSysSolverQUDANEF(A, state,SysSolverQUDANEFParams(xml_in, path));
		}

		//! Register all the factories
		bool registerAll() 
		{
			bool success = true; 
			if (! registered)
			{
				success &= Chroma::TheLinOpFermSystemSolverArrayFactory::Instance().registerObject(name, createFerm);
				registered = true;
			}
			return success;
		}
	}

	SystemSolverResults_t LinOpSysSolverQUDANEF::qudaInvert(const multi1d<T>& chi_s, multi1d<T>& psi_s) const{

		SystemSolverResults_t ret;
			
		//size of field to copy. it is only half the field, since preconditioning is running
		const multi1d<int>& latdims = Layout::subgridLattSize();
		int halfsize=latdims[0]*latdims[1]*latdims[2]*latdims[3]/2;
		int fermsize=halfsize*Nc*Ns*2;
	
		REAL* spinorIn = new REAL[quda_inv_param.Ls*fermsize];
		REAL* spinorOut = new REAL[quda_inv_param.Ls*fermsize];
		memset(reinterpret_cast<char*>(spinorIn), 0, fermsize*quda_inv_param.Ls*sizeof(REAL));
		memset(reinterpret_cast<char*>(spinorOut), 0, fermsize*quda_inv_param.Ls*sizeof(REAL));
		
		//copy negative parity
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
		for(unsigned int s=0; s<quda_inv_param.Ls; s++){
			memcpy(reinterpret_cast<char*>(&spinorIn[fermsize*s]),reinterpret_cast<char*>(const_cast<double*>(&(chi_s[s].elem(rb[1].start()).elem(0).elem(0).real()))),fermsize*sizeof(REAL));
		}
#else
		//not yet
		//for(unsigned int s=0; s<quda_inv_param.Ls; s++){
		//	spinorIn[s]=QDPCache::Instance().getDevicePtr( chi_s[s].getId() );
		//	spinorOut[s]=QDPCache::Instance().getDevicePtr( psi_s[s].getId() );
		//}
#endif
		
		// Do the solve here 
		StopWatch swatch1; 
		swatch1.reset();
		swatch1.start();
		invertQuda(reinterpret_cast<void*>(spinorOut), reinterpret_cast<void*>(spinorIn), (QudaInvertParam*)&quda_inv_param);
		swatch1.stop();
		
		//copy result
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
		for(unsigned int s=0; s<quda_inv_param.Ls; s++){
			memset(reinterpret_cast<char*>(&(psi_s[s].elem(all.start()).elem(0).elem(0).real())),0,fermsize*2*sizeof(REAL));
			memcpy(reinterpret_cast<char*>(const_cast<double*>(&(psi_s[s].elem(rb[1].start()).elem(0).elem(0).real()))),reinterpret_cast<char*>(&spinorOut[fermsize*s]),fermsize*sizeof(REAL));
		}
#else
		//not yet implemented
		//for(unsigned int s=0; s<quda_inv_param.Ls; s++){
		//	spinorIn[s]=QDPCache::Instance().getDevicePtr( chi_s[s].getId() );
		//	spinorOut[s]=QDPCache::Instance().getDevicePtr( psi_s[s].getId() );
		//}
#endif
		
		//clean up
		delete [] spinorIn;
		delete [] spinorOut;

		QDPIO::cout << "Cuda Space Required" << endl;
		QDPIO::cout << "\t Spinor:" << quda_inv_param.spinorGiB << " GiB" << endl;
		QDPIO::cout << "\t Gauge :" << q_gauge_param.gaugeGiB << " GiB" << endl;
		QDPIO::cout << "QUDA_" << solver_string << "_NEF_SOLVER: time=" << quda_inv_param.secs << " s" ;
		QDPIO::cout << "\tPerformance="<<  quda_inv_param.gflops/quda_inv_param.secs<<" GFLOPS" ; 
		QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() <<" s"<<endl;

		ret.n_count =quda_inv_param.iter;

		return ret;

	}
  

}
