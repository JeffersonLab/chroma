#include "avp_ssef_solver.h"

using namespace QDP;
namespace Chroma 
{ 
  namespace AVPSolver 
  { 
    
#include <dwf-ssef.h>

    MIT_ssef_DWF_Fermion* SSEDWFSolverF::loadFermionRHS(const void* OuterFermion) const {
      return MIT_ssef_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderRHS);
    }

    MIT_ssef_DWF_Fermion* SSEDWFSolverF::loadFermionGuess(const void *OuterFermion) const {
      return MIT_ssef_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderGuess);
    }

    MIT_ssef_DWF_Fermion* SSEDWFSolverF::allocateFermion(void) const {
      return MIT_ssef_DWF_allocate_fermion();
    }

    void SSEDWFSolverF::saveFermionSolver(void *OuterFermion, 
					  MIT_ssef_DWF_Fermion* CGFermion) const {
      MIT_ssef_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterSolver, CGFermion);
    }

    void SSEDWFSolverF::saveFermionOperator(void *OuterFermion, 
					    MIT_ssef_DWF_Fermion* CGFermion) const {
      MIT_ssef_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterOperator, CGFermion);
    }
    
    void SSEDWFSolverF::deleteFermion(MIT_ssef_DWF_Fermion* ptr) const {
      MIT_ssef_DWF_delete_fermion(ptr);
    }

    
    int SSEDWFSolverF::cgInternal(MIT_ssef_DWF_Fermion       *psi,
				  double        *out_eps,
				  int           *out_iter,
				  double        M,
				  double        m_f,
				  const MIT_ssef_DWF_Fermion *x0,
				  const MIT_ssef_DWF_Fermion *eta,
				  double        eps,
				  int           min_iter,
				  int           max_iter)  const 
    {
      QDPIO::cout << "Entering MIT_ssef_DWF_cg_solver" << endl;
      return MIT_ssef_DWF_cg_solver(psi, out_eps, out_iter, g, M, m_f,
				    x0, eta, eps, min_iter, max_iter);
    }
    
    void SSEDWFSolverF::loadGauge(const void *u,
				  const void *v) { 
      g=MIT_ssef_DWF_load_gauge(u, v, NULL, &AVPSolverFunctions::gaugeReader);
    }
     
    void SSEDWFSolverF::deleteGauge(void) {
      MIT_ssef_DWF_delete_gauge(g);
    }

     
    // Init the system -- Constructor call?
    int SSEDWFSolverF::init(const int lattice[5],
			    void *(*allocator)(size_t size),
			    void (*deallocator)(void *)) {
      return MIT_ssef_DWF_init(lattice, allocator, deallocator);
    }
     
    // Finalize - destructor call
    void SSEDWFSolverF::fini(void) {
      MIT_ssef_DWF_fini();
    }
  }
}


