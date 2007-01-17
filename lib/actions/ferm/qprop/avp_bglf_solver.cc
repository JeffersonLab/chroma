// -*- C++ -*-
/*! \file
 *  \brief DWF/Bluegene single-prec solver
 */

#include "avp_bglf_solver.h"
#include <dwf-bluelightf.h>

namespace Chroma 
{ 
  //! Bluegene single-prec solver
  /*!
   * \ingroup qprop
   *
   * @{
   */
  namespace AVPSolver 
  { 
    MIT_bluelightf_DWF_Fermion* BGLDWFSolverF::loadFermionRHS(const void* OuterFermion) const {
      return MIT_bluelightf_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderRHS);
    }

    MIT_bluelightf_DWF_Fermion* BGLDWFSolverF::loadFermionGuess(const void *OuterFermion) const {
      return MIT_bluelightf_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderGuess);
    }

    MIT_bluelightf_DWF_Fermion* BGLDWFSolverF::allocateFermion(void) const {
      return MIT_bluelightf_DWF_allocate_fermion();
    }

    void BGLDWFSolverF::saveFermionSolver(void *OuterFermion, 
					  MIT_bluelightf_DWF_Fermion* CGFermion) const {
      MIT_bluelightf_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterSolver, CGFermion);
    }

    void BGLDWFSolverF::saveFermionOperator(void *OuterFermion, 
					    MIT_bluelightf_DWF_Fermion* CGFermion) const {
      MIT_bluelightf_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterOperator, CGFermion);
    }
    
    void BGLDWFSolverF::deleteFermion(MIT_bluelightf_DWF_Fermion* ptr) const {
      MIT_bluelightf_DWF_delete_fermion(ptr);
    }

    
    int BGLDWFSolverF::cgInternal(MIT_bluelightf_DWF_Fermion       *psi,
				  double        *out_eps,
				  int           *out_iter,
				  double        M,
				  double        m_f,
				  const MIT_bluelightf_DWF_Fermion *x0,
				  const MIT_bluelightf_DWF_Fermion *eta,
				  double        eps,
				  int           min_iter,
				  int           max_iter)  const 
    {
      QDPIO::cout << "Entering MIT_bluelightf_DWF_cg_solver" << endl;
      return MIT_bluelightf_DWF_cg_solver(psi, out_eps, out_iter, g, M, m_f,
					  x0, eta, eps, min_iter, max_iter);
    }
    
    void BGLDWFSolverF::loadGauge(const void *u,
				  const void *v) { 
      g=MIT_bluelightf_DWF_load_gauge(u, v, NULL, &AVPSolverFunctions::gaugeReader);
    }
     
    void BGLDWFSolverF::deleteGauge(void) {
      MIT_bluelightf_DWF_delete_gauge(g);
    }

     
    // Init the system -- Constructor call?
    int BGLDWFSolverF::init(const int lattice[5],
			    void *(*allocator)(size_t size),
			    void (*deallocator)(void *)) {
      return MIT_bluelightf_DWF_init(lattice, allocator, deallocator);
    }
     
    // Finalize - destructor call
    void BGLDWFSolverF::fini(void) {
      MIT_bluelightf_DWF_fini();
    }
  }

  /*! @} */   // end of group qprop
}


