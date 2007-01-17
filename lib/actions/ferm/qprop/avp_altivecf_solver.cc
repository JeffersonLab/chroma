// -*- C++ -*-
/*! \file
 *  \brief DWF/Bluegene altivec solver
 */

#include "avp_altivecf_solver.h"
#include <dwf-altivecf.h>

namespace Chroma 
{ 
  //! Bluegene altivec single-prec solver
  /*!
   * \ingroup qprop
   *
   * @{
   */
  namespace AVPSolver 
  { 
  
    MIT_altivecf_DWF_Fermion* AltiVecDWFSolverF::loadFermionRHS(const void* OuterFermion) const {
      return MIT_altivecf_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderRHS);
    }

    MIT_altivecf_DWF_Fermion* AltiVecDWFSolverF::loadFermionGuess(const void *OuterFermion) const {
      return MIT_altivecf_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderGuess);
    }

    MIT_altivecf_DWF_Fermion* AltiVecDWFSolverF::allocateFermion(void) const {
      return MIT_altivecf_DWF_allocate_fermion();
    }

    void AltiVecDWFSolverF::saveFermionSolver(void *OuterFermion, 
					      MIT_altivecf_DWF_Fermion* CGFermion) const {
      MIT_altivecf_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterSolver, CGFermion);
    }

    void AltiVecDWFSolverF::saveFermionOperator(void *OuterFermion, 
						MIT_altivecf_DWF_Fermion* CGFermion) const {
      MIT_altivecf_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterOperator, CGFermion);
    }
    
    void AltiVecDWFSolverF::deleteFermion(MIT_altivecf_DWF_Fermion* ptr) const {
      MIT_altivecf_DWF_delete_fermion(ptr);
    }

    
    int AltiVecDWFSolverF::cgInternal(MIT_altivecf_DWF_Fermion       *psi,
				      double        *out_eps,
				      int           *out_iter,
				      double        M,
				      double        m_f,
				      const MIT_altivecf_DWF_Fermion *x0,
				      const MIT_altivecf_DWF_Fermion *eta,
				      double        eps,
				      int           min_iter,
				      int           max_iter)  const 
    {
      QDPIO::cout << "Entering MIT_altivecf_DWF_cg_solver" << endl;
      return MIT_altivecf_DWF_cg_solver(psi, out_eps, out_iter, g, M, m_f,
					x0, eta, eps, min_iter, max_iter);
    }
    
    void AltiVecDWFSolverF::loadGauge(const void *u,
				      const void *v) { 
      g=MIT_altivecf_DWF_load_gauge(u, v, NULL, &AVPSolverFunctions::gaugeReader);
    }
     
    void AltiVecDWFSolverF::deleteGauge(void) {
      MIT_altivecf_DWF_delete_gauge(g);
    }

     
    // Init the system -- Constructor call?
    int AltiVecDWFSolverF::init(const int lattice[5],
				void *(*allocator)(size_t size),
				void (*deallocator)(void *)) {
      return MIT_altivecf_DWF_init(lattice, allocator, deallocator);
    }
     
    // Finalize - destructor call
    void AltiVecDWFSolverF::fini(void) {
      MIT_altivecf_DWF_fini();
    }
  }

  /*! @} */   // end of group qprop

}


