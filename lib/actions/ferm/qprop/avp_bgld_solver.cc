/*! \file
 *  \brief DWF/Bluegene double-prec solver
 */

#include "actions/ferm/qprop/avp_bgld_solver.h"
   
#include <dwf-bluelightd.h>


using namespace QDP;
namespace Chroma 
{ 
  //! Bluegene double-prec solver
  /*!
   * \ingroup qprop
   *
   * @{
   */
  namespace AVPSolver 
  { 
    MIT_bluelightd_DWF_Fermion* BGLDWFSolverD::loadFermionRHS(const void* OuterFermion) const {
      return MIT_bluelightd_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderRHS);
    }

    MIT_bluelightd_DWF_Fermion* BGLDWFSolverD::loadFermionGuess(const void *OuterFermion) const {
	return MIT_bluelightd_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderGuess);
    }

    MIT_bluelightd_DWF_Fermion* BGLDWFSolverD::allocateFermion(void) const {
      return MIT_bluelightd_DWF_allocate_fermion();
    }

    void BGLDWFSolverD::saveFermionSolver(void *OuterFermion, 
					  MIT_bluelightd_DWF_Fermion* CGFermion) const {
      MIT_bluelightd_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterSolver, CGFermion);
    }

    void BGLDWFSolverD::saveFermionOperator(void *OuterFermion, 
					    MIT_bluelightd_DWF_Fermion* CGFermion) const {
      MIT_bluelightd_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterOperator, CGFermion);
    }
    
    void BGLDWFSolverD::deleteFermion(MIT_bluelightd_DWF_Fermion* ptr) const {
      MIT_bluelightd_DWF_delete_fermion(ptr);
    }

    
    int BGLDWFSolverD::cgInternal(MIT_bluelightd_DWF_Fermion       *psi,
				  double        *out_eps,
				  int           *out_iter,
				  double        M,
				  double        m_f,
				  const MIT_bluelightd_DWF_Fermion *x0,
				  const MIT_bluelightd_DWF_Fermion *eta,
				  double        eps,
				  int           min_iter,
				  int           max_iter)  const 
    {
      QDPIO::cout << "Entering MIT_bluelightd_DWF_cg_solver" << endl;
      return MIT_bluelightd_DWF_cg_solver(psi, out_eps, out_iter, g, M, m_f,
				    x0, eta, eps, min_iter, max_iter);
    }
    
    void BGLDWFSolverD::loadGauge(const void *u,
				  const void *v) { 
      g=MIT_bluelightd_DWF_load_gauge(u, v, NULL, &AVPSolverFunctions::gaugeReader);
    }
     
    void BGLDWFSolverD::deleteGauge(void) {
      MIT_bluelightd_DWF_delete_gauge(g);
    }

     
    // Init the system -- Constructor call?
    int BGLDWFSolverD::init(const int lattice[5],
			    void *(*allocator)(size_t size),
			    void (*deallocator)(void *)) {
      return MIT_bluelightd_DWF_init(lattice, allocator, deallocator);
    }
     
    // Finalize - destructor call
    void BGLDWFSolverD::fini(void) {
      MIT_bluelightd_DWF_fini();
    }
  }

  /*! @} */   // end of group qprop

}
