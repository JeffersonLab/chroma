#ifndef AVP_SSE_SOLVER_H
#define AVP_SSE_SOLVER_H

#include "avp_inverter_interface.h"

extern "C" { 
  struct MIT_ssef_DWF_Gauge;
  struct MIT_ssef_DWF_Fermion;
  struct MIT_ssed_DWF_Gauge;
  struct MIT_ssed_DWF_Fermion;
};

using namespace QDP;
namespace Chroma { 
  namespace AVPSolver { 
    
#include <dwf-ssef.h>

    class SSEDWFSolverF : public AVPSolverInterface< MIT_ssef_DWF_Gauge, MIT_ssef_DWF_Fermion > {
public:
    protected:
      Fermion* loadFermionRHS(const void* OuterFermion) const {
	return MIT_ssef_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderRHS);
      }

      Fermion* loadFermionGuess(const void *OuterFermion) const {
	return MIT_ssef_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderGuess);
      }

      Fermion* allocateFermion(void) const {
	return MIT_ssef_DWF_allocate_fermion();
      }

      void saveFermionSolver(void *OuterFermion, 
		       Fermion* CGFermion) const {
	MIT_ssef_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterSolver, CGFermion);
      }

     void saveFermionOperator(void *OuterFermion, 
		       Fermion* CGFermion) const {
	MIT_ssef_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterOperator, CGFermion);
      }

     void deleteFermion(Fermion* ptr) const {
       MIT_ssef_DWF_delete_fermion(ptr);
     }


     int cgInternal(Fermion       *psi,
		    double        *out_eps,
		    int           *out_iter,
		    double        M,
		    double        m_f,
		    const Fermion *x0,
		    const Fermion *eta,
		    double        eps,
		    int           min_iter,
		    int           max_iter)  const {
       return MIT_ssef_DWF_cg_solver(psi, out_eps, out_iter, g, M, m_f,
				     x0, eta, eps, min_iter, max_iter);
     }

    public:
     

     void loadGauge(const void *u,
		      const void *v) { 
       g=MIT_ssef_DWF_load_gauge(u, v, NULL, &AVPSolverFunctions::gaugeReader);
     }
     
     void deleteGauge(void) {
       MIT_ssef_DWF_delete_gauge(g);
     }

     
     // Init the system -- Constructor call?
     int init(const int lattice[5],
	       void *(*allocator)(size_t size),
	       void (*deallocator)(void *)) {
       return MIT_ssef_DWF_init(lattice, allocator, deallocator);
     }
     
     // Finalize - destructor call
     void fini(void) {
       MIT_ssef_DWF_fini();
     }
     
    private:
     Gauge *g;
    };
  


#include <dwf-ssed.h>    
    class SSEDWFSolverD : public AVPSolverInterface< MIT_ssed_DWF_Gauge, MIT_ssed_DWF_Fermion > {
public:
    protected:
      Fermion* loadFermionRHS(const void* OuterFermion) const {
	return MIT_ssed_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderRHS);
      }

      Fermion* loadFermionGuess(const void *OuterFermion) const {
	return MIT_ssed_DWF_load_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionReaderGuess);
      }

      Fermion* allocateFermion(void) const {
	return MIT_ssed_DWF_allocate_fermion();
      }

      void saveFermionSolver(void *OuterFermion, 
		       Fermion* CGFermion) const {
	MIT_ssed_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterSolver, CGFermion);
      }

     void saveFermionOperator(void *OuterFermion, 
		       Fermion* CGFermion) const {
	MIT_ssed_DWF_save_fermion(OuterFermion, NULL, &AVPSolverFunctions::fermionWriterOperator, CGFermion);
      }

     void deleteFermion(Fermion* ptr) const {
       MIT_ssed_DWF_delete_fermion(ptr);
     }


     int cgInternal(Fermion       *psi,
		    double        *out_eps,
		    int           *out_iter,
		    double        M,
		    double        m_f,
		    const Fermion *x0,
		    const Fermion *eta,
		    double        eps,
		    int           min_iter,
		    int           max_iter)  const {
       return MIT_ssed_DWF_cg_solver(psi, out_eps, out_iter, g, M, m_f,
				     x0, eta, eps, min_iter, max_iter);
     }

    public:
     

     void loadGauge(const void *u,
		      const void *v) { 
       g=MIT_ssed_DWF_load_gauge(u, v, NULL, &AVPSolverFunctions::gaugeReader);
     }
     
     void deleteGauge(void) {
       MIT_ssed_DWF_delete_gauge(g);
     }

     
     // Init the system -- Constructor call?
     int init(const int lattice[5],
	       void *(*allocator)(size_t size),
	       void (*deallocator)(void *)) {
       return MIT_ssed_DWF_init(lattice, allocator, deallocator);
     }
     
     // Finalize - destructor call
     void fini(void) {
       MIT_ssed_DWF_fini();
     }
     
    private:
     Gauge *g;
    };
  };




};

#endif
