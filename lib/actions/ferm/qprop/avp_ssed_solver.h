// -*- C++ -*-
// $Id: avp_ssed_solver.h,v 3.4 2007-01-17 02:39:27 edwards Exp $
/*! \file
 *  \brief DWF/SSE double-prec solver
 */

#ifndef AVP_SSED_SOLVER_H
#define AVP_SSED_SOLVER_H

#include "actions/ferm/qprop/avp_inverter_interface.h"

extern "C" {
  struct MIT_ssed_DWF_Gauge;
  struct MIT_ssed_DWF_Fermion;
};


using namespace QDP;
namespace Chroma 
{ 
  //! AVP's DWF Solver interface
  /*!
   * \ingroup qprop
   *
   * @{
   */
  namespace AVPSolver 
  { 
    
    class SSEDWFSolverD : public AVPSolverInterface< MIT_ssed_DWF_Gauge, MIT_ssed_DWF_Fermion > 
    {
    public:
    protected:
      MIT_ssed_DWF_Fermion* loadFermionRHS(const void* OuterFermion) const; 
      MIT_ssed_DWF_Fermion* loadFermionGuess(const void *OuterFermion) const;
      MIT_ssed_DWF_Fermion* allocateFermion(void) const;
      void saveFermionSolver(void *OuterFermion, 
		       MIT_ssed_DWF_Fermion* CGFermion) const;

      void saveFermionOperator(void *OuterFermion, 
			       MIT_ssed_DWF_Fermion* CGFermion) const;
      void deleteFermion(MIT_ssed_DWF_Fermion* ptr) const;
      int cgInternal(MIT_ssed_DWF_Fermion       *psi,
		     double        *out_eps,
		     int           *out_iter,
		     double        M,
		     double        m_f,
		     const MIT_ssed_DWF_Fermion *x0,
		     const MIT_ssed_DWF_Fermion *eta,
		     double        eps,
		     int           min_iter,
		     int           max_iter)  const;

    public:
     void loadGauge(const void *u,
		    const void *v);
     
     void deleteGauge(void);

     int init(const int lattice[5],
	      void *(*allocator)(size_t size),
	      void (*deallocator)(void *));
     
     // Finalize - destructor call
     void fini(void);
    private:
     MIT_ssed_DWF_Gauge *g;
    };
  }

  /*! @} */   // end of group qprop
}

#endif
