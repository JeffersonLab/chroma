// -*- C++ -*-
// $Id: avp_bgld_solver.h,v 3.3 2007-01-17 02:39:27 edwards Exp $
/*! \file
 *  \brief DWF/Bluegene double-prec solver
 */

#ifndef AVP_BGLD_SOLVER_H
#define AVP_BGLD_SOLVER_H

#include "chromabase.h"
#include "actions/ferm/qprop/avp_inverter_interface.h"

extern "C" {
  struct MIT_bluelightd_DWF_Gauge;
  struct MIT_bluelightd_DWF_Fermion;
};


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
    //! DWF double prec-solver for Bluegene
    class BGLDWFSolverD : public AVPSolverInterface< MIT_bluelightd_DWF_Gauge, MIT_bluelightd_DWF_Fermion > 
    {
    public:
    protected:
      MIT_bluelightd_DWF_Fermion* loadFermionRHS(const void* OuterFermion) const; 
      MIT_bluelightd_DWF_Fermion* loadFermionGuess(const void *OuterFermion) const;
      MIT_bluelightd_DWF_Fermion* allocateFermion(void) const;
      void saveFermionSolver(void *OuterFermion, 
			     MIT_bluelightd_DWF_Fermion* CGFermion) const;

      void saveFermionOperator(void *OuterFermion, 
			       MIT_bluelightd_DWF_Fermion* CGFermion) const;
      void deleteFermion(MIT_bluelightd_DWF_Fermion* ptr) const;
      int cgInternal(MIT_bluelightd_DWF_Fermion       *psi,
		     double        *out_eps,
		     int           *out_iter,
		     double        M,
		     double        m_f,
		     const MIT_bluelightd_DWF_Fermion *x0,
		     const MIT_bluelightd_DWF_Fermion *eta,
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
      MIT_bluelightd_DWF_Gauge *g;
    };
  }

  /*! @} */   // end of group qprop
}

#endif
