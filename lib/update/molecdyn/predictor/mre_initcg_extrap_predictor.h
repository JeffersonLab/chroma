// -*- C++ -*-
// $Id: mre_initcg_extrap_predictor.h,v 1.1 2009-04-17 02:22:46 bjoo Exp $
/*! \file
 * \brief Minimal residual predictor - which takes also information
 *        from an EIG CG e-vec basis...
 *
 * Predictors for HMC
 */

#ifndef __mre_initcg_hmc_extrap_predictor_h__
#define __mre_initcg_hmc_extrap_predictor_h__

#include "chromabase.h"
#include "handle.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/circular_buffer.h"
#include "actions/ferm/invert/containers.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  
  /*! @ingroup predictor */
  namespace MREInitCG4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Minimal residual predictor
  /*! @ingroup predictor */
  class MREInitCG4DChronoPredictor  
    : public AbsChronologicalPredictor4D<LatticeFermion> 
  {
  private:
    Handle< CircularBuffer<LatticeFermion> > chrono_buf;

    void find_extrap_solution(LatticeFermion& psi, 
			      const LinearOperator<LatticeFermion>& A,
			      const LatticeFermion& chi);

    std::string opt_eigen_id;
    int Neig;   // No of vecs produced by EigCG per solve
    
  public:
    
    MREInitCG4DChronoPredictor(unsigned int max_chrono, const std::string& eigen_id, unsigned int max_evec) : chrono_buf(new CircularBuffer<LatticeFermion>(max_chrono)),  opt_eigen_id(eigen_id), Neig(max_evec) {}
    
    // Destructor is automagic
    ~MREInitCG4DChronoPredictor(void) {}

    // Do the hard work
    void operator()(LatticeFermion& psi, 
		    const LinearOperator<LatticeFermion>& A,
		    const LatticeFermion& chi);
    
    // No internal state so reset is a nop
    void reset(void) {
      chrono_buf->reset();
    }

    // Ignore new vector
    void newVector(const LatticeFermion& psi) 
    {
      START_CODE();

      QDPIO::cout << "MREPredictor: registering new solution. " << endl;
      if( chrono_buf->sizeMax() > 0 ) { 
        chrono_buf->push(psi);
      }

      
      QDPIO::cout << "MREPredictor: number of vectors stored is = " << chrono_buf->size() << endl;
        
      END_CODE();
    }

  };

  
  
} // End Namespace Chroma

#endif 
