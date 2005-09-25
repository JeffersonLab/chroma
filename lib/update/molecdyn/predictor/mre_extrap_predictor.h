// -*- C++ -*-
// $Id: mre_extrap_predictor.h,v 2.0 2005-09-25 21:04:43 edwards Exp $
/*! \file
 * \brief Minimal residual predictor
 *
 * Predictors for HMC
 */

#ifndef __mre_extrap_predictor_h__
#define __mre_extrap_predictor_h__

#include "chromabase.h"
#include "handle.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/circular_buffer.h"

namespace Chroma 
{ 
  
  /*! @ingroup predictor */
  namespace MinimalResidualExtrapolation4DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };

  //! Minimal residual predictor
  /*! @ingroup predictor */
  class MinimalResidualExtrapolation4DChronoPredictor  
    : public AbsChronologicalPredictor4D<LatticeFermion> {

  private:
    Handle< CircularBuffer<LatticeFermion> > chrono_buf;
    
    void find_extrap_solution(LatticeFermion& psi, 
			      const LinearOperator<LatticeFermion>& A,
			      const LatticeFermion& chi);
  public:
    
    MinimalResidualExtrapolation4DChronoPredictor(unsigned int max_chrono) : chrono_buf(new CircularBuffer<LatticeFermion>(max_chrono)) {}
    
    // Destructor is automagic
    ~MinimalResidualExtrapolation4DChronoPredictor(void) {}

    // Do the hard work
    void operator()(LatticeFermion& psi, 
		    const LinearOperator<LatticeFermion>& A,
		    const LatticeFermion& chi);
    
    // No internal state so reset is a nop
    void reset(void) {
      QDPIO::cout << "Resetting Chrono Predictor" << endl;
      chrono_buf->reset();
    }

    // Ignore new vector
    void newVector(const LatticeFermion& psi) {
      QDPIO::cout << "MREPredictor: registering new solution. " << endl;
      chrono_buf->push(psi);
      QDPIO::cout << "MREPredictor: number of vectors stored is = " << chrono_buf->size() << endl;
    }

  };

  

  /*! @ingroup predictor */
  namespace MinimalResidualExtrapolation5DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };
  
  //! Minimal residual predictor
  /*! @ingroup predictor */
  class MinimalResidualExtrapolation5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> {
    
  private: 
    Handle< CircularBufferArray<LatticeFermion>  > chrono_buf;
    const int N5;

    void find_extrap_solution(multi1d<LatticeFermion>& psi, 
			      const LinearOperator<multi1d<LatticeFermion> >& A,
			      const multi1d<LatticeFermion>& chi);

  public:
    MinimalResidualExtrapolation5DChronoPredictor(const int N5_, const unsigned int max_chrono) : chrono_buf( new CircularBufferArray<LatticeFermion>(max_chrono, N5_) ), N5(N5_) {}

      
    ~MinimalResidualExtrapolation5DChronoPredictor(void) {}
    
    // Copying
    MinimalResidualExtrapolation5DChronoPredictor(const MinimalResidualExtrapolation5DChronoPredictor& p) : chrono_buf(p.chrono_buf), N5(p.N5) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi,
		    const LinearOperator<multi1d<LatticeFermion> >& A,
		    const multi1d<LatticeFermion>& chi); 

    // No internal state so reset is a Nop
    void reset(void) {
      QDPIO::cout << "Resetting Chrono Predictor" << endl;
      chrono_buf->reset();

    }

    // Ignore new vector
    // Ignore new vector
    void newVector(const multi1d<LatticeFermion>& psi) {
      QDPIO::cout << "MRE Predictor: registering new solution. " << endl;
      chrono_buf->push(psi);
      QDPIO::cout << "MRE Predictor: number of vectors stored is = " << chrono_buf->size() << endl;
    }


  };
  
}; // End Namespace Chroma

#endif 
