// -*- C++ -*-
// $Id: mre_extrap_predictor.h,v 3.5 2009-06-01 20:39:37 bjoo Exp $
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
  namespace MinimalResidualExtrapolation4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Minimal residual predictor
  /*! @ingroup predictor */
  class MinimalResidualExtrapolation4DChronoPredictor  
    : public AbsTwoStepChronologicalPredictor4D<LatticeFermion> 
  {
  private:
    Handle< CircularBuffer<LatticeFermion> > chrono_bufX;
    Handle< CircularBuffer<LatticeFermion> > chrono_bufY;
    
    void find_extrap_solution(LatticeFermion& psi, 
			      const LinearOperator<LatticeFermion>& A,
			      const LatticeFermion& chi,
			      const Handle<CircularBuffer<LatticeFermion> >& chrono_buf,
			      enum PlusMinus isign);
  public:
    
    MinimalResidualExtrapolation4DChronoPredictor(unsigned int max_chrono) : 
      chrono_bufX(new CircularBuffer<LatticeFermion>(max_chrono)),
      chrono_bufY(new CircularBuffer<LatticeFermion>(max_chrono)) {}
    
    // Destructor is automagic
    ~MinimalResidualExtrapolation4DChronoPredictor(void) {}

    // Do the hard work
    void predictX(LatticeFermion& X, 
		    const LinearOperator<LatticeFermion>& A,
		    const LatticeFermion& chi);

    void predictY(LatticeFermion& Y, 
		    const LinearOperator<LatticeFermion>& A,
		    const LatticeFermion& chi);
    
    // No internal state so reset is a nop
    void reset(void) {
      chrono_bufX->reset();
      chrono_bufY->reset();
    }


    void newXVector(const LatticeFermion& X) 
    {
      START_CODE();

      QDPIO::cout << "MREPredictor: registering new X solution. " << endl;
      chrono_bufX->push(X);
      QDPIO::cout << "MREPredictor: number of X vectors stored is = " << chrono_bufX->size() << endl;
    
      END_CODE();
    }


    void newYVector(const LatticeFermion& Y) 
    {
      START_CODE();

      QDPIO::cout << "MREPredictor: registering new Y solution. " << endl;
      chrono_bufY->push(Y);
      QDPIO::cout << "MREPredictor: number of Y vectors stored is = " << chrono_bufY->size() << endl;
    
      END_CODE();
    }

  };

  

  /*! @ingroup predictor */
  namespace MinimalResidualExtrapolation5DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  
  //! Minimal residual predictor
  /*! @ingroup predictor */
  class MinimalResidualExtrapolation5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> {
    
  private: 
    Handle< CircularBufferArray<LatticeFermion>  > chrono_buf;
    const int N5;

    void find_extrap_solution(multi1d<LatticeFermion>& psi, 
			      const LinearOperatorArray<LatticeFermion>& A,
			      const multi1d<LatticeFermion>& chi);

  public:
    MinimalResidualExtrapolation5DChronoPredictor(const int N5_, const unsigned int max_chrono) : chrono_buf( new CircularBufferArray<LatticeFermion>(max_chrono, N5_) ), N5(N5_) {}

      
    ~MinimalResidualExtrapolation5DChronoPredictor(void) {}
    
    // Copying
    MinimalResidualExtrapolation5DChronoPredictor(const MinimalResidualExtrapolation5DChronoPredictor& p) : chrono_buf(p.chrono_buf), N5(p.N5) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi,
		    const LinearOperatorArray<LatticeFermion>& A,
		    const multi1d<LatticeFermion>& chi); 

    // No internal state so reset is a Nop
    void reset(void) {
      chrono_buf->reset();

    }

    // Ignore new vector
    // Ignore new vector
    void newVector(const multi1d<LatticeFermion>& psi) 
    {
      START_CODE();

      QDPIO::cout << "MRE Predictor: registering new solution. " << endl;
      chrono_buf->push(psi);
      QDPIO::cout << "MRE Predictor: number of vectors stored is = " << chrono_buf->size() << endl;
    
      END_CODE();
    }

  };
  
} // End Namespace Chroma

#endif 
