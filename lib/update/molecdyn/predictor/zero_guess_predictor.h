// -*- C++ -*-
/*! \file
 * \brief Zero initial guess predictor
 *
 * Predictors for HMC
 */

#ifndef __zero_guess_predictor_h__
#define __zero_guess_predictor_h__

#include "chromabase.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"


namespace Chroma 
{ 
  
  /*! @ingroup predictor */
  namespace ZeroGuess4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Zero initial guess predictor
  /*! @ingroup predictor */
  class ZeroGuess4DChronoPredictor : 
    public AbsTwoStepChronologicalPredictor4D<LatticeFermion> 
  {
  public:

    // Destructor is automagic
    ~ZeroGuess4DChronoPredictor(void) {}

    // Zero out psi -- it is a zero guess after all
    void predictX(LatticeFermion& X,
		  const LinearOperator<LatticeFermion> &A,
		  const LatticeFermion& chi) 
    {
      START_CODE();

      QDPIO::cout << "ZeroGuessPredictor: zeroing initial guess" << std::endl;
      X = zero;

      END_CODE();
    }

    void predictY(LatticeFermion& Y,
		  const LinearOperator<LatticeFermion> &A,
		  const LatticeFermion& chi) 
    {
      START_CODE();

      QDPIO::cout << "ZeroGuessPredictor: zeroing initial guess" << std::endl;
      Y = zero;

      END_CODE();
    }
    


    // No internal state so reset is a nop
    void reset(void) {
    }

    // Ignore new std::vector
    void newXVector(const LatticeFermion& psi) {
      QDPIO::cout << "ZeroGuessPredictor: registering new solution (not)" << std::endl;
    }

    void newYVector(const LatticeFermion& psi) {
      QDPIO::cout << "ZeroGuessPredictor: registering new solution (not)" << std::endl;
    }

  };

  

  /*! @ingroup predictor */
  namespace ZeroGuess5DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }
  
  //! Zero initial guess predictor
  /*! @ingroup predictor */
  class ZeroGuess5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> 
  {
  public:
    ~ZeroGuess5DChronoPredictor(void) {}

    // Creation
    ZeroGuess5DChronoPredictor(const int N5_) : N5(N5_) {}

    // Copying
    ZeroGuess5DChronoPredictor(const ZeroGuess5DChronoPredictor& p) : 
      N5(p.N5) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi,
		    const LinearOperatorArray<LatticeFermion>& A,
		    const multi1d<LatticeFermion>& chi) 
    { 
      START_CODE();

      if (A.size() != N5)
      {
	QDPIO::cerr << "ZeroGuess5D: mismatched sizes A.size=" << A.size() 
		    << "  and N5=" << N5 << std::endl;
	QDP_abort(1);
      }
      psi.resize(N5);
      psi = zero;
    
      END_CODE();
    }
    

    // No internal state so reset is a Nop
    void reset(void) {
    }

    // Ignore new std::vector
    void newVector(const multi1d<LatticeFermion>& psi) {
      QDPIO::cout << "ZeroGuessPredictor: registering new solution (not)" << std::endl;
    }
    
  private:
    const int N5;
  };
  
} // End Namespace Chroma

#endif 
