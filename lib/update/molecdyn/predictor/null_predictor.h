// -*- C++ -*-
// $Id: null_predictor.h,v 3.1 2009-06-04 20:29:13 bjoo Exp $
/*! \file
 * \brief Null predictor: Leaves input x0 unchanged
 *
 * Predictors for HMC
 */

#ifndef __null_predictor_h__
#define __null_predictor_h__

#include "chromabase.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"


namespace Chroma 
{ 
  
  /*! @ingroup predictor */
  namespace Null4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Zero initial guess predictor
  /*! @ingroup predictor */
  class Null4DChronoPredictor : 
    public AbsTwoStepChronologicalPredictor4D<LatticeFermion> 
  {
  public:

    // Destructor is automagic
    ~Null4DChronoPredictor(void) {}

    // Zero out psi -- it is a zero guess after all
    void predictX(LatticeFermion& X,
		  const LinearOperator<LatticeFermion> &A,
		  const LatticeFermion& chi) 
    {
      START_CODE();

      QDPIO::cout << "Null Predictor Predict X: Leaving guess unchanged" << endl;

      END_CODE();
    }

    void predictY(LatticeFermion& Y,
		  const LinearOperator<LatticeFermion> &A,
		  const LatticeFermion& chi) 
    {
      START_CODE();
      QDPIO::cout << "Null Predictor Predict Y: Leaving guess unchanged" << endl;

      END_CODE();
    }
    


    // No internal state so reset is a nop
    void reset(void) {
    }

    // Ignore new vector
    void newXVector(const LatticeFermion& psi) {
      // Nothing
    }

    void newYVector(const LatticeFermion& psi) {
      // Nothing
    }

  };

  

  /*! @ingroup predictor */
  namespace Null5DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }
  
  //! Zero initial guess predictor
  /*! @ingroup predictor */
  class Null5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> 
  {
  public:
    ~Null5DChronoPredictor(void) {}

    // Creation
    Null5DChronoPredictor(const int N5_) {}

    // Copying
    Null5DChronoPredictor(const Null5DChronoPredictor& p)
     {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi,
		    const LinearOperatorArray<LatticeFermion>& A,
		    const multi1d<LatticeFermion>& chi) 
    { 
    }
    

    // No internal state so reset is a Nop
    void reset(void) {
    }

    // Ignore new vector
    void newVector(const multi1d<LatticeFermion>& psi) {
    }
    
  };
  
} // End Namespace Chroma

#endif 
