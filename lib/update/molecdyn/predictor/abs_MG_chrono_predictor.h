// abs_MG_chrono_predictor.h
//Abstract class for MG chrono predictor - Arjun. 
/*
 *
 * Predictors for HMC
 */

#ifndef __abs_MG_chrono_predictor_h__
#define __abs_MG_chrono_predictor_h__

#include "chromabase.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"


namespace Chroma 
{ 
  
  template<typename T>
  class AbsMGChronologicalPredictor4D : public AbsChronologicalPredictor4D<T> {
  public:
    
    // Virtual destructor to help with cleanup
    virtual ~AbsMGChronologicalPredictor4D(void) {}

    
    virtual void getSubspace() = 0;


    virtual void resetSubspace(int counter) = 0;

    virtual void reset(void) = 0;


    virtual void operator()(T& psi, 
			    const LinearOperator<T>& A, 
			    const T& chi) = 0;

    virtual void newVector(const T& psi) = 0;

  };
  
} // End Namespace Chroma

#endif 
