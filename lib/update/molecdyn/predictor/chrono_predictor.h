#ifndef chrono_predictor_h
#define chrono_predictor_h

#include "chromabase.h"


namespace Chroma {
  
  // Abstract interface for a Chronological Solution predictor
  template<typename T>
  class AbsChronologicalPredictor4D {
  public:
    
    // Virtual destructor to help with cleanup
    virtual ~AbsChronologicalPredictor4D(void) {}

    // Set psi to be the next initial guess
    virtual void operator()(T& psi) = 0;

    // Reset internal state (call this if the gauge field or 
    // pseudofermion fields change)
    virtual void reset(void) = 0;

    // Present new vector for use in future chronological
    // Predictors
    virtual void newVector(const T& psi) = 0;
  };
  

  template<typename T>
    class AbsChronologicalPredictor5D {
    public:

    virtual ~AbsChronologicalPredictor5D(void) {}

    // Set psi to be the next initial guess
    virtual void operator()(multi1d<T>& psi) = 0;

    // Reset internal state (call this if the gauge field or 
    // pseudofermion fields change)
    virtual void reset(void) = 0;
    
    // Present new vector for use in future chronological
    // Predictors
    virtual void newVector(const multi1d<T>& psi) = 0;
  };

}; // End namespace
#endif
