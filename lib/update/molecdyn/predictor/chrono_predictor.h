#ifndef chrono_predictor_h
#define chrono_predictor_h

#include "chromabase.h"
#include "linearop.h"

using namespace QDP;
using namespace Chroma;

namespace Chroma {
  
  // Abstract interface for a Chronological Solution predictor
  template<typename T>
  class AbsChronologicalPredictor4D {
  public:
    
    // Virtual destructor to help with cleanup
    virtual ~AbsChronologicalPredictor4D(void) {}

    // Set psi to be the next initial guess
    //
    // I have expanded the interface to allow us to 
    // pass the Matrix M, and the RHS chi as well as phi
    // 
    // We are trying to solve the system: 
    //            A  psi = chi
    //
    // for a CG situation A = MdagM
    //
    // and we are trying to get a guess for phi which 
    // minimises the initial residual.
    virtual void operator()(T& psi, 
			    const LinearOperator<T>& A, 
			    const T& chi) = 0;

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
    // I have expanded the interface to allow us to 
    // pass the Matrix M, and the RHS chi as well as phi
    // 
    // We are trying to solve the system
    //            A psi = chi
    // 
    // (for now in a CG based situation A = MdagM)
    //
    // and we are trying to get a guess for phi which 
    // minimises the initial residual.
    virtual void operator()(multi1d<T>& psi,
			    const LinearOperator<multi1d<T> >& A, 
			    const multi1d<T>& chi) = 0;
			    

    // Reset internal state (call this if the gauge field or 
    // pseudofermion fields change)
    virtual void reset(void) = 0;
    
    // Present new vector for use in future chronological
    // Predictors
    virtual void newVector(const multi1d<T>& psi) = 0;
  };

}; // End namespace
#endif
