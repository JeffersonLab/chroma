#ifndef __zero_guess_predictor_h__
#define __zero_guess_predictor_h__

#include "chromabase.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"

using namespace std;
using namespace QDP;
using namespace Chroma;

namespace Chroma { 
  
  namespace ZeroGuess4DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };

  class ZeroGuess4DChronoPredictor : 
    public AbsChronologicalPredictor4D<LatticeFermion> {
    
    public:

    // Destructor is automagic
    ~ZeroGuess4DChronoPredictor(void) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(LatticeFermion& psi) {
      QDPIO::cout << "ZeroGuessPredictor: zeroing initial guess" << endl;
      psi = zero;
    }
    
    // No internal state so reset is a nop
    void reset(void) {
      QDPIO::cout << "Resetting Chrono Predictor" << endl;
    }

    // Ignore new vector
    void newVector(const LatticeFermion& psi) {
      QDPIO::cout << "ZeroGuessPredictor: registering new solution (not)" << endl;
    }

  };

  

  namespace ZeroGuess5DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };
  
  class ZeroGuess5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> {

    public:
    ~ZeroGuess5DChronoPredictor(void) {}

    // Creation
    ZeroGuess5DChronoPredictor(const int N5_) : N5(N5_) {}

    // Copying
    ZeroGuess5DChronoPredictor(const ZeroGuess5DChronoPredictor& p) : 
      N5(p.N5) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi) { 

      QDPIO::cout << "ZeroGuessPredictor: zeroing initial guess" << endl;
      psi.resize(N5);
      psi = zero;


    }
    

    // No internal state so reset is a Nop
    void reset(void) {
      QDPIO::cout << "Resetting Chrono Predictor" << endl;
    }

    // Ignore new vector
    void newVector(const multi1d<LatticeFermion>& psi) {
      QDPIO::cout << "ZeroGuessPredictor: registering new solution (not)" << endl;
    }
    
    private:
    const int N5;
  };
  
}; // End Namespace Chroma

#endif 
