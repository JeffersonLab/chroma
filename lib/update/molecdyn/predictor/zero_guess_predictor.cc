#include "chromabase.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma { 

  namespace ZeroGuess4DChronoPredictorEnv {

    // Create a new 4D Zero Guess Predictor
    // No params to read -- but preserve form
    AbsChronologicalPredictor4D<LatticeFermion>* createPredictor(XMLReader& xml,
								 const std::string& path) {

      // No params to read
      return new ZeroGuess4DChronoPredictor;
    }

     const std::string name = "ZERO_GUESS_4D_PREDICTOR";

    // Register it
    const bool registered = The4DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);
  
  };


  namespace ZeroGuess5DChronoPredictorEnv {
    // Create a new 5D Zero Guess Predictor
    // No params to read 
    AbsChronologicalPredictor5D<LatticeFermion>* createPredictor(const int N5,
								 XMLReader& xml,
								 const std::string& path) {
      return new ZeroGuess5DChronoPredictor(N5);
    };

    const std::string name = "ZERO_GUESS_5D_PREDICTOR";

    const bool registered = The5DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);

  };

};
