#include "chromabase.h"
#include "update/molecdyn/last_solution_predictor.h"

using namespace std;
using namespace QDP;
using namespace Chroma;

namespace Chroma { 

  namespace LastSolution4DChronoPredictorEnv {

    // Create a new 4D Zero Guess Predictor
    // No params to read -- but preserve form
    AbsChronologicalPredictor4D<LatticeFermion>* createPredictor(XMLReader& xml,
								 const std::string& path) {

      // No params to read
      return new LastSolution4DChronoPredictor;
    }

     const std::string name = "LAST_SOLUTION_4D_PREDICTOR";

    // Register it
    const bool registered = The4DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);
  
  };


  namespace LastSolution5DChronoPredictorEnv {
    // Create a new 5D Zero Guess Predictor
    // No params to read 
    AbsChronologicalPredictor5D<LatticeFermion>* createPredictor(const int N5,
								 XMLReader& xml,
								 const std::string& path) {
      return new LastSolution5DChronoPredictor(N5);
    };

    const std::string name = "LAST_SOLUTION_5D_PREDICTOR";

    const bool registered = The5DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);

  };

};
