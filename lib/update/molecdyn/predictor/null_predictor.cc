#include "chromabase.h"
#include "update/molecdyn/predictor/null_predictor.h"


namespace Chroma 
{ 

  namespace Null4DChronoPredictorEnv 
  {
    namespace
    {
      // Create a new 4D Zero Guess Predictor
      // No params to read -- but preserve form
      AbsChronologicalPredictor4D<LatticeFermion>* createPredictor(XMLReader& xml,
								   const std::string& path) 
      {
	// No params to read
	return new Null4DChronoPredictor;
      }

      //! Local registration flag
      bool registered = false;
    }
    
    const std::string name = "NULL_4D_PREDICTOR";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= The4DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);
  	registered = true;
      }
      return success;
    }

  }


  namespace Null5DChronoPredictorEnv 
  {
    namespace
    {
      // Create a new 5D Zero Guess Predictor
      // No params to read 
      AbsChronologicalPredictor5D<LatticeFermion>* createPredictor(const int N5,
								   XMLReader& xml,
								   const std::string& path) 
      {
	return new Null5DChronoPredictor(N5);
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "NULL_5D_PREDICTOR";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= The5DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);
	registered = true;
      }
      return success;
    }
  }

}
