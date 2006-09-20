#include "chromabase.h"
#include "update/molecdyn/predictor/linear_extrap_predictor.h"



namespace Chroma 
{ 

  namespace LinearExtrapolation4DChronoPredictorEnv 
  {
    namespace
    {
      // Create a new 4D Zero Guess Predictor
      // No params to read -- but preserve form
      AbsChronologicalPredictor4D<LatticeFermion>* createPredictor(XMLReader& xml,
								   const std::string& path) 
      {
	// No params to read
	return new LinearExtrapolation4DChronoPredictor;
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LINEAR_EXTRAPOLATION_4D_PREDICTOR";

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


  namespace LinearExtrapolation5DChronoPredictorEnv 
  {
    namespace
    {
      // Create a new 5D Zero Guess Predictor
      // No params to read 
      AbsChronologicalPredictor5D<LatticeFermion>* createPredictor(const int N5,
								   XMLReader& xml,
								   const std::string& path) 
      {
	return new LinearExtrapolation5DChronoPredictor(N5);
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LINEAR_EXTRAPOLATION_5D_PREDICTOR";

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
