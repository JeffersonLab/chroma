#include "chromabase.h"
#include "update/molecdyn/predictor/MG_predictor.h"


namespace Chroma 
{ 

  namespace MG4DChronoPredictorEnv 
  {
    namespace
    {
      AbsChronologicalPredictor4D<LatticeFermion>* createPredictor(XMLReader& xml,
								   const std::string& path) 
      {
	//Params read here.
	std::string subspace_name = "Default";
	int refresh_rate = 0;
	try 
	{
	  XMLReader paramtop(xml, path);
	  read( paramtop, "./SubspaceName", subspace_name);
	  read( paramtop, "./RefreshRate", refresh_rate);
	}
	catch( const std::string& e ) { 
	  QDPIO::cerr << "Caught exception reading XML: " << e << std::endl;
	  QDP_abort(1);
	}
	return new MG4DChronoPredictor<LatticeFermion>(subspace_name, refresh_rate);
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MG_4D_PREDICTOR";

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
}
