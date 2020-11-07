#include "chromabase.h"
#include "update/molecdyn/predictor/quda_predictor.h"


namespace Chroma 
{ 

  namespace QUDA4DChronoPredictorEnv 
  {
    namespace
    {
    // Create a new 4D Zero Guess Predictor
    // No params to read -- but preserve form
    AbsChronologicalPredictor4D<LatticeFermion>* createPredictor(XMLReader& xml,
    		const std::string& path)  {
    	unsigned int max_chrono = 1;
    	QudaPrecisionType prec = DEFAULT;

    	try {

    		XMLReader paramtop(xml, path);
    		read( paramtop, "./MaxChrono", max_chrono);
    		if ( paramtop.count("./Precision") == 1) {
    			read(paramtop, "./Precision", prec);
    		}
    	}
    	catch( const std::string& e ) {
    		QDPIO::cerr << "Caught exception reading XML: " << e << std::endl;
    		QDP_abort(1);
    	}

    	// No params to read
    	return new QUDA4DChronoPredictor(max_chrono,prec);
    }

      //! Local registration flag
      bool registered = false;

      //! Local Chrono index
      int chrono_index = 0;
      
    }
    
    const std::string name = "QUDA_4D_PREDICTOR";

    int getAndIncrGlobalQUDAChronoIndex()
    {
      int ret_val = chrono_index;
      ++chrono_index;
      return ret_val;
    }
    
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
