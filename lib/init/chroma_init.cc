#include "init/chroma_init.h"

using namespace QDP;
using namespace Chroma;
using namespace std;


namespace Chroma { 


  bool ChromaLinkage(void);


  void ChromaInitialize(int* argc, char ***argv) 
  {
    QDP_initialize(argc, argv);
    
    // Process Chroma type arguments
    bool input_found = false;
    bool output_found = false;
    
    std::string input_filename;
    std::string output_filename;
    
    
    for( int i=0; i < *argc; i++) {
      
      // Search for -i or --input-xml
      std::string argv_i = std::string( (*argv)[i] );
      
      // Search for -i or --input-xml
      if( argv_i == std::string("-i") ) {
	if( i + 1 < *argc ) {
	  input_filename = std::string( (*argv)[i+1] );
	  input_found = true;
	  // Skip over next
	  i++;
	}
	else {
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -i specified. " << endl;
	  QDP_abort(1);
	}
      }
      
      if( argv_i == std::string("-o") ) {
	if( i + 1 < *argc ) {
	  output_filename = std::string( (*argv)[i+1] );
	  output_found = true;
	  // Skip over next
	  i++;
	}
	else {
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -o specified. " << endl;
	  QDP_abort(1);
	}
      }
      
    }
    
    if ( ! input_found ) { 
      input_filename = "./DATA";
    }
    
    if( ! output_found ) { 
      output_filename = "./XMLDAT";
    }
    
    QDPIO::cout << "ChromaInitialize: Reading XML Parameters from : " 
		<< input_filename << endl;
    QDPIO::cout << "ChromaInitialize: Writing XML Log file to     : " 
		<< output_filename << endl;
    
    TheXMLInputReader::Instance().open(input_filename);
    TheXMLOutputWriter::Instance().open(output_filename);
    
    // Linkage
    bool foo = ChromaLinkage();

    // Init done
  }

  void ChromaFinalize(void) {
    // Finish
    // Chroma Finalize Stuff
    TheXMLInputReader::Instance().close();
    TheXMLOutputWriter::Instance().close();
    QDP_finalize();
  }

  void ChromaAbort(int i) { 
    // Cleanup
    TheXMLInputReader::Instance().close();
    TheXMLOutputWriter::Instance().flush();
    TheXMLOutputWriter::Instance().close();
    QDP_abort(i);
  }


  //! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
  bool ChromaLinkage(void)
  {
    bool foo = true;
    
    // Gauge Monomials
    foo &= GaugeMonomialEnv::registered;
    
    // 4D Ferm Monomials
    foo &= UnprecTwoFlavorWilsonTypeFermMonomialEnv::registered;
    foo &= EvenOddPrecTwoFlavorWilsonTypeFermMonomialEnv::registered;
    
    // 5D Ferm Monomials
    foo &= UnprecTwoFlavorWilsonTypeFermMonomial5DEnv::registered;
    foo &= EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DEnv::registered;
    
    // MD Integrators
    foo &= LatColMatPQPLeapfrogIntegratorEnv::registered;
    
    // Chrono predictor
    foo &= ZeroGuess4DChronoPredictorEnv::registered;
    foo &= ZeroGuess5DChronoPredictorEnv::registered;
    foo &= LastSolution4DChronoPredictorEnv::registered;  
    foo &= LastSolution5DChronoPredictorEnv::registered;
    
    return foo;
  }


};
