// $Id: chroma_init.cc,v 1.4 2005-01-17 21:08:29 edwards Exp $


//#include "chroma.h"
#include "init/chroma_init.h"
#include "io/xmllog_io.h"
#include <string>

namespace Chroma { 

  // WARNING THIS CURRENTLY CAUSES A SEGFAULT IN PROPAGATOR SO I AM 
  // NOT USING IT. LINKAGE HACK STILL NEEDS TO STAY.
  // bool ChromaLinkage(void);


  void ChromaInitialize(int* argc, char ***argv) 
  {
    QDP_initialize(argc, argv);
    

    bool output_found = false;
    std::string output_filename;
    
    
    for( int i=0; i < *argc; i++) {

      // Get argv[i] into a string
      std::string argv_i = std::string( (*argv)[i] );

      /*
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
      */
 
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
    
    /*
    if ( ! input_found ) { 
      input_filename = "./DATA";
    }
    */

    if( ! output_found ) { 
      output_filename = "./XMLLOG";
    }
    

    /*
    QDPIO::cout << "ChromaInitialize: Reading XML Parameters from : " 
		<< input_filename << endl;
		TheXMLInputReader::Instance().open(input_filename);
    */
    QDPIO::cout << "ChromaInitialize: Writing XML Log file to     : " 
		<< output_filename << endl;
    
    TheXMLOutputWriter::Instance().open(output_filename);


    
    // Linkage
    // WARNING THIS CURRENTLY CAUSES A SEGFAULT IN PROPAGATOR SO I AM 
    // NOT USING IT. LINKAGE HACK STILL NEEDS TO STAY.
    //    bool foo = ChromaLinkage();

    // Init done
  }

  void ChromaFinalize(void) {
    // Finish
    // Chroma Finalize Stuff
    /* TheXMLInputReader::Instance().close(); */
    TheXMLOutputWriter::Instance().close();
    QDP_finalize();
  }

  void ChromaAbort(int i) { 
    // Cleanup
    /*    TheXMLInputReader::Instance().close(); */
    TheXMLOutputWriter::Instance().flush();
    TheXMLOutputWriter::Instance().close();
    QDP_abort(i);
  }


  //! To insure linking of code, place the registered code flags here
  /*! This is the bit of code that dictates what fermacts are in use */

  // WARNING THIS CURRENTLY CAUSES A SEGFAULT IN PROPAGATOR SO I AM 
  // NOT USING IT. LINKAGE HACK STILL NEEDS TO STAY.
  /*
  bool ChromaLinkage(void)
  {
    bool foo = true;
    
    // Propagator Fermact Actions 
    // 4D actions
    foo &= EvenOddPrecWilsonFermActEnv::registered;
    foo &= UnprecWilsonFermActEnv::registered;
    foo &= OvlapPartFrac4DFermActEnv::registered;
    foo &= EvenOddPrecParWilsonFermActEnv::registered;
    foo &= UnprecParWilsonFermActEnv::registered;
    foo &= UnprecDWFTransfFermActEnv::registered;

    // 5D actions
    foo &= EvenOddPrecDWFermActArrayEnv::registered;
    foo &= UnprecDWFermActArrayEnv::registered;
    foo &= EvenOddPrecNEFFermActArrayEnv::registered;
    foo &= UnprecNEFFermActArrayEnv::registered;
    foo &= UnprecOvlapContFrac5DFermActArrayEnv::registered;
    foo &= UnprecHTContFrac5DFermActArrayEnv::registered;
    foo &= EvenOddPrecOvlapContFrac5DFermActArrayEnv::registered;
    foo &= UnprecOvDWFermActArrayEnv::registered;
    foo &= EvenOddPrecOvDWFermActArrayEnv::registered;
    foo &= UnprecOvExtFermActArrayEnv::registered;
    foo &= UnprecZoloNEFFermActArrayEnv::registered;
    foo &= EvenOddPrecZoloNEFFermActArrayEnv::registered;
    foo &= EvenOddPrecKNOFermActArrayEnv::registered;
    foo &= UnprecDWFTransfFermActEnv::registered;
    
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
  */

};
