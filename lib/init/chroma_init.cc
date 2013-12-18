// $Id: chroma_init.cc,v 1.5 2005/01/24 21:31:32 edwards Exp 
/*! \file
 *  \brief Initialization of Chroma
 */

#include "chroma_config.h"

#if defined(BUILD_JIT_CLOVER_TERM)
#include "../actions/ferm/linop/clover_term_ptx_w.h"
#endif

#include "init/chroma_init.h"
#include "io/xmllog_io.h"

#ifdef BUILD_QUDA
#include <quda.h>
#endif

namespace Chroma 
{

  //! Private (anonymous) namespace
  namespace
  {
    //! The "DATA" filename
    string input_filename = "DATA";

    //! The "XMLDAT" filename
    string output_filename = "XMLDAT";

    //! The "XMLLOG" filename
    string log_filename = "XMLLOG";

    //! The "current working directory" -- prepended to above if set
    string cwd = ".";

    //! Has the xml output instance been created?
    bool xmlOutputP = false;

    //! Has the xml log instance been created?
    bool xmlLogP = false;

    //! Has the input instance been created?
    // bool xmlInputP = false;


    // Internal
    string constructFileName(const string& filename)
    {
      string ret_val;
      if ( filename[0] == '.' || filename[0]=='/'  ) {
	// Fully qualified pathname
	ret_val = filename;
      }
      else { 
	// Prepend CWD
	ret_val = getCWD() + "/" + filename;
      }
      return ret_val;
    }


  }; // End anonymous namespace

  //! Get input file name
  string getXMLInputFileName() {return constructFileName(input_filename);}

  //! Get output file name
  string getXMLOutputFileName() {return constructFileName(output_filename);}

  //! Get log file name
  string getXMLLogFileName() {return constructFileName(log_filename);}

  //! Get current working directory
  string getCWD() {return cwd;}


  //! Set input file name
  void setXMLInputFileName(const string& name) {input_filename = name;}

  //! Set output file name
  void setXMLOutputFileName(const string& name) {output_filename = name;}

  //! Set output logfile name
  void setXMLLogFileName(const string& name) {log_filename = name;}

  //! Set current working directory
  void setCWD(const string& name) {cwd = name;}


  //! Chroma initialisation routine
  void initialize(int* argc, char ***argv) 
  {
#ifndef QDP_IS_QDPJIT
    if (! QDP_isInitialized())
      QDP_initialize(argc, argv);
#else
    if (! QDP_isInitialized())
      QDP_initialize_CUDA(argc, argv);
#endif

    for(int i=0; i < *argc; i++) 
    {
      // Get argv[i] into a string
      string argv_i = string( (*argv)[i] );

      // Search for -i or --chroma-i
      if( argv_i == string("-h") || argv_i == string("--help") ) 
      {
	QDPIO::cerr << "Usage: " << (*argv)[0] << "  <options>" << endl
		    << "   -h           help\n"
		    << "   --help       help\n"
		    << "   -i           [" << getXMLInputFileName() << "]  xml input file name\n"
		    << "   --chroma-i   [" << getXMLInputFileName() << "]  xml input file name\n"
		    << "   -o           [" << getXMLOutputFileName() << "]  xml output file name\n"
		    << "   --chroma-p   [" << getXMLOutputFileName() << "]  xml output file name\n"
		    << "   -l           [" << getXMLLogFileName() << "]  xml output file name\n"
		    << "   --chroma-l   [" << getXMLLogFileName() << "]  xml log file name\n"
		    << "   -cwd         [" << getCWD() << "]  xml log file name\n"
		    << "   --chroma-cwd [" << getCWD() << "]  xml log file name\n"

		    
		    << endl;
	QDP_abort(0);
      }

      // Search for -i or --chroma-i
      if( argv_i == string("-i") || argv_i == string("--chroma-i") ) 
      {
	if( i + 1 < *argc ) 
	{
	  setXMLInputFileName(string( (*argv)[i+1] ));
	  // Skip over next
	  i++;
	}
	else 
	{
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -i specified. " << endl;
	  QDP_abort(1);
	}
      }
 
      // Search for -o or --chroma-o
      if( argv_i == string("-o") || argv_i == string("--chroma-o") ) 
      {
	if( i + 1 < *argc ) {
	  setXMLOutputFileName(string( (*argv)[i+1] ));
	  // Skip over next
	  i++;
	}
	else {
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -o specified. " << endl;
	  QDP_abort(1);
	}
      }
      
      // Search for -l or --chroma-l
      if( argv_i == string("-l") || argv_i == string("--chroma-l") ) 
      {
	if( i + 1 < *argc ) {
	  setXMLLogFileName(string( (*argv)[i+1] ));
	  // Skip over next
	  i++;
	}
	else {
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -l specified. " << endl;
	  QDP_abort(1);
	}
      }
      
      // Search for -cwd or --chroma-cwd
      if( argv_i == string("-cwd") || argv_i == string("--chroma-cwd") ) 
      {
	if( i + 1 < *argc ) {
	  setCWD(string( (*argv)[i+1] ));
	  // Skip over next
	  i++;
	}
	else {
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -cwd specified. " << endl;
	  QDP_abort(1);
	}
      }

    }


#ifdef QDP_IS_QDPJIT
#ifdef BUILD_QUDA
  std::cout << "Setting CUDA device" << endl;
  int cuda_device = QDP_setGPU();
  std::cout << "Setting QUDA verbosity to silent" << endl;
  setVerbosityQuda(QUDA_SILENT, "", stdout);
  //std::cout << "Setting QUDA verbosity to summarize" << endl;
  //setVerbosityQuda(QUDA_SUMMARIZE, "", stdout);
  std::cout << "Initializing QMP part" << endl;
  QDP_initialize_QMP(argc, argv);
  std::cout << "Initializing QUDA device (using CUDA device no. " << cuda_device << ")" << endl;
  initQudaDevice(cuda_device);
  std::cout << "Initializing QDP-JIT GPUs" << endl;
  QDP_startGPU();
  std::cout << "Initializing QUDA memory" << endl;
  initQudaMemory();
#else
  std::cout << "Setting device" << endl;
  QDP_setGPU();
  std::cout << "Initializing QMP part" << endl;
  QDP_initialize_QMP(argc, argv);
  QDPIO::cout << "Initializing start GPUs" << endl;
  QDP_startGPU();
#endif
#else
#ifdef BUILD_QUDA
  std::cout << "Initializing QUDA" << endl;
  initQuda(-1);
#endif
#endif




  }


  //! Chroma finalization routine
  void finalize(void)
  {

#ifdef BUILD_QUDA
    endQuda();
#endif

#if defined(BUILD_JIT_CLOVER_TERM)
    QDP_info_primary("Time for packForQUDA: %f sec",PackForQUDATimer::Instance().get() / 1.0e6);
#endif

    if (! QDP_isInitialized())
      return;

    
    /*
    if( xmlInputP ) { 
      Chroma::getXMLInputInstance().close();
    }
    */
    if( xmlOutputP ) { 
      Chroma::getXMLOutputInstance().close();
    }
    if( xmlLogP ) {
      Chroma::getXMLLogInstance().close();
    }

    QDP_finalize();
  }


  //! Chroma abort routine
  void abort(int i) 
  {
    QDP_abort(i);
  }


  //! Get xml output instance
  XMLFileWriter& getXMLOutputInstance()
  {
    if (! xmlOutputP)
    {
      try { 
	TheXMLOutputWriter::Instance().open(getXMLOutputFileName());
      }
      catch(...) { 
	QDPIO::cerr << "Unable to open " << getXMLOutputFileName() << endl;
	QDP_abort(1);
      }

      xmlOutputP = true;
    }

    return TheXMLOutputWriter::Instance();
  }
  
  //! Get xml log instance
  XMLFileWriter& getXMLLogInstance()
  {
    if (! xmlLogP)
    {
      try { 
	TheXMLLogWriter::Instance().open(getXMLLogFileName());
      }
      catch(...) {
	QDPIO::cerr << "Unable to open " << getXMLLogFileName() << endl;
	QDP_abort(1);
      }

      xmlLogP = true;
    }

    return TheXMLLogWriter::Instance();
  }
  
  /*
  //! Get xml input instance
  XMLReader& getXMLInputInstance()
  {
    if(! xmlInputP ) {
      try {
	TheXMLInputReader::Instance().open(getXMLInputFileName());
      }
      catch() { 
	QDPIO::cerr << "Unable to open " << pathname << endl;
	QDP_abort(1);
      }
      xmlInputP = true;
    }
    
    return TheXMLInputReader::Instance();
  }
  */


}
