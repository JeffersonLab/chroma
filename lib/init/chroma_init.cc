// $Id: chroma_init.cc,v 1.5 2005/01/24 21:31:32 edwards Exp 
/*! \file
 *  \brief Initialization of Chroma
 */

#include "chroma_config.h"

#include "init/chroma_init.h"
#include "io/xmllog_io.h"

#if defined(BUILD_JIT_CLOVER_TERM)
#if defined(QDPJIT_IS_QDPJITPTX)
#include "../actions/ferm/linop/clover_term_ptx_w.h"
#endif
#endif

#ifdef BUILD_QUDA
#include <quda.h>
#endif

namespace Chroma 
{

  //! Private (anonymous) namespace
  namespace
  {
    //! The "DATA" filename
    std::string input_filename = "DATA";

    //! The "XMLDAT" filename
    std::string output_filename = "XMLDAT";

    //! The "XMLLOG" filename
    std::string log_filename = "XMLLOG";

    //! The "current working directory" -- prepended to above if set
    std::string cwd = ".";

    //! Has the xml output instance been created?
    bool xmlOutputP = false;

    //! Has the xml log instance been created?
    bool xmlLogP = false;

    //! Has the input instance been created?
    // bool xmlInputP = false;


    // Internal
    std::string constructFileName(const std::string& filename)
    {
      std::string ret_val;
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
  std::string getXMLInputFileName() {return constructFileName(input_filename);}

  //! Get output file name
  std::string getXMLOutputFileName() {return constructFileName(output_filename);}

  //! Get log file name
  std::string getXMLLogFileName() {return constructFileName(log_filename);}

  //! Get current working directory
  std::string getCWD() {return cwd;}


  //! Set input file name
  void setXMLInputFileName(const std::string& name) {input_filename = name;}

  //! Set output file name
  void setXMLOutputFileName(const std::string& name) {output_filename = name;}

  //! Set output logfile name
  void setXMLLogFileName(const std::string& name) {log_filename = name;}

  //! Set current working directory
  void setCWD(const std::string& name) {cwd = name;}


  //! Chroma initialisation routine
  void initialize(int* argc, char ***argv) 
  {
#if defined QDPJIT_IS_QDPJITPTX || defined QDPJIT_IS_QDPJITNVVM
    if (! QDP_isInitialized())
      QDP_initialize_CUDA(argc, argv);
#else
    if (! QDP_isInitialized())
      QDP_initialize(argc, argv);
#endif

    for(int i=0; i < *argc; i++) 
    {
      // Get argv[i] into a std::string
      std::string argv_i = std::string( (*argv)[i] );

      // Search for -i or --chroma-i
      if( argv_i == std::string("-h") || argv_i == std::string("--help") ) 
      {
	QDPIO::cerr << "Usage: " << (*argv)[0] << "  <options>" << std::endl
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

		    
		    << std::endl;
	QDP_abort(0);
      }

      // Search for -i or --chroma-i
      if( argv_i == std::string("-i") || argv_i == std::string("--chroma-i") ) 
      {
	if( i + 1 < *argc ) 
	{
	  setXMLInputFileName(std::string( (*argv)[i+1] ));
	  // Skip over next
	  i++;
	}
	else 
	{
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -i specified. " << std::endl;
	  QDP_abort(1);
	}
      }
 
      // Search for -o or --chroma-o
      if( argv_i == std::string("-o") || argv_i == std::string("--chroma-o") ) 
      {
	if( i + 1 < *argc ) {
	  setXMLOutputFileName(std::string( (*argv)[i+1] ));
	  // Skip over next
	  i++;
	}
	else {
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -o specified. " << std::endl;
	  QDP_abort(1);
	}
      }
      
      // Search for -l or --chroma-l
      if( argv_i == std::string("-l") || argv_i == std::string("--chroma-l") ) 
      {
	if( i + 1 < *argc ) {
	  setXMLLogFileName(std::string( (*argv)[i+1] ));
	  // Skip over next
	  i++;
	}
	else {
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -l specified. " << std::endl;
	  QDP_abort(1);
	}
      }
      
      // Search for -cwd or --chroma-cwd
      if( argv_i == std::string("-cwd") || argv_i == std::string("--chroma-cwd") ) 
      {
	if( i + 1 < *argc ) {
	  setCWD(std::string( (*argv)[i+1] ));
	  // Skip over next
	  i++;
	}
	else {
	  // i + 1 is too big
	  QDPIO::cerr << "Error: dangling -cwd specified. " << std::endl;
	  QDP_abort(1);
	}
      }

    }


#if defined QDPJIT_IS_QDPJITPTX || defined QDPJIT_IS_QDPJITNVVM
#ifdef BUILD_QUDA
  std::cout << "Setting CUDA device" << std::endl;
  int cuda_device = QDP_setGPU();
  //std::cout << "Setting QUDA verbosity to silent" << std::endl;
  //setVerbosityQuda(QUDA_SILENT, "", stdout);
  //std::cout << "Setting QUDA verbosity to summarize" << std::endl;
  //setVerbosityQuda(QUDA_SUMMARIZE, "", stdout);
  std::cout << "Initializing QMP part" << std::endl;
  QDP_initialize_QMP(argc, argv);
  QDPIO::cout << "Initializing QUDA device (using CUDA device no. " << cuda_device << ")" << std::endl;
  initQudaDevice(cuda_device);
  QDPIO::cout << "Initializing QDP-JIT GPUs" << std::endl;
  QDP_startGPU();
  QDPIO::cout << "Initializing QUDA memory" << std::endl;
  initQudaMemory();
#else
  std::cout << "Setting device" << std::endl;
  QDP_setGPU();
  std::cout << "Initializing QMP part" << std::endl;
  QDP_initialize_QMP(argc, argv);
  QDPIO::cout << "Initializing start GPUs" << std::endl;
  QDP_startGPU();
#endif
#else
#ifdef BUILD_QUDA
  std::cout << "Initializing QUDA" << std::endl;
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
#if defined(QDPJIT_IS_QDPJITPTX)
    QDP_info_primary("Time for packForQUDA: %f sec",PackForQUDATimer::Instance().get() / 1.0e6);
#endif
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
	QDPIO::cerr << "Unable to open " << getXMLOutputFileName() << std::endl;
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
	QDPIO::cerr << "Unable to open " << getXMLLogFileName() << std::endl;
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
	QDPIO::cerr << "Unable to open " << pathname << std::endl;
	QDP_abort(1);
      }
      xmlInputP = true;
    }
    
    return TheXMLInputReader::Instance();
  }
  */


}
