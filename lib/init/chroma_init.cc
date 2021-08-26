/*! \file
 *  \brief Initialization of Chroma
 */

#include "chroma_config.h"

#include "singleton.h"
#include "init/chroma_init.h"
#include "io/xmllog_io.h"

#include "qdp_init.h"

#if defined(BUILD_JIT_CLOVER_TERM)
#if defined(QDPJIT_IS_QDPJITPTX)
#include "../actions/ferm/linop/clover_term_ptx_w.h"
#endif
#endif

#ifdef BUILD_QUDA
#include <quda.h>
#include <unistd.h>
#ifndef QDP_IS_QDPJIT
#include "init/local_rank.h"
#include <quda_api.h>
#include <device.h>
#endif
#endif

// Indlude it no-matter what
#ifdef BUILD_QPHIX
#include "../qphix_singleton.h"
#ifdef CHROMA_BUILDING_QPHIX_DSLASH
#include "qdp_datalayout.h"
#include "qphix/geometry.h"
#include "actions/ferm/invert/qphix/qphix_vec_traits.h"
#endif
#endif

#ifdef BUILD_MGPROTO
#include "MG_config.h"
#include "utils/memory.h"
#include "utils/initialize.h"
#ifdef MG_ENABLE_TIMERS
#include "utils/timer.h"
#endif
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


  } // End anonymous namespace

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
#if defined QDP_IS_QDPJIT
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
		    << "   -l           [" << getXMLLogFileName() << "]  xml log file name\n"
		    << "   --chroma-l   [" << getXMLLogFileName() << "]  xml log file name\n"
		    << "   -cwd         [" << getCWD() << "]  xml working directory\n"
		    << "   --chroma-cwd [" << getCWD() << "]  xml working directory\n"

		    
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

    // Good luck following the flow of the conditional compilation macros
#if defined QDP_IS_QDPJIT
#  ifdef BUILD_QUDA
    std::cout << "Setting CUDA device" << std::endl;
#    ifndef QDP_USE_COMM_SPLIT_INIT
    int cuda_device = QDP_setGPU();
#    endif
    std::cout << "Initializing QMP part" << std::endl;
    QDP_initialize_QMP(argc, argv);
#    ifdef QDP_USE_COMM_SPLIT_INIT
    int cuda_device = QDP_setGPUCommSplit();
#    endif
    setVerbosityQuda(QUDA_SUMMARIZE, "", stdout);

    QDPIO::cout << "Initializing QUDA device (using CUDA device no. " << cuda_device << ")"
		<< std::endl;

    initQudaDevice(cuda_device);
    QDPIO::cout << "Initializing QDP-JIT GPUs" << std::endl;
#    ifdef QDP_FIX_GPU_SETTING
    QDP_startGPU(cuda_device);
#    else
    QDP_startGPU();
#    endif
    QDPIO::cout << "Initializing QUDA memory" << std::endl;
    initQudaMemory();

#  else // BUILD_QUDA
    std::cout << "Setting device" << std::endl;
#    ifndef QDP_USE_COMM_SPLIT_INIT
    int dev = QDP_setGPU();
#    endif
    std::cout << "Initializing QMP part" << std::endl;
    QDP_initialize_QMP(argc, argv);
#    ifdef QDP_USE_COMM_SPLIT_INIT
    QDP_setGPUCommSplit();
#    endif
    QDPIO::cout << "Initializing start GPUs" << std::endl;
#    ifdef QDP_FIX_GPU_SETTING
    QDP_startGPU(dev);
#    else
    QDP_startGPU();
#    endif
#  endif // BUILD_CUDA

#else // defined QDP_IS_QDPJIT
#  ifdef BUILD_QUDA
   {
     char hostname[128];
     ::gethostname(hostname, sizeof(hostname)); 
     std::cout << "Initializing QUDA on local rank: " << localRank() << " on host: " << hostname <<  " with initQuda(-1)" <<  std::endl;
     initQuda(-1);
   }
#  endif
#endif

#ifdef BUILD_QPHIX
  QDPIO::cout << "Initializing QPhiX CLI Args" << std::endl;
  QPhiX::QPhiXCLIArgs& QPhiXArgs = TheQPhiXParams::Instance();
  QPhiXArgs.init((*argc),(*argv));
  QDPIO::cout << "QPhiX CLI Args Initialized" << std::endl;
  QDPIO::cout << " QPhiX: By="<< QPhiXArgs.getBy() << std::endl;
  QDPIO::cout << " QPhiX: Bz="<< QPhiXArgs.getBz() << std::endl;
  QDPIO::cout << " QPhiX: Pxy="<< QPhiXArgs.getPxy() << std::endl;
  QDPIO::cout << " QPhiX: Pxyz="<< QPhiXArgs.getPxyz() << std::endl;
  QDPIO::cout << " QPhiX: NCores="<< QPhiXArgs.getNCores() << std::endl;
  QDPIO::cout << " QPhiX: Sy="<< QPhiXArgs.getSy() << std::endl;
  QDPIO::cout << " QPhiX: Sz="<< QPhiXArgs.getSz() << std::endl;
  QDPIO::cout << " QPhiX: MinCt="<< QPhiXArgs.getMinCt() << std::endl;

  int num_threads = QPhiXArgs.getNCores()*QPhiXArgs.getSy()*QPhiXArgs.getSz();
  if( qdpNumThreads() != num_threads ) {
	  QDPIO::cerr << "ChromaInit: qdpNumThreads is different from NCores*Sy*Sz" << std::endl;
	  QDP_abort(1);
  }

 #if defined(QDP_IS_QDPJIT)
	// For QDP-JIT We insist on matchig QDP-JIT layout and QPhiX Layout for now

  	// Want ocsri layout. So pos_o = 0, pos_s=2, pos_c=1, pos_r=3, pos_i=4
  	if (! QDP_assert_jit_datalayout(0,2,1,3,4) ) {
  		QDPIO::cerr << "ChromaInit: DataLayout Ordering Mismatch. Wanted ocsri layout, but have ";
  		QDP_print_jit_datalayout();
  		QDP_abort(1);
  	}

  	// Check QDP Inner length is the SOALEN
  	int64_t layout_inner_size=QDP::getDataLayoutInnerSize();
  	if( layout_inner_size != CHROMA_QPHIX_SOALEN ) {
  		QDPIO::cout << "ChromaInit: Our SOA Length is " << CHROMA_QPHIX_SOALEN << " but QDP-JIT has inner="<< layout_inner_size << std::endl;
  		QDP_abort(1);
  	}
#endif
#endif

#ifdef BUILD_MGPROTO
  	// Initialzie MG Proto memory
  	QDPIO::cout << "Initializing MG_proto memory system" << std::endl;
  	MG::InitMemory(argc,argv);
#ifdef BUILD_QPHIX
  	QDPIO::cout << "Initializing QPhiX CLI Args for MG_proto" << std::endl;
  	MG::InitCLIArgs(argc,argv);
#endif
#endif

  }


  //! Chroma finalization routine
  void finalize(void)
  {

#ifdef BUILD_QUDA
    endQuda();
#endif

#ifdef BUILD_MGPROTO
#ifdef MG_ENABLE_TIMERS
    MG::Timer::TimerAPI::reportAllTimer();
#endif

    MG::FinalizeMemory();
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

    // Destroy singletons
    destroySingletons();

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
