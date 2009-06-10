#include "chromabase.h"

namespace Chroma {
  namespace BiCGStabKernels {
    
    REAL64* _reduction_space;

    void  initScalarSiteKernels()
    {
      _reduction_space = new REAL64 [ 3*qdpNumThreads() ];
    }

    void finishScalarSiteKernels()
    {
      delete [] _reduction_space;
    }


    REAL64* getNormSpace() { 
      return _reduction_space;
    }

 
  }
}
