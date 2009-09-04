#include "chromabase.h"

namespace Chroma {
  namespace BiCGStabKernels {
    
    REAL64* _reduction_space;
    REAL64* _reduction_space_un;

    void  initScalarSiteKernels()
    {
      // Need Space for 5 innerprods, (10 doubles)
      // and            2 norms       ( 2 doubles)

#if 0
       QDPIO::cout << "Initing Reduction Space: There are " << qdpNumThreads() << " threads" << endl;
#endif
      _reduction_space_un = new REAL64 [ 12*qdpNumThreads()+16 ];
      _reduction_space = (REAL64*)((((ptrdiff_t)(_reduction_space_un)) + 15L)&(-16L));
    }

    void finishScalarSiteKernels()
    {
      delete [] _reduction_space_un;
    }


    REAL64* getNormSpace() { 
      return _reduction_space;
    }

 
  }
}
