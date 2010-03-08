// -*- C++ -*-
#ifndef __PQ_TRAITS_H__
#define __PQ_TRAITS_H__

#include "chromabase.h"

using namespace QDP;

namespace Chroma {

  template<typename T>
    struct PQTraits {
      typedef  T  Base_t;
    };

  template<typename T>
    struct PQTraits< multi1d<T> > {
      typedef T  Base_t;
    };

};



#endif
