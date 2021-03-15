// -*- C++ -*-
/*! \file
 *  \brief Named object function std::map
 */

#ifndef __named_objmap_h__
#define __named_objmap_h__

#include "singleton.h"
#include "named_obj.h"
#include "chromabase.h"

namespace Chroma
{

  // Turn into a Singleton that it not going to be destroy, never!
  /*! \ingroup inlineio */
  typedef Chroma::SingletonHolder<NamedObjectMap> TheNamedObjMap;
  
} // end namespace Chroma


#endif
