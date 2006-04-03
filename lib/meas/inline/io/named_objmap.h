// -*- C++ -*-
// $Id: named_objmap.h,v 3.0 2006-04-03 04:59:03 edwards Exp $
/*! \file
 *  \brief Named object function map
 */

#ifndef __named_objmap_h__
#define __named_objmap_h__

#include "singleton.h"
#include "named_obj.h"
#include "chromabase.h"

namespace Chroma
{

  // Turn into a Singleton. Create with CreateUsingNew
  // Has NoDestroy lifetime, as it may be needed for 
  // the destruction policy is No Destroy, so the 
  // Singleton is not cleaned up on exit. This is so 
  // that static objects can refer to it with confidence
  // in their own destruction, not having to worry that
  // atexit() may have destroyed the allocator before
  // the static objects need to feed memory. 
  /*! \ingroup inlineio */
  typedef SingletonHolder<NamedObjectMap,
			  QDP::CreateUsingNew,
			  QDP::NoDestroy,
			  QDP::SingleThreaded> TheNamedObjMap;
  
} // end namespace Chroma


#endif
