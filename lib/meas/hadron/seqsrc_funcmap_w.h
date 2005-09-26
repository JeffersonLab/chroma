// -*- C++ -*-
// $Id: seqsrc_funcmap_w.h,v 2.1 2005-09-26 04:48:35 edwards Exp $
/*! \file
 *  \brief Sequential source function map
 */

#ifndef __seqsrc_funcmap_h__
#define __seqsrc_funcmap_h__

#include "singleton.h"
#include "funcmap.h"

namespace Chroma
{

  //! Sequential Source function map 
  /*! \ingroup hadron */
  typedef SingletonHolder< 
    FunctionMap<LatticePropagator, 
		std::string,
		TYPELIST_1(const multi1d<LatticePropagator>&),
		LatticePropagator (*)(const multi1d<LatticePropagator>&),
		StringFunctionMapError> >
  TheSeqSourceFuncMap;


  namespace SeqSourceCallMapEnv
  { 
    extern bool registered;   // forward decl
  };

} // end namespace Chroma


#endif
