// -*- C++ -*-
// $Id: writeobj_funcmap.h,v 1.1 2005-09-23 03:43:10 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#ifndef __writeobj_funcmap_h__
#define __writeobj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Write object function map
  /*! \ingroup inlineio */
  typedef SingletonHolder< 
    FunctionMap<void,
		std::string,
		TYPELIST_4(const string&,
			   const string&, 
			   QDP_volfmt_t, QDP_serialparallel_t),
		void (*)(const string& buffer_id,
			 const string& filename, 
			 QDP_volfmt_t volfmt, QDP_serialparallel_t serpar),
		StringFunctionMapError> >
  TheWriteObjFuncMap;


  //! Write object function map
  /*! \ingroup inlineio */
  namespace WriteObjCallMapEnv
  { 
    extern bool registered;   // forward decl
  };

} // end namespace Chroma


#endif
