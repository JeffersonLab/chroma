// -*- C++ -*-
// $Id: szin_read_obj_funcmap.h,v 2.1 2005-11-01 22:00:01 edwards Exp $
/*! \file
 *  \brief Read object function map
 */

#ifndef __szin_read_obj_funcmap_h__
#define __szin_read_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Read object function map
  /*! \ingroup inlineio */
  typedef SingletonHolder< 
    FunctionMap<void,
		std::string,
		TYPELIST_2(const string&,
			   const string&),
		void (*)(const string& buffer_id,
			 const string& filename),
		StringFunctionMapError> >
  TheSZINReadObjFuncMap;


  //! Read object function map
  /*! \ingroup inlineio */
  namespace SZINReadObjCallMapEnv
  { 
    extern bool registered;   // forward decl
  };

} // end namespace Chroma


#endif
