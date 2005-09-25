// -*- C++ -*-
// $Id: szin_write_obj_funcmap.h,v 2.0 2005-09-25 21:04:39 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#ifndef __szin_write_obj_funcmap_h__
#define __szin_write_obj_funcmap_h__

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
		TYPELIST_5(const string&,
			   const string&, 
			   int, int, int),
		void (*)(const string& buffer_id,
			 const string& filename, 
			 int j_decay, int t_start, int t_end),
		StringFunctionMapError> >
  TheSZINWriteObjFuncMap;


  //! Write object function map
  /*! \ingroup inlineio */
  namespace SZINWriteObjCallMapEnv
  { 
    extern bool registered;   // forward decl
  };

} // end namespace Chroma


#endif
