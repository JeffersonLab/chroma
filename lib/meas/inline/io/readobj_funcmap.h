// -*- C++ -*-
// $Id: readobj_funcmap.h,v 1.1 2005-09-23 03:43:10 edwards Exp $
/*! \file
 *  \brief Read object function map
 */

#ifndef __readobj_funcmap_h__
#define __readobj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"

namespace Chroma
{

  //! Read object function map
  /*! \ingroup inlineio */
  typedef SingletonHolder< 
    FunctionMap<void, 
		std::string,
		TYPELIST_3(const string&,
			   const string&, 
			   QDP_serialparallel_t),
		void (*)(const string& buffer_id,
			 const string& filename, 
			 QDP_serialparallel_t serpar),
		StringFunctionMapError> >
  TheReadObjFuncMap;


  //! Read object function map
  /*! \ingroup inlineio */
  namespace ReadObjCallMapEnv
  { 
    extern bool registered;   // forward decl
  };

} // end namespace Chroma


#endif
