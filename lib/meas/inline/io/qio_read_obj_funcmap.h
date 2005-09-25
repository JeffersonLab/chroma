// -*- C++ -*-
// $Id: qio_read_obj_funcmap.h,v 2.0 2005-09-25 21:04:39 edwards Exp $
/*! \file
 *  \brief Read object function map
 */

#ifndef __qio_read_obj_funcmap_h__
#define __qio_read_obj_funcmap_h__

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
  TheQIOReadObjFuncMap;


  //! Read object function map
  /*! \ingroup inlineio */
  namespace QIOReadObjCallMapEnv
  { 
    extern bool registered;   // forward decl
  };

} // end namespace Chroma


#endif
