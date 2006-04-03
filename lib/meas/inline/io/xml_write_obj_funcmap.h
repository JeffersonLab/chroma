// -*- C++ -*-
// $Id: xml_write_obj_funcmap.h,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#ifndef __xml_write_obj_funcmap_h__
#define __xml_write_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Write object function map
  /*! \ingroup inlineio */
  namespace XMLWriteObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Write object function map
    /*! \ingroup inlineio */
    typedef SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void,
		  std::string,
		  TYPELIST_2(const string&,
			     const string&),
		  void (*)(const string& buffer_id,
			   const string& filename),
		  StringFunctionMapError> >
    TheXMLWriteObjFuncMap;

    extern bool registered;   // forward decl
  }

} // end namespace Chroma


#endif
