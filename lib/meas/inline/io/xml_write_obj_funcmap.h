// -*- C++ -*-
// $Id: xml_write_obj_funcmap.h,v 2.2 2006-03-24 22:16:40 edwards Exp $
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
